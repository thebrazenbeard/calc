"""Independent exact finite-measure-theory identity and boundary tests."""

from fractions import Fraction as F
import unittest

from calc.measure import FiniteMeasure


class FiniteMeasureTests(unittest.TestCase):
    def test_canonical_immutable_copy_and_total(self):
        supplied = {"z": F(1, 6), "a": 2, "b": 0}
        measure = FiniteMeasure(supplied)
        supplied["a"] = 100
        self.assertEqual(measure.atoms, ("a", "b", "z"))
        self.assertEqual(measure.weights, (F(2), F(0), F(1, 6)))
        self.assertEqual(measure.total_mass, F(13, 6))
        self.assertEqual(measure.support, ("a", "z"))
        with self.assertRaises(AttributeError):
            measure.atoms = ("q",)

    def test_event_mass_is_setwise_additive(self):
        measure = FiniteMeasure({"a": F(2, 3), "b": F(1, 9), "c": 0})
        self.assertEqual(measure.mass(("a", "a")), F(2, 3))
        self.assertEqual(measure.mass(("a", "b")), F(7, 9))
        self.assertEqual(measure.mass(()), F(0))
        self.assertEqual(measure.mass(measure.atoms), measure.total_mass)

    def test_linearity_of_integration(self):
        measure = FiniteMeasure({"a": F(2, 3), "b": F(1, 6), "c": F(4, 5)})
        f = {"a": F(3, 5), "b": F(-7, 2), "c": F(1, 3)}
        g = {"a": F(4, 7), "b": F(2, 5), "c": -1}
        left = measure.integrate({key: 2 * f[key] - 3 * g[key] for key in f})
        self.assertEqual(left, 2 * measure.integrate(f) - 3 * measure.integrate(g))
        self.assertEqual(measure.integrate({key: 1 for key in f}), measure.total_mass)

    def test_monotonicity_and_nonnegative_squares(self):
        measure = FiniteMeasure({"a": 1, "b": F(3, 7), "c": 0})
        f = {"a": -3, "b": 4, "c": -99}
        positive = {key: f[key] + 5 if key != "c" else 0 for key in f}
        self.assertGreaterEqual(measure.integrate(positive), measure.integrate(f))
        self.assertGreaterEqual(measure.integrate({x: f[x] ** 2 for x in f}), 0)

    def test_probability_normalization_and_variance(self):
        mu = FiniteMeasure({"a": 1, "b": 3, "c": 0})
        p = mu.normalize()
        self.assertEqual(p.total_mass, F(1))
        self.assertEqual(p.weights, (F(1, 4), F(3, 4), F(0)))
        self.assertEqual(p.expectation({"a": 1, "b": 5, "c": 1000}), F(4))
        self.assertEqual(p.variance({"a": 1, "b": 5, "c": 1000}), F(3))
        self.assertEqual(mu.expectation({"a": 1, "b": 5, "c": 1000}), F(4))

    def test_variance_identity_is_exact(self):
        mu = FiniteMeasure({"x": F(2, 7), "y": F(3, 11), "z": F(4, 13)})
        f = {"x": F(-5, 3), "y": F(7, 10), "z": F(11, 2)}
        mean = mu.expectation(f)
        second = mu.expectation({x: v * v for x, v in f.items()})
        self.assertEqual(mu.variance(f), second - mean * mean)

    def test_normalization_preserves_valid_large_exact_output(self):
        # Combining valid bounded input denominators can require >2048 bits
        # for normalized probabilities; these are derived results, not a new
        # untrusted input exceeding the original input-size ceiling.
        mu = FiniteMeasure({"a": F(1, 2**2040 - 1),
                            "b": F(1, 2**2040 - 3),
                            "c": F(1, 2**2040 - 7)})
        p = mu.normalize()
        self.assertEqual(p.total_mass, F(1))
        self.assertEqual(sum(p.weights, F(0)), F(1))
        self.assertGreater(max(weight.denominator.bit_length() for weight in p.weights), 2048)

    def test_radon_nikodym_density_reconstructs_target_integrals(self):
        ref = FiniteMeasure({"a": F(2, 3), "b": F(4, 5), "c": 0})
        target = FiniteMeasure({"a": F(1, 9), "b": F(3, 10), "c": 0})
        density = ref.density_of(target)
        self.assertEqual(density, {"a": F(1, 6), "b": F(3, 8), "c": F(0)})
        self.assertEqual(ref.integrate(density), target.total_mass)
        observable = {"a": F(6, 7), "b": -5, "c": 31337}
        self.assertEqual(
            ref.integrate({x: observable[x] * density[x] for x in ref.atoms}),
            target.integrate(observable),
        )

    def test_radon_nikodym_null_reference_requires_null_target(self):
        ref = FiniteMeasure({"a": 2, "b": 0})
        with self.assertRaisesRegex(ValueError, "absolutely continuous"):
            ref.density_of(FiniteMeasure({"a": 1, "b": 1}))
        self.assertEqual(ref.density_of(FiniteMeasure({"a": 1, "b": 0}))["b"], 0)

    def test_conditional_expectation_weighted_mean(self):
        mu = FiniteMeasure({"a": 1, "b": 3, "c": 2, "d": 0})
        f = {"a": -1, "b": 7, "c": 11, "d": 999}
        conditional = mu.conditional_expectation(f, (("b", "a"), ("c",), ("d",)))
        self.assertEqual(conditional, {"a": F(5), "b": F(5), "c": F(11), "d": F(0)})
        self.assertEqual(mu.integrate(conditional), mu.integrate(f))

    def test_conditional_tower_property_for_nested_partitions(self):
        mu = FiniteMeasure({"a": 2, "b": 3, "c": 5, "d": 7})
        f = {"a": F(1, 3), "b": -4, "c": F(13, 7), "d": 6}
        fine = (("a",), ("b", "c"), ("d",))
        coarse = (("a", "b", "c"), ("d",))
        fine_expectation = mu.conditional_expectation(f, fine)
        coarse_direct = mu.conditional_expectation(f, coarse)
        coarse_via_fine = mu.conditional_expectation(fine_expectation, coarse)
        self.assertEqual(coarse_direct, coarse_via_fine)

    def test_conditional_partition_integral_identity(self):
        mu = FiniteMeasure({"a": F(3, 5), "b": 0, "c": F(1, 7), "d": 1})
        f = {"a": F(11, 13), "b": 1000, "c": -2, "d": 3}
        partition = (("a", "b", "c"), ("d",))
        c = mu.conditional_expectation(f, partition)
        for cell in partition:
            indicator = {atom: c[atom] if atom in cell else 0 for atom in mu.atoms}
            original = {atom: f[atom] if atom in cell else 0 for atom in mu.atoms}
            self.assertEqual(mu.integrate(indicator), mu.integrate(original))

    def test_null_measure_cannot_normalize_or_expect(self):
        mu = FiniteMeasure({"a": 0, "b": 0})
        with self.assertRaises(ValueError):
            mu.normalize()
        with self.assertRaises(ValueError):
            mu.expectation({"a": 1, "b": 2})
        self.assertEqual(mu.conditional_expectation({"a": 1, "b": 2}, (("a", "b"),)),
                         {"a": F(0), "b": F(0)})

    def test_input_types_and_domain_validation(self):
        with self.assertRaises(TypeError):
            FiniteMeasure({1: 2})
        with self.assertRaises(TypeError):
            FiniteMeasure({"a": 1.0})
        with self.assertRaises(TypeError):
            FiniteMeasure({"a": True})
        with self.assertRaises(ValueError):
            FiniteMeasure({"a": -1})
        with self.assertRaises(ValueError):
            FiniteMeasure({})
        measure = FiniteMeasure({"a": 1, "b": 2})
        with self.assertRaises(ValueError):
            measure.integrate({"a": 1})
        with self.assertRaises(ValueError):
            measure.integrate({"a": 1, "b": 2, "c": 3})
        with self.assertRaises(ValueError):
            measure.mass(("c",))
        with self.assertRaises(ValueError):
            measure.density_of(FiniteMeasure({"x": 2}))

    def test_invalid_partitions_rejected(self):
        mu = FiniteMeasure({"a": 1, "b": 2, "c": 3})
        f = {"a": 1, "b": 2, "c": 3}
        bad = ((), (("a",),), (("a", "b"), ("b", "c")),
               (("a", "a"), ("b", "c")), (("a", "b", "z"), ("c",)),
               (("a", "b", "c"), ()), (("a", "b", "c"), ("c",)))
        for partition in bad:
            with self.subTest(partition=partition), self.assertRaises(ValueError):
                mu.conditional_expectation(f, partition)
        with self.assertRaises(TypeError):
            mu.conditional_expectation(f, ("abc",))

    def test_bounded_weights_and_finite_iterables(self):
        with self.assertRaises(ValueError):
            FiniteMeasure({f"a{i}": 1 for i in range(129)})
        with self.assertRaises(ValueError):
            FiniteMeasure({"a": 1 << 2049})
        with self.assertRaises(TypeError):
            FiniteMeasure({"a": F(1), "b": "1"})
        mu = FiniteMeasure({"a": 1})
        with self.assertRaises(ValueError):
            mu.mass(("a" for _ in range(1_000_000_000)))
        with self.assertRaises(ValueError):
            mu.conditional_expectation({"a": 1}, (("a",) for _ in range(1_000_000_000)))


if __name__ == "__main__":
    unittest.main()
