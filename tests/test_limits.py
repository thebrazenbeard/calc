"""Independent rational and numerical examples for bounded limit tools."""

from fractions import Fraction
import math
import random
import unittest

from calc.algebra import Polynomial
from calc.limits import (
    EpsilonDeltaCertificate,
    OneSidedSamples,
    polynomial_limit_certificate,
    one_sided_samples,
)


class PolynomialEpsilonDeltaTests(unittest.TestCase):
    def assert_witnesses(self, polynomial, certificate):
        for multiplier in (Fraction(1, 2), Fraction(3, 4), Fraction(9, 10)):
            for sign in (-1, 1):
                x = certificate.point + sign * multiplier * certificate.delta
                actual = abs(polynomial(x) - certificate.value)
                self.assertLess(actual, certificate.epsilon)
                self.assertLessEqual(actual, certificate.error_bound_at(x))
                self.assertLess(certificate.error_bound_at(x), certificate.epsilon)

    def test_affine_exact_certificate(self):
        polynomial = Polynomial((1, 2))
        result = polynomial_limit_certificate(polynomial, 2, Fraction(1, 5))
        self.assertIsInstance(result, EpsilonDeltaCertificate)
        self.assertEqual(result.value, 5)
        self.assertEqual(result.lipschitz_bound, 2)
        self.assertEqual(result.delta, Fraction(1, 20))
        self.assert_witnesses(polynomial, result)

    def test_quadratic_local_bound_is_not_global_derivative_claim(self):
        polynomial = Polynomial((0, 0, 1))
        result = polynomial_limit_certificate(
            polynomial, 1, Fraction(1, 8), radius=1
        )
        self.assertEqual(result.value, 1)
        self.assertEqual(result.lipschitz_bound, 4)
        self.assertEqual(result.delta, Fraction(1, 64))
        self.assert_witnesses(polynomial, result)

    def test_cubic_negative_center_and_rational_coefficients(self):
        polynomial = Polynomial((Fraction(-3, 7), 0, 0, Fraction(2, 3)))
        result = polynomial_limit_certificate(
            polynomial, Fraction(-1, 2), Fraction(1, 7), radius=Fraction(1, 2)
        )
        self.assertEqual(result.lipschitz_bound, 2)
        self.assertEqual(result.value, Fraction(-43, 84))
        self.assertEqual(result.delta, Fraction(1, 28))
        self.assert_witnesses(polynomial, result)

    def test_constant_and_zero_polynomial(self):
        for polynomial in (Polynomial((Fraction(5, 7),)), Polynomial(())):
            result = polynomial_limit_certificate(
                polynomial, Fraction(7, 2), Fraction(1, 100), radius=Fraction(13, 4)
            )
            self.assertEqual(result.delta, Fraction(13, 4))
            self.assertEqual(result.lipschitz_bound, 0)
            self.assert_witnesses(polynomial, result)

    def test_sequence_input_and_immutable_record(self):
        result = polynomial_limit_certificate((Fraction(1, 2), 1), 0, 1)
        self.assertEqual(result.value, Fraction(1, 2))
        self.assertEqual(result.delta, Fraction(1, 2))
        with self.assertRaises((AttributeError, TypeError)):
            result.delta = Fraction(1, 3)

    def test_outside_certificate_is_not_accepted_as_bound(self):
        result = polynomial_limit_certificate((0, 1), 2, 1)
        for x in (result.point + result.delta, result.point - result.delta):
            with self.assertRaises(ValueError):
                result.error_bound_at(x)

    def test_invalid_exact_and_resource_inputs(self):
        cases = [
            ((1, 2), 0.5, 1),
            ((1, 2), 0, 0),
            ((1, 2), 0, -1),
            ((1, True), 0, 1),
            ((1, 2), 0, 1.0),
            ((1, 2), 0, 1, 0),
            ((1, 2), 0, 1, -1),
            ((1 << 9000,), 0, 1),
            (tuple([1] * 35), 0, 1),
        ]
        for parts in cases:
            with self.subTest(parts=str(parts)[:90]):
                with self.assertRaises((TypeError, ValueError)):
                    polynomial_limit_certificate(*parts[:3], radius=parts[3] if len(parts) == 4 else 1)

    def test_reject_unbounded_iterable_without_exhausting_it(self):
        def endless():
            while True:
                yield 1

        with self.assertRaises(ValueError):
            polynomial_limit_certificate(endless(), 0, 1)

    def test_seeded_polynomial_certificates_against_separate_horner_evaluation(self):
        rng = random.Random(20261101)
        for _ in range(120):
            degree = rng.randrange(0, 7)
            coefficients = tuple(Fraction(rng.randrange(-9, 10), rng.randrange(1, 8))
                                 for _ in range(degree + 1))
            polynomial = Polynomial(coefficients)
            point = Fraction(rng.randrange(-5, 6), rng.randrange(1, 7))
            epsilon = Fraction(rng.randrange(1, 5), rng.randrange(5, 20))
            radius = Fraction(rng.randrange(1, 5), rng.randrange(1, 6))
            certificate = polynomial_limit_certificate(
                coefficients, point, epsilon, radius=radius
            )
            self.assertEqual(certificate.value, polynomial(point))
            self.assertGreater(certificate.delta, 0)
            self.assertLessEqual(certificate.delta, radius)
            self.assert_witnesses(polynomial, certificate)


class OneSidedDiagnosticTests(unittest.TestCase):
    def test_jump_values_are_visible_without_a_false_limit_verdict(self):
        result = one_sided_samples(lambda x: -1 if x < 0 else 2, 0, initial_step=0.5, count=6)
        self.assertIsInstance(result, OneSidedSamples)
        self.assertEqual(result.left, (-1.0,) * 6)
        self.assertEqual(result.right, (2.0,) * 6)
        self.assertEqual(result.offsets, (0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625))
        self.assertFalse(hasattr(result, "converges"))

    def test_removable_point_need_not_be_evaluated(self):
        def removable(x):
            if x == 1.0:
                raise AssertionError("center must never be evaluated")
            return (x * x - 1) / (x - 1)

        result = one_sided_samples(removable, 1, initial_step=0.25, count=8)
        for offset, left, right in zip(result.offsets, result.left, result.right):
            self.assertAlmostEqual(left, 2 - offset, places=12)
            self.assertAlmostEqual(right, 2 + offset, places=12)

    def test_continuous_polynomial_samples_are_only_samples(self):
        result = one_sided_samples(lambda x: 3 * x + 4, -2, initial_step=0.1, count=5)
        self.assertTrue(all(math.isfinite(v) for v in result.left + result.right))
        self.assertAlmostEqual(result.left[-1], -2 - 3 * result.offsets[-1])
        self.assertAlmostEqual(result.right[-1], -2 + 3 * result.offsets[-1])

    def test_reject_invalid_callbacks_steps_counts_and_roundoff(self):
        invalids = [
            (None, 0, 1, 3),
            (lambda x: x, float("nan"), 1, 3),
            (lambda x: x, 0, 0, 3),
            (lambda x: x, 0, float("inf"), 3),
            (lambda x: x, 0, 1, 0),
            (lambda x: x, 0, 1, 65),
            (lambda x: x, 0, 1, True),
            (lambda x: x, 1, 1e-300, 3),
            (lambda x: x, 1, 1.0, 64),
            (lambda x: float("nan"), 0, 1, 3),
            (lambda x: complex(x), 0, 1, 3),
        ]
        for function, point, step, count in invalids:
            with self.subTest(point=point, step=step, count=count):
                with self.assertRaises((TypeError, ValueError)):
                    one_sided_samples(function, point, initial_step=step, count=count)


if __name__ == "__main__":
    unittest.main()
