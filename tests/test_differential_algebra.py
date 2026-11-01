"""Independent exact witnesses for the truncated differential polynomial ring."""

from fractions import Fraction
from itertools import product
import random
import unittest

from calc.differential_algebra import JetPolynomial


class JetPolynomialTests(unittest.TestCase):
    def test_constructor_canonicalizes_to_immutable_rational_terms(self):
        original = {(1, 0, 0): 3, (0, 1, 0): 0}
        polynomial = JetPolynomial(2, original)
        original[(1, 0, 0)] = 19
        self.assertEqual(polynomial.as_dict(), {(1, 0, 0): Fraction(3)})
        self.assertIsInstance(polynomial.coefficients, tuple)
        self.assertEqual(JetPolynomial(2).coefficients, ())

    def test_jet_generators_are_exact_and_ordered(self):
        for index in range(5):
            symbol = JetPolynomial.jet(index, max_jet=4)
            self.assertEqual(symbol.as_dict(), {
                tuple(int(j == index) for j in range(5)): Fraction(1)
            })
        self.assertEqual(JetPolynomial.constant(Fraction(2, 3), max_jet=4)
                         .evaluate((0, 0, 0, 0, 0)), Fraction(2, 3))

    def test_total_derivative_advances_each_represented_jet(self):
        for index in range(4):
            symbol = JetPolynomial.jet(index, max_jet=4)
            self.assertEqual(symbol.total_derivative(),
                             JetPolynomial.jet(index + 1, max_jet=4))
        with self.assertRaisesRegex(ValueError, "higher unrepresented jet"):
            JetPolynomial.jet(4, max_jet=4).total_derivative()

    def test_total_derivative_of_constants_and_polynomial_power(self):
        u = JetPolynomial.jet(0, max_jet=4)
        v = JetPolynomial.jet(1, max_jet=4)
        self.assertEqual(JetPolynomial.constant(7, max_jet=4)
                         .total_derivative(), JetPolynomial(4))
        self.assertEqual((u ** 3).total_derivative(), 3 * u ** 2 * v)
        self.assertEqual(((u + 2) ** 4).total_derivative(),
                         4 * (u + 2) ** 3 * v)

    def test_total_derivative_of_mixed_monomials(self):
        u, v, w, z = (JetPolynomial.jet(i, max_jet=4) for i in range(4))
        self.assertEqual((u * v).total_derivative(), v * v + u * w)
        self.assertEqual((v * w).total_derivative(), w * w + v * z)
        self.assertEqual((u ** 2 * v ** 3).total_derivative(),
                         2 * u * v ** 4 + 3 * u ** 2 * v ** 2 * w)

    def test_leibniz_linearity_and_commutativity_exact(self):
        u, v, w = (JetPolynomial.jet(i, max_jet=5) for i in range(3))
        a = Fraction(2, 7) * (u ** 2 + 3 * v) - w
        b = Fraction(-5, 3) * (u * v + 1)
        self.assertEqual((a * b).total_derivative(),
                         a.total_derivative() * b + a * b.total_derivative())
        self.assertEqual((a + b).total_derivative(),
                         a.total_derivative() + b.total_derivative())
        self.assertEqual(a * b, b * a)
        self.assertEqual(a * (b + u), a * b + a * u)

    def test_seeded_exact_leibniz_witnesses(self):
        random_gen = random.Random(2701)
        u, v = (JetPolynomial.jet(i, max_jet=4) for i in range(2))
        for _ in range(70):
            a, b, c, d, e, f = (Fraction(random_gen.randint(-4, 4),
                                           random_gen.randint(1, 5))
                                for _ in range(6))
            left = a * u * u + b * u * v + c
            right = d * u + e * v ** 2 + f
            direct = (left * right).total_derivative()
            oracle = left.total_derivative() * right + left * right.total_derivative()
            self.assertEqual(direct, oracle)

    def test_formal_partial_derivative_and_total_derivative_distinct(self):
        u = JetPolynomial.jet(0, max_jet=3)
        v = JetPolynomial.jet(1, max_jet=3)
        f = 3 * u ** 2 * v + 7 * v
        self.assertEqual(f.partial_derivative(0), 6 * u * v)
        self.assertEqual(f.partial_derivative(1), 3 * u ** 2 + 7)
        self.assertEqual(f.partial_derivative(3), JetPolynomial(3))
        self.assertNotEqual(f.partial_derivative(0), f.total_derivative())

    def test_repeated_total_derivative_of_base_jet(self):
        result = JetPolynomial.jet(0, max_jet=4)
        for index in range(1, 5):
            result = result.total_derivative()
            self.assertEqual(result, JetPolynomial.jet(index, max_jet=4))

    def test_exact_evaluation_against_independent_substitution(self):
        u, v, w, z = (JetPolynomial.jet(i, max_jet=3) for i in range(4))
        polynomial = 3 * u ** 2 * v + 2 * w - Fraction(1, 3) * z
        points = (Fraction(2), Fraction(3, 2), Fraction(4), Fraction(9))
        expected = 3 * 4 * Fraction(3, 2) + 2 * 4 - Fraction(9, 3)
        self.assertEqual(polynomial.evaluate(points), expected)

    def test_infinite_iterable_evaluation_is_bounded(self):
        u = JetPolynomial.jet(0, max_jet=2)
        with self.assertRaises(ValueError):
            u.evaluate((1 for _ in range(100000000)))
        with self.assertRaises(ValueError):
            u.evaluate((1, 2))

    def test_scalar_multiplication_division_and_zero(self):
        u = JetPolynomial.jet(0, max_jet=2)
        self.assertEqual((u * Fraction(3, 4)) / Fraction(3, 4), u)
        self.assertEqual(u - u, JetPolynomial(2))
        self.assertEqual(1 - u + u, JetPolynomial.constant(1, max_jet=2))
        self.assertEqual(u ** 0, JetPolynomial.constant(1, max_jet=2))
        with self.assertRaises(ZeroDivisionError):
            _ = u / 0

    def test_jet_order_and_exponent_validation(self):
        for invalid in (-1, 9, True, 1.5):
            with self.subTest(max_jet=invalid), self.assertRaises((TypeError, ValueError)):
                JetPolynomial(invalid)
        for terms in ({(1, 0): 2}, {(1, 0, 0, 0): 1},
                      {(True, 0, 0): 1}, {(-1, 0, 0): 1}):
            with self.subTest(terms=terms), self.assertRaises((TypeError, ValueError)):
                JetPolynomial(2, terms)
        with self.assertRaises(ValueError):
            JetPolynomial.jet(4, max_jet=3)

    def test_exact_field_and_ambient_ring_validation(self):
        with self.assertRaises(TypeError):
            JetPolynomial(2, {(1, 0, 0): 1.0})
        with self.assertRaises(TypeError):
            JetPolynomial(2, {(1, 0, 0): True})
        with self.assertRaises(ValueError):
            _ = JetPolynomial.jet(0, max_jet=2) + JetPolynomial.jet(0, max_jet=3)
        with self.assertRaises(TypeError):
            _ = JetPolynomial.jet(0, max_jet=2) + 1.5

    def test_degree_and_rational_bit_limits(self):
        with self.assertRaises(ValueError):
            JetPolynomial(1, {(17, 0): 1})
        with self.assertRaises(ValueError):
            JetPolynomial(1, {(0, 0): 1 << 4096})
        u = JetPolynomial.jet(0, max_jet=1)
        with self.assertRaises(ValueError):
            _ = u ** 17
        with self.assertRaises(ValueError):
            _ = (u ** 16) * u
        with self.assertRaises(TypeError):
            _ = u ** 1.5

    def test_monomial_count_and_multiplication_work_limits(self):
        points = list(product(range(2), repeat=9))
        p = JetPolynomial(8, {entry: 1 for entry in points[:142]})
        with self.assertRaisesRegex(ValueError, "work limit"):
            _ = p * p
        with self.assertRaises(ValueError):
            JetPolynomial(8, {entry: 1 for entry in [*points, *product(range(3), repeat=7)]})

    def test_truncation_does_not_silently_drop_highest_jet(self):
        top = JetPolynomial.jet(2, max_jet=2)
        lower = JetPolynomial.jet(0, max_jet=2)
        with self.assertRaisesRegex(ValueError, "higher unrepresented"):
            (top * lower).total_derivative()
        self.assertEqual((lower * lower).total_derivative(),
                         2 * lower * JetPolynomial.jet(1, max_jet=2))


if __name__ == "__main__":
    unittest.main()
