"""Analytic rational-function leading-term and order-symbol witnesses."""

from fractions import Fraction
from itertools import repeat
import random
import unittest

from calc.algebra import Polynomial
from calc.asymptotics import RationalAsymptotic, compare_growth


class RationalAsymptoticTests(unittest.TestCase):
    def test_leading_at_positive_infinity(self):
        f = RationalAsymptotic((1, 0, 3), (2, 1))
        self.assertEqual(f.leading_term(), (1, Fraction(3)))
        self.assertEqual(f.leading_term(at="zero_plus"), (0, Fraction(1, 2)))

    def test_leading_at_zero_uses_lowest_nonzero_order(self):
        f = RationalAsymptotic((0, 2, 3), (0, 0, 4))
        self.assertEqual(f.leading_term(at="zero_plus"), (-1, Fraction(1, 2)))
        self.assertEqual(f.leading_term(), (0, Fraction(3, 4)))

    def test_order_comparisons_at_infinity(self):
        x = RationalAsymptotic((0, 1))
        x2 = RationalAsymptotic((0, 0, 1))
        higher = compare_growth(x2, x)
        self.assertEqual(higher.ratio_power, 1)
        self.assertFalse(higher.big_o)
        self.assertFalse(higher.little_o)
        self.assertFalse(higher.theta)
        self.assertIsNone(higher.ratio_limit)
        lower = compare_growth(x, x2)
        self.assertTrue(lower.big_o)
        self.assertTrue(lower.little_o)
        self.assertFalse(lower.theta)
        self.assertEqual(lower.ratio_limit, 0)

    def test_order_comparisons_near_zero_reverse(self):
        x = RationalAsymptotic((0, 1))
        x2 = RationalAsymptotic((0, 0, 1))
        smaller = compare_growth(x2, x, at="zero_plus")
        self.assertTrue(smaller.big_o)
        self.assertTrue(smaller.little_o)
        self.assertEqual(smaller.ratio_power, 1)
        self.assertEqual(smaller.ratio_limit, 0)
        larger = compare_growth(x, x2, at="zero_plus")
        self.assertFalse(larger.big_o)
        self.assertIsNone(larger.ratio_limit)

    def test_asymptotic_equality_is_ratio_limit_one(self):
        first = RationalAsymptotic((1, -1), (1, 1))
        second = RationalAsymptotic((-1,))
        result = compare_growth(first, second)
        self.assertTrue(result.big_o and result.theta)
        self.assertTrue(result.asymptotically_equivalent)
        self.assertFalse(result.little_o)
        self.assertEqual(result.ratio_limit, Fraction(1))

    def test_theta_can_have_nonunit_or_negative_ratio(self):
        first = RationalAsymptotic((0, -3, 2))
        second = RationalAsymptotic((0, 0, -4))
        result = compare_growth(first, second)
        self.assertTrue(result.theta)
        self.assertEqual(result.ratio_limit, Fraction(-1, 2))
        self.assertFalse(result.asymptotically_equivalent)

    def test_exact_fraction_leading_coefficient(self):
        first = RationalAsymptotic((Fraction(3, 7), 0, Fraction(5, 8)), (Fraction(4, 9),))
        self.assertEqual(first.leading_term(), (2, Fraction(45, 32)))
        self.assertEqual(first.leading_term(at="zero_plus"), (0, Fraction(27, 28)))

    def test_common_factor_does_not_need_cancellation(self):
        numerator = Polynomial((-1, 0, 1))
        denominator = Polynomial((-1, 1))
        first = RationalAsymptotic(numerator, denominator)
        second = RationalAsymptotic((1, 1))
        self.assertTrue(compare_growth(first, second).asymptotically_equivalent)
        self.assertEqual(first.leading_term(at="zero_plus"), (0, Fraction(1)))

    def test_zero_numerator_and_zero_comparator_contract(self):
        zero = RationalAsymptotic(())
        constant = RationalAsymptotic((1,))
        outcome = compare_growth(zero, constant)
        self.assertTrue(outcome.big_o and outcome.little_o)
        self.assertFalse(outcome.theta or outcome.asymptotically_equivalent)
        self.assertEqual(outcome.ratio_limit, 0)
        with self.assertRaisesRegex(ValueError, "leading term"):
            zero.leading_term()
        with self.assertRaisesRegex(ValueError, "reference"):
            compare_growth(constant, zero)
        with self.assertRaisesRegex(ValueError, "reference"):
            compare_growth(zero, zero)

    def test_denominator_poles_do_not_impair_eventual_leading_term(self):
        first = RationalAsymptotic((1, 0, 1), (1, -1))
        self.assertEqual(first.leading_term(), (1, Fraction(-1)))
        self.assertEqual(first.leading_term(at="zero_plus"), (0, Fraction(1)))

    def test_invalid_inputs_and_direction(self):
        with self.assertRaises(ZeroDivisionError):
            RationalAsymptotic((1,), ())
        with self.assertRaises(TypeError):
            RationalAsymptotic((1.0,))
        with self.assertRaises(TypeError):
            RationalAsymptotic((True,))
        with self.assertRaises(TypeError):
            RationalAsymptotic("1")
        with self.assertRaises(ValueError):
            RationalAsymptotic((1,)).leading_term(at="sideways")
        with self.assertRaises(ValueError):
            compare_growth(RationalAsymptotic((1,)), RationalAsymptotic((1,)), at="sideways")
        with self.assertRaises(TypeError):
            compare_growth((1,), RationalAsymptotic((1,)))

    def test_resource_bounds(self):
        with self.assertRaises(ValueError):
            RationalAsymptotic(repeat(1))
        with self.assertRaises(ValueError):
            RationalAsymptotic((1 << 5000,))
        with self.assertRaises(ValueError):
            RationalAsymptotic((Fraction(1, 1 << 5000),))
        self.assertEqual(RationalAsymptotic((0,) * 64 + (1,)).leading_term(), (64, 1))

    def test_seeded_analytic_degree_oracles(self):
        rng = random.Random(519205)
        for _ in range(200):
            pa = [rng.randint(-5, 5) for _ in range(rng.randint(1, 9))]
            qa = [rng.randint(-5, 5) for _ in range(rng.randint(1, 8))]
            pb = [rng.randint(-5, 5) for _ in range(rng.randint(1, 8))]
            qb = [rng.randint(-5, 5) for _ in range(rng.randint(1, 8))]
            # Explicit independent nonzero endpoint witnesses.
            for row in (pa, qa, pb, qb):
                if not any(row):
                    row[-1] = 1
            f = RationalAsymptotic(pa, qa)
            g = RationalAsymptotic(pb, qb)
            for approach in ("infinity", "zero_plus"):
                if approach == "infinity":
                    indices = [max(i for i, v in enumerate(row) if v) for row in (pa, qa, pb, qb)]
                else:
                    indices = [min(i for i, v in enumerate(row) if v) for row in (pa, qa, pb, qb)]
                exponent = indices[0] - indices[1] - indices[2] + indices[3]
                leading = (Fraction(pa[indices[0]], qa[indices[1]]) /
                           Fraction(pb[indices[2]], qb[indices[3]]))
                outcome = compare_growth(f, g, at=approach)
                with self.subTest(approach=approach, pa=pa, qa=qa, pb=pb, qb=qb):
                    self.assertEqual(outcome.ratio_power, exponent)
                    self.assertEqual(outcome.ratio_leading_coefficient, leading)
                    decays = exponent < 0 if approach == "infinity" else exponent > 0
                    self.assertEqual(outcome.little_o, decays)
                    self.assertEqual(outcome.big_o, decays or exponent == 0)
                    self.assertEqual(outcome.theta, exponent == 0)


if __name__ == "__main__":
    unittest.main()
