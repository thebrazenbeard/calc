"""Exact rational witnesses for outward-rounded binary64 interval arithmetic."""

import math
import random
import sys
import unittest
from fractions import Fraction

from calc.intervals import Interval


def exact(value: float) -> Fraction:
    return Fraction.from_float(value)


class IntervalTests(unittest.TestCase):
    def assert_encloses(self, actual: Interval, low: Fraction, high: Fraction) -> None:
        if not math.isinf(actual.lower):
            self.assertLessEqual(exact(actual.lower), low)
        else:
            self.assertLess(actual.lower, 0)
        if not math.isinf(actual.upper):
            self.assertGreaterEqual(exact(actual.upper), high)
        else:
            self.assertGreater(actual.upper, 0)

    def test_singleton_and_exact_integer_inputs(self):
        self.assertEqual(Interval.point(2), Interval(2.0, 2.0))
        self.assertEqual(Interval(-2, 3), Interval(-2.0, 3.0))
        self.assertTrue(Interval(2, 3).contains(2.5))
        self.assertFalse(Interval(2, 3).contains(4))

    def test_addition_and_subtraction_exact_enclosure(self):
        a, b = Interval(0.1, 0.1), Interval(0.2, 0.2)
        rational_sum = exact(0.1) + exact(0.2)
        self.assert_encloses(a + b, rational_sum, rational_sum)
        rational_difference = exact(0.1) - exact(0.2)
        self.assert_encloses(a - b, rational_difference, rational_difference)
        self.assertLess((a + b).lower, (a + b).upper)

    def test_multiplication_extrema_and_signs(self):
        a, b = Interval(-2.5, 3.75), Interval(-4, -0.125)
        combinations = [exact(x) * exact(y) for x in (-2.5, 3.75) for y in (-4, -0.125)]
        self.assert_encloses(a * b, min(combinations), max(combinations))
        self.assertEqual(Interval(-1, 2) * Interval(0, 0), Interval(-0.0, 0.0))

    def test_division_extrema_and_zero(self):
        a, b = Interval(-2, 3), Interval(0.1, 0.2)
        values = [exact(x) / exact(y) for x in (-2.0, 3.0) for y in (0.1, 0.2)]
        self.assert_encloses(a / b, min(values), max(values))
        with self.assertRaises(ZeroDivisionError):
            _ = Interval(1, 2) / Interval(-1, 1)
        with self.assertRaises(ZeroDivisionError):
            _ = Interval(1, 2) / Interval(0, 1)

    def test_overflow_is_unbounded_enclosure_not_false_finite(self):
        result = Interval.point(1e308) * Interval.point(2.0)
        self.assertEqual(result.lower, sys.float_info.max)
        self.assertEqual(result.upper, math.inf)
        self.assert_encloses(result, exact(1e308) * 2, exact(1e308) * 2)
        negative = -result
        self.assertEqual(negative.lower, -math.inf)
        self.assertEqual(negative.upper, -sys.float_info.max)
        with self.assertRaises(ValueError):
            _ = result + Interval.point(1)

    def test_underflow_encloses_subnormal_exact_product(self):
        sub = math.ulp(0.0)
        positive = Interval.point(sub) * Interval.point(sub)
        self.assertEqual(positive.lower, 0.0)
        self.assertEqual(positive.upper, sub)
        self.assert_encloses(positive, exact(sub) ** 2, exact(sub) ** 2)
        negative = Interval.point(-sub) * Interval.point(sub)
        self.assertEqual(negative.lower, -sub)
        self.assertEqual(negative.upper, 0.0)
        self.assert_encloses(negative, -(exact(sub) ** 2), -(exact(sub) ** 2))

    def test_signed_zero_directed_bounds(self):
        zero = Interval.point(0.0) + Interval.point(-0.0)
        self.assertEqual(math.copysign(1, zero.lower), -1)
        self.assertEqual(math.copysign(1, zero.upper), 1)

    def test_intersection_and_hull(self):
        self.assertEqual(Interval(1, 5).intersection(Interval(3, 8)), Interval(3, 5))
        self.assertIsNone(Interval(1, 2).intersection(Interval(3, 4)))
        self.assertEqual(Interval(1, 2).hull(Interval(-3, 1)), Interval(-3, 2))
        self.assertEqual(Interval(1, 2).intersection(Interval(2, 3)), Interval(2, 2))

    def test_bad_constructor_and_scalar_inputs(self):
        with self.assertRaises(ValueError):
            Interval(2, 1)
        with self.assertRaises(ValueError):
            Interval(float("nan"), 1)
        with self.assertRaises(TypeError):
            Interval(True, 1)
        with self.assertRaises(ValueError):
            Interval(2**53 + 1, 2**53 + 2)
        with self.assertRaises(TypeError):
            Interval.point(Fraction(1, 3))

    def test_randomized_binary64_rational_enclosure(self):
        rng = random.Random(42017)
        for _ in range(150):
            a, b, c, d = [rng.uniform(-100, 100) for _ in range(4)]
            x, y = Interval(min(a, b), max(a, b)), Interval(min(c, d), max(c, d))
            xl, xh = exact(x.lower), exact(x.upper)
            yl, yh = exact(y.lower), exact(y.upper)
            self.assert_encloses(x + y, xl + yl, xh + yh)
            self.assert_encloses(x - y, xl - yh, xh - yl)
            candidates = (xl * yl, xl * yh, xh * yl, xh * yh)
            self.assert_encloses(x * y, min(candidates), max(candidates))
            if not y.contains(0):
                quotients = (xl / yl, xl / yh, xh / yl, xh / yh)
                self.assert_encloses(x / y, min(quotients), max(quotients))


if __name__ == "__main__":
    unittest.main()

