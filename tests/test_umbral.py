"""Independent identities for bounded exact classical finite operator calculus."""

from fractions import Fraction
import unittest

from calc.umbral import (
    falling_factorial, rising_factorial, stirling_first, stirling_second,
    forward_difference, monomial_to_falling, falling_to_monomial,
    newton_forward_interpolate,
)


def evaluate(coefficients, x):
    return sum((a * Fraction(x) ** i for i, a in enumerate(coefficients)), Fraction(0))


class FactorialBasisTests(unittest.TestCase):
    def test_integer_and_fraction_factorials(self):
        self.assertEqual(falling_factorial(5, 3), 60)
        self.assertEqual(rising_factorial(3, 4), 360)
        self.assertEqual(falling_factorial(Fraction(-1, 2), 3), Fraction(-15, 8))
        self.assertEqual(falling_factorial(99, 0), 1)
        self.assertEqual(rising_factorial(0, 4), 0)

    def test_stirling_tables_independent_classic_rows(self):
        self.assertEqual(tuple(stirling_first(4, k) for k in range(5)), (0, -6, 11, -6, 1))
        self.assertEqual(tuple(stirling_second(4, k) for k in range(5)), (0, 1, 7, 6, 1))
        self.assertEqual(stirling_first(0, 0), 1)
        self.assertEqual(stirling_second(0, 0), 1)
        self.assertEqual(stirling_first(5, 6), 0)

    def test_falling_as_signed_stirling_expansion(self):
        for degree in range(9):
            representation = falling_to_monomial(tuple([0] * degree + [1]))
            for x in (Fraction(-3, 2), Fraction(0), Fraction(7, 3)):
                self.assertEqual(evaluate(representation, x), falling_factorial(x, degree))


class DifferenceOperatorTests(unittest.TestCase):
    def test_forward_difference_cubic(self):
        cubic = (0, 0, 0, 1)
        self.assertEqual(forward_difference(cubic), (Fraction(1), Fraction(3), Fraction(3)))
        self.assertEqual(forward_difference(cubic, order=2), (Fraction(6), Fraction(6)))
        self.assertEqual(forward_difference(cubic, order=3), (Fraction(6),))
        self.assertEqual(forward_difference(cubic, order=4), (Fraction(0),))

    def test_step_and_analytic_quotient_identity(self):
        p = (Fraction(2, 3), -5, Fraction(1, 7), 4)
        for h in (Fraction(1, 2), Fraction(3, 2), -2):
            delta = forward_difference(p, step=h)
            for x in (Fraction(-2, 3), 0, Fraction(7, 5)):
                self.assertEqual(evaluate(delta, x), evaluate(p, x + h) - evaluate(p, x))

    def test_delta_lowers_falling_factorial_basis(self):
        for n in range(1, 11):
            f = falling_to_monomial(tuple([0] * n + [1]))
            d = forward_difference(f)
            smaller = falling_to_monomial(tuple([0] * (n - 1) + [n]))
            self.assertEqual(d, smaller)

    def test_zero_and_order_zero(self):
        self.assertEqual(forward_difference([]), (Fraction(0),))
        self.assertEqual(forward_difference([0, 0, 3, 0], order=0), (Fraction(0), Fraction(0), Fraction(3)))


class BasisTransformTests(unittest.TestCase):
    def test_classic_monomial_fourth_power(self):
        self.assertEqual(monomial_to_falling((0, 0, 0, 0, 1)),
                         (Fraction(0), Fraction(1), Fraction(7), Fraction(6), Fraction(1)))

    def test_inversion_over_many_polynomials(self):
        for degree in range(12):
            monomial = tuple(Fraction((-1) ** k * (k + 2), k + 1)
                             for k in range(degree + 1))
            falling = monomial_to_falling(monomial)
            self.assertEqual(falling_to_monomial(falling), monomial)

    def test_falling_monomial_roundtrip(self):
        for k in range(12):
            coefficients = tuple(Fraction(j * j - 5, 3) for j in range(k + 1))
            self.assertEqual(monomial_to_falling(falling_to_monomial(coefficients)),
                             coefficients)


class NewtonInterpolationTests(unittest.TestCase):
    def test_quadratic_on_offset_integer_grid(self):
        values = [(-2 + 3 * j) ** 2 + 4 * (-2 + 3 * j) - 2 for j in range(4)]
        self.assertEqual(newton_forward_interpolate(values, start=-2, step=3),
                         (Fraction(-2), Fraction(4), Fraction(1)))

    def test_cubic_and_fractional_uniform_grid(self):
        p = (Fraction(1, 3), Fraction(-2, 5), 0, Fraction(3, 4))
        start, h = Fraction(-3, 4), Fraction(2, 3)
        data = [evaluate(p, start + k * h) for k in range(4)]
        self.assertEqual(newton_forward_interpolate(data, start=start, step=h), p)

    def test_singleton_and_empty(self):
        self.assertEqual(newton_forward_interpolate([Fraction(4, 7)]), (Fraction(4, 7),))
        with self.assertRaises(ValueError):
            newton_forward_interpolate([])

    def test_interpolation_respects_data(self):
        values = (3, 8, -1, 20, 100)
        coefficients = newton_forward_interpolate(values, start=4, step=2)
        for index, y in enumerate(values):
            self.assertEqual(evaluate(coefficients, 4 + 2 * index), y)


class StrictValidationTests(unittest.TestCase):
    def test_no_inexact_or_bool_inputs(self):
        for func, args in ((falling_factorial, (0.1, 2)),
                           (rising_factorial, (True, 2)),
                           (forward_difference, ((1.0, 2),)),
                           (newton_forward_interpolate, ((0, 1),))):
            if func is newton_forward_interpolate:
                with self.assertRaises(TypeError):
                    func(*args, step=0.5)
            else:
                with self.assertRaises(TypeError):
                    func(*args)

    def test_orders_bounds_and_step(self):
        for bad in (-1, True, 33, Fraction(1, 2)):
            with self.assertRaises((TypeError, ValueError)):
                stirling_second(bad, 0)
        with self.assertRaises(ValueError):
            forward_difference((1,), step=0)
        with self.assertRaises(ValueError):
            newton_forward_interpolate([1, 2], step=0)

    def test_input_degree_limit(self):
        with self.assertRaises(ValueError):
            monomial_to_falling([1] * 34)
        with self.assertRaises(ValueError):
            newton_forward_interpolate([1] * 35)


if __name__ == "__main__":
    unittest.main()
