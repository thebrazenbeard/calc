"""Independent mathematical checks for select advanced-calculus operations."""

import math
import unittest

from calc.advanced import (
    euler_lagrange_residual,
    finite_difference,
    q_derivative,
    riemann_liouville_power_integral,
    wedge,
)


class FractionalCalculusTests(unittest.TestCase):
    def test_fractional_integral_of_one(self):
        # I^(1/2) 1 (x) = 2 sqrt(x) / sqrt(pi).
        self.assertAlmostEqual(
            riemann_liouville_power_integral(0.5, 0.0, 4.0),
            4.0 / math.sqrt(math.pi),
            places=12,
        )

    def test_integer_order_reduces_to_integral(self):
        # I^1 (t) (x) = x^2 / 2 with lower endpoint zero.
        self.assertAlmostEqual(riemann_liouville_power_integral(1, 1, 6), 18)

    def test_invalid_domain(self):
        for args in ((0, 0, 1), (-1, 0, 1), (0.5, -1, 1), (0.5, 0, -1)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                riemann_liouville_power_integral(*args)

    def test_divergent_zero_endpoint(self):
        with self.assertRaises(ValueError):
            riemann_liouville_power_integral(0.25, -0.5, 0)


class DiscreteCalculusTests(unittest.TestCase):
    def test_first_forward_difference(self):
        self.assertEqual(finite_difference(lambda x: x * x, 3, step=2), 8)

    def test_second_forward_difference(self):
        self.assertEqual(finite_difference(lambda x: x * x, 3, step=2, order=2), 2)

    def test_zero_order(self):
        self.assertEqual(finite_difference(lambda x: 2 * x + 1, 3, order=0), 7)

    def test_invalid_step_and_order(self):
        with self.assertRaises(ValueError):
            finite_difference(lambda x: x, 2, step=0)
        with self.assertRaises(ValueError):
            finite_difference(lambda x: x, 2, order=-1)
        with self.assertRaises(TypeError):
            finite_difference(lambda x: x, 2, order=1.5)

    def test_q_derivative_polynomial(self):
        # D_q(x^2) = (q + 1) x.
        self.assertEqual(q_derivative(lambda x: x * x, 3, q=2), 9)

    def test_q_derivative_domain(self):
        for x, q in ((3, 1), (0, 2), (3, 0)):
            with self.subTest(x=x, q=q), self.assertRaises(ValueError):
                q_derivative(lambda value: value**2, x, q=q)


class ExteriorCalculusTests(unittest.TestCase):
    def test_anticommuting_basis(self):
        self.assertEqual(wedge({(0,): 2}, {(1,): 3}), {(0, 1): 6})
        self.assertEqual(wedge({(1,): 3}, {(0,): 2}), {(0, 1): -6})

    def test_repeated_basis_is_zero(self):
        self.assertEqual(wedge({(1,): 2}, {(1,): 3}), {})

    def test_two_forms_graded_commute(self):
        self.assertEqual(wedge({(0, 1): 1}, {(2,): 5}), {(0, 1, 2): 5})
        self.assertEqual(wedge({(2,): 5}, {(0, 1): 1}), {(0, 1, 2): 5})

    def test_associativity(self):
        a, b, c = {(2,): 2}, {(0,): 3}, {(1,): 5}
        self.assertEqual(wedge(wedge(a, b), c), wedge(a, wedge(b, c)))

    def test_scalar_wedge_and_cancellation(self):
        self.assertEqual(wedge({(): 3}, {(2,): 4}), {(2,): 12})
        self.assertEqual(wedge({(0,): 1, (1,): 1}, {(0,): 1, (1,): 1}), {})
        self.assertEqual(wedge({(): 0, (0,): 2}, {(1,): 1}), {(0, 1): 2})

    def test_invalid_basis(self):
        with self.assertRaises(ValueError):
            wedge({(1, 0): 1}, {(2,): 1})
        with self.assertRaises(ValueError):
            wedge({(1, 1): 1}, {(2,): 1})
        with self.assertRaises(ValueError):
            wedge({(0,): 1, (0, 1): 1}, {(2,): 1})


class VariationalCalculusTests(unittest.TestCase):
    def test_harmonic_oscillator_stationary_path(self):
        lagrangian = lambda t, q, v: 0.5 * (v * v - q * q)
        self.assertAlmostEqual(
            euler_lagrange_residual(lagrangian, math.cos, 0.3), 0.0, delta=1e-4
        )

    def test_nonstationary_path(self):
        lagrangian = lambda t, q, v: 0.5 * (v * v - q * q)
        # q = 1: d/dt dL/dv - dL/dq = +1, not a stationary path.
        self.assertAlmostEqual(
            euler_lagrange_residual(lagrangian, lambda t: 1, 0.3), 1.0,
            delta=1e-4,
        )

    def test_invalid_step(self):
        with self.assertRaises(ValueError):
            euler_lagrange_residual(lambda t, q, v: q, lambda t: t, 1, step=-1)


if __name__ == "__main__":
    unittest.main()
