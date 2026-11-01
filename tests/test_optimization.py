"""Independent mathematical and failure-case witnesses for root/minimization."""

import math
import unittest

from calc.optimization import (
    ConvergenceError, bisection, safeguarded_newton, golden_section_minimize,
)


class RootFindingTests(unittest.TestCase):
    def test_bisection_cubic(self):
        root = bisection(lambda x: x**3 - 2, 1, 2)
        self.assertAlmostEqual(root, 2 ** (1 / 3), delta=1e-11)

    def test_bisection_endpoints_and_reverse(self):
        self.assertEqual(bisection(lambda x: x - 4, 4, 6), 4)
        self.assertEqual(bisection(lambda x: x - 4, 2, 4), 4)
        self.assertAlmostEqual(bisection(lambda x: x - 4, 6, 2), 4)

    def test_no_bracket_rejected(self):
        with self.assertRaises(ValueError):
            bisection(lambda x: x*x + 1, -1, 1)
        with self.assertRaises(ValueError):
            safeguarded_newton(lambda x: x*x + 1, lambda x: 2*x, -1, 1)

    def test_discontinuity_is_not_falsely_root(self):
        jump = lambda x: -1.0 if x < 0.3 else 1.0
        with self.assertRaises(ConvergenceError):
            bisection(jump, 0, 1)
        with self.assertRaises(ConvergenceError):
            safeguarded_newton(jump, lambda _: 1.0, 0, 1)

    def test_tiny_function_scale_does_not_short_circuit_horizontal_accuracy(self):
        f = lambda x: 1e-15 * (x - 123)
        self.assertAlmostEqual(bisection(f, 0, 200), 123, delta=1e-9)
        self.assertAlmostEqual(safeguarded_newton(f, lambda _: 1e-15, 0, 200), 123, delta=1e-9)

    def test_nonfinite_function_and_invalid_parameters(self):
        with self.assertRaises(ValueError):
            bisection(lambda _: math.nan, 0, 1)
        with self.assertRaises(ValueError):
            bisection(lambda x: x - 0.5, 0, 1, abs_tol=0)
        with self.assertRaises(ValueError):
            bisection(lambda x: x - 0.5, 0, 1, max_iter=0)
        with self.assertRaises(ValueError):
            bisection(lambda x: x - 0.5, 0, 1, f_tol=-1)

    def test_max_iteration_failure(self):
        with self.assertRaises(ConvergenceError):
            bisection(lambda x: x - 0.1234567, 0, 1, max_iter=1)

    def test_safeguarded_newton_polynomial(self):
        root = safeguarded_newton(lambda x: x**3 - 2, lambda x: 3*x*x, 0, 2, x0=1)
        self.assertAlmostEqual(root, 2 ** (1 / 3), delta=1e-11)

    def test_safeguarded_newton_zero_derivative_falls_back(self):
        root = safeguarded_newton(lambda x: x**3 - 2, lambda x: 3*x*x, 0, 2, x0=0)
        self.assertAlmostEqual(root, 2 ** (1 / 3), delta=1e-11)

    def test_safeguarded_newton_rejects_initial_outside(self):
        with self.assertRaises(ValueError):
            safeguarded_newton(lambda x: x-1, lambda _: 1, 0, 2, x0=10)
        with self.assertRaises(ValueError):
            safeguarded_newton(lambda x: x, lambda _: 1, 0, 2, x0=10)

    def test_safeguarded_newton_nonfinite_derivative_falls_back(self):
        root = safeguarded_newton(lambda x: x - 0.25, lambda _: math.nan, 0, 1)
        self.assertAlmostEqual(root, 0.25, delta=1e-11)


class MinimizationTests(unittest.TestCase):
    def test_parabola(self):
        x = golden_section_minimize(lambda x: (x - 1.75)**2, -3, 4)
        self.assertAlmostEqual(x, 1.75, delta=1e-8)

    def test_flat_objective_rounding_limitation(self):
        # Adding 2 causes roundoff to erase sub-1e-8 differences in x;
        # minimizer accuracy cannot exceed the objective's resolution.
        x = golden_section_minimize(lambda x: (x - 1.75)**2 + 2, -3, 4)
        self.assertAlmostEqual(x, 1.75, delta=5e-8)

    def test_boundary_minimum(self):
        x = golden_section_minimize(lambda x: x, 1, 3)
        self.assertAlmostEqual(x, 1, delta=1e-8)

    def test_invalid_interval_and_nan(self):
        with self.assertRaises(ValueError):
            golden_section_minimize(lambda x: x*x, 3, 1)
        with self.assertRaises(ValueError):
            golden_section_minimize(lambda _: math.inf, 0, 1)
        with self.assertRaises(ValueError):
            golden_section_minimize(lambda x: x*x, 0, 1, max_iter=0)

    def test_minimization_convergence_limit(self):
        with self.assertRaises(ConvergenceError):
            golden_section_minimize(lambda x: x*x, -2, 1, abs_tol=1e-16, max_iter=1)


if __name__ == "__main__":
    unittest.main()
