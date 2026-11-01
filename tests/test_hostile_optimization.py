"""Independent adversarial contracts for scalar root/minimization methods.

The expected values come from algebra and monotonicity, not the solver code.
"""

import math
import unittest

from calc.optimization import (
    ConvergenceError,
    bisection,
    golden_section_minimize,
    safeguarded_newton,
)


class HostileRootContracts(unittest.TestCase):
    def test_bisection_extreme_symmetric_bracket(self):
        # Opposite-sign finite endpoints whose naive width overflows float.
        root = bisection(lambda x: x / 1e308, -1e308, 1e308)
        self.assertEqual(root, 0.0)

    def test_bisection_nonzero_discontinuity_not_root(self):
        def sign_jump(x):
            return -0.25 if x < math.pi / 10 else 0.75
        with self.assertRaises(ConvergenceError):
            bisection(sign_jump, 0, 1, f_tol=1e-8)

    def test_newton_exceptional_derivative_falls_back(self):
        def unusable_derivative(_):
            raise ZeroDivisionError("synthetic derivative failure")
        root = safeguarded_newton(
            lambda x: x - 0.125, unusable_derivative, 0, 1
        )
        self.assertAlmostEqual(root, 0.125, delta=1e-11)

    def test_quadratic_double_root_requires_sign_bracket(self):
        # An even-multiplicity interior zero does not supply a sign bracket.
        with self.assertRaises(ValueError):
            bisection(lambda x: (x - 0.5) ** 2, 0, 1)

    def test_nonfinite_function_after_initial_valid_bracket_rejected(self):
        def f(x):
            if 0 < x < 1:
                return math.nan
            return x - 0.5
        with self.assertRaises(ValueError):
            bisection(f, 0, 1)

    def test_newton_zero_derivative_at_start_preserves_bracket(self):
        # A zero slope at x0 is not a root; bisection must take over.
        root = safeguarded_newton(
            lambda x: x**3 - 0.25, lambda x: 3*x*x,
            0, 1, x0=0,
        )
        self.assertAlmostEqual(root, 0.25 ** (1 / 3), delta=1e-10)


class HostileMinimizationContracts(unittest.TestCase):
    def test_shifted_parabola_and_boundary(self):
        x = golden_section_minimize(lambda t: 7*(t + 1.25)**2 - 2, -3, 10)
        self.assertAlmostEqual(x, -1.25, delta=1e-7)
        x_boundary = golden_section_minimize(lambda t: -3*t + 4, -1, 8)
        self.assertAlmostEqual(x_boundary, 8, delta=1e-7)

    def test_finite_endpoints_never_expose_nonfinite_query(self):
        # Previous direct reproduction passed +/-inf and NaN to the user's
        # function and nevertheless returned a nominal argmin.
        seen = []

        def constant_objective(x):
            seen.append(x)
            self.assertTrue(math.isfinite(x), "objective evaluated at nonfinite coordinate")
            return 0.0

        try:
            argmin = golden_section_minimize(
                constant_objective, -1e308, 1e308, max_iter=30
            )
        except (ValueError, ConvergenceError):
            # Rejection of an unrepresentably wide bracket is also valid.
            pass
        else:
            self.assertTrue(math.isfinite(argmin))
        self.assertTrue(all(math.isfinite(x) for x in seen))

    def test_invalid_boolean_parameters(self):
        with self.assertRaises(TypeError):
            golden_section_minimize(lambda x: x*x, 0, 1, abs_tol=True)
        with self.assertRaises(ValueError):
            bisection(lambda x: x - 0.5, 0, 1, max_iter=True)


if __name__ == "__main__":
    unittest.main()
