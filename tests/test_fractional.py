"""Analytic witnesses and negative tests for narrow fractional operators."""

import math
import unittest

from calc.fractional import (
    caputo_power_derivative,
    grunwald_letnikov_derivative,
    riemann_liouville_power_derivative,
)


class PowerFormulaTests(unittest.TestCase):
    def test_caputo_constant_is_zero(self):
        self.assertEqual(caputo_power_derivative(0.5, 0, 2.0), 0.0)
        self.assertEqual(caputo_power_derivative(1, 0, 2.0), 0.0)

    def test_rl_constant_is_not_zero(self):
        self.assertAlmostEqual(riemann_liouville_power_derivative(0.5, 0, 4),
                               1 / (2 * math.sqrt(math.pi)), places=12)

    def test_nonconstant_power_rules_agree(self):
        expected = 4 / math.sqrt(math.pi)  # D^1/2 t at t=4
        self.assertAlmostEqual(caputo_power_derivative(0.5, 1, 4), expected, places=12)
        self.assertAlmostEqual(riemann_liouville_power_derivative(0.5, 1, 4), expected, places=12)

    def test_another_exact_gamma_witness(self):
        # Gamma(3)/Gamma(5/2) = 8/(3 sqrt(pi)) at x=1.
        expected = 8 / (3 * math.sqrt(math.pi))
        self.assertAlmostEqual(caputo_power_derivative(0.5, 2, 1), expected, places=12)

    def test_integer_order_is_ordinary_derivative(self):
        for method in (caputo_power_derivative, riemann_liouville_power_derivative):
            self.assertAlmostEqual(method(1, 3, 2), 12, places=12)
            self.assertEqual(method(1, 0, 0), 0)

    def test_huge_power_gamma_ratio_avoids_cancellation(self):
        for method in (caputo_power_derivative, riemann_liouville_power_derivative):
            self.assertAlmostEqual(method(1, 1e12, 1), 1e12, delta=1e-2)
            # Γ(beta+1)/Γ(beta+1/2) ~ beta**(1/2) as beta→∞.
            value = method(0.5, 1e12, 1)
            self.assertAlmostEqual(value, 1e6, delta=1e-2)

    def test_endpoint_finite_limit(self):
        self.assertEqual(caputo_power_derivative(0.5, 0, 0), 0)
        self.assertEqual(caputo_power_derivative(0.5, 2, 0), 0)
        self.assertAlmostEqual(caputo_power_derivative(0.5, 0.5, 0),
                               math.sqrt(math.pi) / 2, places=12)
        self.assertAlmostEqual(riemann_liouville_power_derivative(0.5, 0.5, 0),
                               math.sqrt(math.pi) / 2, places=12)

    def test_endpoint_singularities_rejected(self):
        with self.assertRaises(ValueError):
            caputo_power_derivative(0.75, 0.5, 0)
        with self.assertRaises(ValueError):
            riemann_liouville_power_derivative(0.5, 0, 0)
        for method in (caputo_power_derivative, riemann_liouville_power_derivative):
            with self.assertRaises(ValueError):
                method(1, 0.5, 0)

    def test_invalid_domain(self):
        for alpha, beta, x in [(0, 1, 1), (1.01, 1, 1), (.5, -1, 1),
                               (.5, 1, -1), (True, 1, 1), (math.inf, 1, 1),
                               (.5, math.nan, 1)]:
            for method in (caputo_power_derivative, riemann_liouville_power_derivative):
                with self.subTest(method=method.__name__, values=(alpha, beta, x)):
                    with self.assertRaises((TypeError, ValueError)):
                        method(alpha, beta, x)


class GrunwaldTests(unittest.TestCase):
    def test_constant_approximates_rl_not_caputo(self):
        value = grunwald_letnikov_derivative(lambda t: 1.0, 1.0, 0.5, steps=4096)
        self.assertAlmostEqual(value, 1 / math.sqrt(math.pi), delta=3e-4)
        self.assertGreater(value, 0.5)

    def test_power_convergence(self):
        exact = 2 / math.gamma(2.5)  # RL half derivative of t^2 at 1
        coarse = grunwald_letnikov_derivative(lambda t: t*t, 1, 0.5, steps=128)
        fine = grunwald_letnikov_derivative(lambda t: t*t, 1, 0.5, steps=2048)
        self.assertLess(abs(fine-exact), abs(coarse-exact))
        self.assertAlmostEqual(fine, exact, delta=1e-3)

    def test_zero_function(self):
        self.assertEqual(grunwald_letnikov_derivative(lambda t: 0, 2, 0.4, steps=16), 0)

    def test_negative_constant_keeps_sign(self):
        value = grunwald_letnikov_derivative(lambda t: -1, 1, .5, steps=2048)
        self.assertAlmostEqual(value, -1 / math.sqrt(math.pi), delta=5e-4)

    def test_input_validation(self):
        invalid = [(1, 0, 10), (0, .5, 20), (1, 1, 20), (1, .5, 0),
                   (1, .5, True), (1, .5, 200001), (math.inf, .5, 20)]
        for x, alpha, n in invalid:
            with self.subTest(x=x, alpha=alpha, n=n):
                with self.assertRaises((TypeError, ValueError)):
                    grunwald_letnikov_derivative(lambda t: t, x, alpha, steps=n)

    def test_nonfinite_function_value(self):
        with self.assertRaises(ValueError):
            grunwald_letnikov_derivative(lambda t: math.inf, 1, .5, steps=10)
        with self.assertRaises(TypeError):
            grunwald_letnikov_derivative(42, 1, .5)

    def test_unrepresentable_step(self):
        with self.assertRaises(ValueError):
            grunwald_letnikov_derivative(lambda t: 1, 5e-324, .5, steps=10)


if __name__ == "__main__":
    unittest.main()
