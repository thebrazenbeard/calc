"""Independent finite time-scale identities and negative cases."""

import math
import unittest

from calc.time_scales import FiniteTimeScale


class FiniteTimeScaleTests(unittest.TestCase):
    def setUp(self):
        self.uniform = FiniteTimeScale((0, 1, 2, 3))

    def test_forward_and_backward_jumps(self):
        self.assertEqual(self.uniform.sigma(1), 2)
        self.assertEqual(self.uniform.rho(1), 0)
        self.assertEqual(self.uniform.sigma(3), 3)
        self.assertEqual(self.uniform.rho(0), 0)

    def test_graininess_and_scattered_predicates(self):
        self.assertEqual(self.uniform.mu(1), 1)
        self.assertEqual(self.uniform.nu(1), 1)
        self.assertEqual(self.uniform.mu(3), 0)
        self.assertEqual(self.uniform.nu(0), 0)
        self.assertTrue(self.uniform.is_right_scattered(1))
        self.assertFalse(self.uniform.is_right_scattered(3))
        self.assertTrue(self.uniform.is_left_scattered(3))
        self.assertFalse(self.uniform.is_left_scattered(0))

    def test_delta_nabla_of_square_on_uniform_grid(self):
        # (t+1)^2-t^2=2t+1, t^2-(t-1)^2=2t-1.
        f = lambda t: t * t
        self.assertEqual(self.uniform.delta_derivative(f, 2), 5)
        self.assertEqual(self.uniform.nabla_derivative(f, 2), 3)

    def test_irregular_grid_changes_denominator(self):
        scale = FiniteTimeScale((0, .5, 2))
        f = lambda t: t*t
        self.assertEqual(scale.delta_derivative(f, 0), .5)
        self.assertEqual(scale.delta_derivative(f, .5), 2.5)
        self.assertEqual(scale.nabla_derivative(f, 2), 2.5)

    def test_q_grid_independent_jackson_witness(self):
        # On positive q-grid t, sigma(t)=q*t, D_delta t^2=(q+1)t.
        scale = FiniteTimeScale((1, 2, 4, 8))
        for t in (1, 2, 4):
            self.assertEqual(scale.delta_derivative(lambda x: x*x, t), 3*t)

    def test_delta_integral_left_samples(self):
        # Integral_0^3 t Delta t = 0 + 1 + 2 = 3 (not 9/2).
        self.assertEqual(self.uniform.delta_integral(lambda t: t, 0, 3), 3)
        self.assertEqual(self.uniform.delta_integral(lambda t: t, 3, 0), -3)

    def test_nabla_integral_right_samples(self):
        # Integral_0^3 t Nabla t = 1 + 2 + 3 = 6.
        self.assertEqual(self.uniform.nabla_integral(lambda t: t, 0, 3), 6)
        self.assertEqual(self.uniform.nabla_integral(lambda t: t, 3, 0), -6)

    def test_discrete_fundamental_theorem(self):
        scale = FiniteTimeScale((1, 2, 4, 8))
        f = lambda t: t*t
        self.assertEqual(scale.delta_integral(lambda t: scale.delta_derivative(f, t), 1, 8), 63)
        self.assertEqual(scale.nabla_integral(lambda t: scale.nabla_derivative(f, t), 1, 8), 63)

    def test_shifted_product_rules(self):
        # Time-scale product rules include a shifted factor, unlike the
        # familiar unshifted continuous Leibniz identity.
        grid = FiniteTimeScale((1, 2, 4))
        f, g = lambda t: t, lambda t: t*t
        for t in (1, 2):
            delta_product = grid.delta_derivative(lambda x: f(x)*g(x), t)
            shifted_rhs = grid.delta_derivative(f, t)*g(t) + f(grid.sigma(t))*grid.delta_derivative(g, t)
            self.assertEqual(delta_product, shifted_rhs)
        for t in (2, 4):
            nabla_product = grid.nabla_derivative(lambda x: f(x)*g(x), t)
            shifted_rhs = grid.nabla_derivative(f, t)*g(t) + f(grid.rho(t))*grid.nabla_derivative(g, t)
            self.assertEqual(nabla_product, shifted_rhs)

    def test_integrals_on_equal_endpoints(self):
        def should_not_be_called(_):
            self.fail("zero interval should not evaluate the integrand")
        self.assertEqual(self.uniform.delta_integral(should_not_be_called, 1, 1), 0)
        self.assertEqual(self.uniform.nabla_integral(should_not_be_called, 1, 1), 0)

    def test_derivative_excludes_terminal_point(self):
        with self.assertRaises(ValueError):
            self.uniform.delta_derivative(lambda t: t, 3)
        with self.assertRaises(ValueError):
            self.uniform.nabla_derivative(lambda t: t, 0)

    def test_singleton_has_jumps_but_no_derivatives(self):
        scale = FiniteTimeScale((4,))
        self.assertEqual(scale.sigma(4), 4)
        self.assertEqual(scale.rho(4), 4)
        with self.assertRaises(ValueError):
            scale.delta_derivative(lambda t: t, 4)
        with self.assertRaises(ValueError):
            scale.nabla_derivative(lambda t: t, 4)

    def test_rejects_unsorted_repeated_nonfinite_and_boolean_grid(self):
        for points in ((), (1, 1), (2, 1), (1, float('nan')),
                       (1, float('inf')), (False, 1), (10**25, 10**25 + 1)):
            with self.subTest(points=points), self.assertRaises((ValueError, TypeError)):
                FiniteTimeScale(points)

    def test_requires_exact_grid_membership(self):
        with self.assertRaises(ValueError):
            self.uniform.sigma(1.01)
        with self.assertRaises((TypeError, ValueError)):
            self.uniform.mu('1')

    def test_invalid_function_values_and_callback(self):
        with self.assertRaises((ValueError, TypeError)):
            self.uniform.delta_derivative(lambda t: math.inf, 1)
        with self.assertRaises((ValueError, TypeError)):
            self.uniform.nabla_integral(lambda t: math.nan, 0, 3)
        with self.assertRaises(TypeError):
            self.uniform.delta_integral(42, 0, 3)

    def test_bounded_output_uses_finite_result(self):
        with self.assertRaises(ValueError):
            self.uniform.delta_derivative(lambda t: 1e308 if t else -1e308, 0)


if __name__ == '__main__':
    unittest.main()
