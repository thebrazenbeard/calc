"""Analytic invariants for the geometric-calculus numerical examples."""

import math
import unittest

from calc.geometry import (
    metric_inner, coordinate_jacobian, pullback_metric, christoffel_symbols,
    exterior_derivative_0, exterior_derivative_1,
)


def cartesian_metric(x):
    return ((1.0, 0.0), (0.0, 1.0))


def polar_metric(x):
    r, _ = x
    return ((1.0, 0.0), (0.0, r * r))


class MetricTests(unittest.TestCase):
    def test_metric_inner_polar(self):
        self.assertAlmostEqual(metric_inner(polar_metric, (3, 0), (0, 2), (0, 5)), 90)
        self.assertAlmostEqual(metric_inner(polar_metric, (3, 0), (1, 0), (0, 1)), 0)

    def test_flat_cartesian_christoffel_zero(self):
        gamma = christoffel_symbols(cartesian_metric, (3, 4))
        for plane in gamma:
            for row in plane:
                for value in row:
                    self.assertEqual(value, 0)

    def test_polar_christoffel_known_formula(self):
        # ds² = dr² + r² dθ²; Γ^r_{θθ}=-r, Γ^θ_{rθ}=Γ^θ_{θr}=1/r.
        gamma = christoffel_symbols(polar_metric, (2, 0.3))
        self.assertAlmostEqual(gamma[0][1][1], -2.0, delta=1e-9)
        self.assertAlmostEqual(gamma[1][0][1], 0.5, delta=1e-9)
        self.assertAlmostEqual(gamma[1][1][0], 0.5, delta=1e-9)
        self.assertAlmostEqual(gamma[0][0][0], 0, delta=1e-10)

    def test_polar_coordinate_singularity(self):
        with self.assertRaises(ValueError):
            christoffel_symbols(polar_metric, (0, 0))

    def test_spherical_surface_connection(self):
        # Unit sphere away from coordinate poles: ds²=dθ²+sin²θ dφ².
        sphere = lambda x: ((1.0, 0.0), (0.0, math.sin(x[0])**2))
        theta = math.pi / 3
        gamma = christoffel_symbols(sphere, (theta, 0.4))
        self.assertAlmostEqual(gamma[0][1][1], -math.sin(theta)*math.cos(theta), delta=1e-9)
        self.assertAlmostEqual(gamma[1][0][1], 1/math.tan(theta), delta=1e-9)

    def test_indefinite_metric_not_assumed_positive(self):
        minkowski = lambda x: ((-1.0, 0.0), (0.0, 1.0))
        self.assertEqual(metric_inner(minkowski, (0, 0), (1, 0), (1, 0)), -1)
        self.assertEqual(christoffel_symbols(minkowski, (1, 2))[0][0][0], 0)

    def test_metric_dimensions_symmetry(self):
        with self.assertRaises(ValueError):
            metric_inner(polar_metric, (2, 0), (1,), (1, 0))
        with self.assertRaises(ValueError):
            metric_inner(lambda x: ((1, 3), (2, 1)), (1, 2), (1, 2), (1, 2))
        with self.assertRaises(ValueError):
            christoffel_symbols(lambda x: ((1, 1), (1, 1)), (1, 2))


class CoordinateTests(unittest.TestCase):
    def test_polar_to_cartesian_jacobian(self):
        transform = lambda x: (x[0] * math.cos(x[1]), x[0] * math.sin(x[1]))
        r, angle = 4.0, math.pi / 3
        jac = coordinate_jacobian(transform, (r, angle))
        expected = ((math.cos(angle), -r * math.sin(angle)),
                    (math.sin(angle), r * math.cos(angle)))
        for row, truth in zip(jac, expected):
            for actual, oracle in zip(row, truth):
                self.assertAlmostEqual(actual, oracle, delta=1e-9)

    def test_polar_metric_from_cartesian_pullback(self):
        transform = lambda x: (x[0] * math.cos(x[1]), x[0] * math.sin(x[1]))
        g = pullback_metric(transform, cartesian_metric, (4, math.pi / 3))
        self.assertAlmostEqual(g[0][0], 1.0, delta=1e-9)
        self.assertAlmostEqual(g[1][1], 16.0, delta=1e-9)
        self.assertAlmostEqual(g[0][1], 0.0, delta=1e-9)

    def test_non_square_coordinate_jacobian(self):
        jac = coordinate_jacobian(lambda x: (x[0], x[0] ** 2), (3,))
        self.assertEqual(len(jac), 2)
        self.assertAlmostEqual(jac[1][0], 6, delta=1e-7)


class ExteriorTests(unittest.TestCase):
    def test_exterior_derivative_scalar(self):
        df = exterior_derivative_0(lambda x: x[0] ** 2 * x[1], (3, 2))
        self.assertAlmostEqual(df[0], 12, delta=1e-8)
        self.assertAlmostEqual(df[1], 9, delta=1e-8)

    def test_one_form_curl_component(self):
        # ω = -y dx + x dy => dω = 2 dx∧dy.
        d_omega = exterior_derivative_1(lambda x: (-x[1], x[0]), (3, -1))
        self.assertAlmostEqual(d_omega[(0, 1)], 2, delta=1e-9)

    def test_exact_one_form_has_zero_exterior_derivative(self):
        # ω=d(x²y)=(2xy,x²), analytically dω=0.
        d_omega = exterior_derivative_1(lambda x: (2*x[0]*x[1], x[0]**2), (3, 2))
        self.assertAlmostEqual(d_omega[(0, 1)], 0, delta=1e-9)

    def test_one_dimensional_two_form_empty(self):
        self.assertEqual(exterior_derivative_1(lambda x: (x[0]**2,), (3,)), {})

    def test_invalid_form_and_steps(self):
        with self.assertRaises(ValueError):
            exterior_derivative_1(lambda x: (1,), (1, 2))
        with self.assertRaises(ValueError):
            exterior_derivative_0(lambda x: x[0], (1,), h=0)
        with self.assertRaises(ValueError):
            coordinate_jacobian(lambda x: (x[0],), (1,), h=float('nan'))


if __name__ == "__main__":
    unittest.main()

