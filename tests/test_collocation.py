"""Independent analytic witnesses for Chebyshev--Lobatto interpolation."""

from __future__ import annotations

import math
import unittest

from calc.collocation import (
    barycentric_interpolate,
    chebyshev_differentiation_matrix,
    chebyshev_lobatto_nodes,
    differentiate_at_nodes,
)


def polynomial(x: float) -> float:
    return x**5 - 2*x*x + 3*x - 1


def polynomial_derivative(x: float) -> float:
    return 5*x**4 - 4*x + 3


class ChebyshevGridTests(unittest.TestCase):
    def test_small_independently_known_grid(self):
        self.assertEqual(chebyshev_lobatto_nodes(1), (1.0, -1.0))
        self.assertEqual(chebyshev_lobatto_nodes(2), (1.0, 0.0, -1.0))
        for actual, expected in zip(chebyshev_lobatto_nodes(3), (1., .5, -.5, -1.)):
            self.assertAlmostEqual(actual, expected, places=14)

    def test_monotone_grid_and_endpoints(self):
        grid = chebyshev_lobatto_nodes(63)
        self.assertEqual((grid[0], grid[-1]), (1., -1.))
        self.assertTrue(all(a > b for a, b in zip(grid, grid[1:])))

    def test_degree_bounds(self):
        for bad in (True, 0, -1, 1.0, 129):
            with self.subTest(bad=bad), self.assertRaises((TypeError, ValueError)):
                chebyshev_lobatto_nodes(bad)


class BarycentricTests(unittest.TestCase):
    def test_interpolation_equals_nodes_by_identity(self):
        nodes = chebyshev_lobatto_nodes(11)
        samples = tuple(math.sin(4*x) for x in nodes)
        for x, y in zip(nodes, samples):
            self.assertEqual(barycentric_interpolate(samples, x), y)

    def test_quintic_reproduction(self):
        nodes = chebyshev_lobatto_nodes(6)
        samples = tuple(polynomial(x) for x in nodes)
        for x in (-.9375, -.5, -.13, 0., .07, .31, .999):
            self.assertAlmostEqual(barycentric_interpolate(samples, x),
                                   polynomial(x), delta=2e-12)

    def test_constant_and_linear(self):
        self.assertAlmostEqual(barycentric_interpolate((4., 4.), .27), 4.)
        self.assertAlmostEqual(barycentric_interpolate((3., -1.), .25), 1.5)

    def test_analytic_exp_converges(self):
        x = .373
        def error(n):
            nodes = chebyshev_lobatto_nodes(n)
            return abs(barycentric_interpolate(tuple(math.exp(t) for t in nodes), x)
                       - math.exp(x))
        coarse, fine = error(4), error(16)
        self.assertGreater(coarse, 1e-6)
        self.assertLess(fine, 1e-12)

    def test_nonfinite_and_bad_data(self):
        for samples, x in (((), 0.), ((1.,), 0.), ((1., 2., float("inf")), 0.),
                           ((1., 2.), float("nan")), ((1., 2.), 1.000001),
                           ((True, 1.), 0.), (("1", 2.), 0.)):
            with self.subTest(samples=samples, x=x), self.assertRaises((TypeError, ValueError)):
                barycentric_interpolate(samples, x)

    def test_too_many_values(self):
        with self.assertRaises(ValueError):
            barycentric_interpolate([1.] * 130, .5)


class DifferentiationTests(unittest.TestCase):
    def test_degree_two_matrix_known_exactly(self):
        matrix = chebyshev_differentiation_matrix(2)
        expected = ((1.5, -2., .5), (.5, 0., -.5), (-.5, 2., -1.5))
        for row, target in zip(matrix, expected):
            for actual, wanted in zip(row, target):
                self.assertAlmostEqual(actual, wanted, delta=1e-14)

    def test_linear_function_degree_one(self):
        self.assertEqual(differentiate_at_nodes((1., -1.)), (1., 1.))

    def test_maximum_supported_degree_linear_witness(self):
        nodes = chebyshev_lobatto_nodes(128)
        result = differentiate_at_nodes(nodes)
        self.assertLess(max(abs(value - 1.0) for value in result), 1e-11)

    def test_quintic_derivative_independent_formula(self):
        nodes = chebyshev_lobatto_nodes(8)
        derivative = differentiate_at_nodes(tuple(polynomial(t) for t in nodes))
        for x, actual in zip(nodes, derivative):
            self.assertAlmostEqual(actual, polynomial_derivative(x), delta=2e-10)

    def test_constant_derivative_zero(self):
        self.assertEqual(differentiate_at_nodes((2.125,) * 32), (0.,) * 32)

    def test_exp_derivative_convergence(self):
        def max_error(n):
            nodes = chebyshev_lobatto_nodes(n)
            actual = differentiate_at_nodes(tuple(math.exp(x) for x in nodes))
            return max(abs(a - math.exp(x)) for a, x in zip(actual, nodes))
        self.assertGreater(max_error(5), 1e-5)
        self.assertLess(max_error(20), 1e-9)

    def test_nonfinite_and_length_boundaries(self):
        for data in ([], [0.], [0., float("nan")], [False, 1.], [1.] * 130):
            with self.subTest(data=data[:5]), self.assertRaises((TypeError, ValueError)):
                differentiate_at_nodes(data)


if __name__ == "__main__":
    unittest.main()
