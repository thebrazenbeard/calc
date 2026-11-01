"""Independent analytic checks for numerical multivariable calculus."""

import math
import unittest

from calc.differential import (
    curl,
    derivative,
    directional_derivative,
    divergence,
    gradient,
    hessian,
    jacobian,
    laplacian,
    partial_derivative,
)


class DifferentialCalculusTests(unittest.TestCase):
    def test_scalar_derivative(self):
        self.assertAlmostEqual(derivative(lambda x: x**3, 2.0), 12.0, places=5)
        self.assertAlmostEqual(derivative(math.sin, math.pi / 3), 0.5, places=6)
        self.assertAlmostEqual(derivative(math.exp, 0.0), 1.0, places=6)

    def test_forward_and_backward(self):
        for method in ("forward", "backward"):
            with self.subTest(method=method):
                self.assertAlmostEqual(
                    derivative(lambda x: x * x, 3.0, method=method, step=1e-6),
                    6.0,
                    places=4,
                )

    def test_scalar_rejects_invalid_inputs(self):
        for invalid in (0.0, -0.1, math.inf, math.nan):
            with self.subTest(step=invalid), self.assertRaises(ValueError):
                derivative(math.sin, 1.0, step=invalid)
        with self.assertRaises(ValueError):
            derivative(math.sin, 1.0, method="imaginary")
        with self.assertRaises(ValueError):
            derivative(math.sin, math.inf)
        with self.assertRaises(ValueError):
            derivative(lambda x: math.nan, 1.0)

    def test_partial_gradient_do_not_mutate_point(self):
        point = [2.0, 1.0]
        func = lambda p: p[0] ** 2 + 3 * p[0] * p[1] + p[1] ** 3
        self.assertAlmostEqual(partial_derivative(func, point, 0), 7.0, places=5)
        self.assertAlmostEqual(partial_derivative(func, point, 1), 9.0, places=5)
        grad = gradient(func, point)
        self.assertEqual(len(grad), 2)
        self.assertAlmostEqual(grad[0], 7.0, places=5)
        self.assertAlmostEqual(grad[1], 9.0, places=5)
        self.assertEqual(point, [2.0, 1.0])

    def test_jacobian(self):
        funcs = (lambda p: p[0] ** 2 + p[1], lambda p: math.exp(p[0] * p[1]))
        result = jacobian(funcs, (1.0, 2.0))
        expected = ((2.0, 1.0), (2.0 * math.exp(2), math.exp(2)))
        for row, expected_row in zip(result, expected):
            for value, answer in zip(row, expected_row):
                self.assertAlmostEqual(value, answer, places=4)

    def test_hessian(self):
        func = lambda p: 3 * p[0] ** 2 + 2 * p[0] * p[1] + 4 * p[1] ** 2
        result = hessian(func, (1.25, -2.5))
        for row, expected in zip(result, ((6.0, 2.0), (2.0, 8.0))):
            for value, answer in zip(row, expected):
                self.assertAlmostEqual(value, answer, places=3)
        self.assertEqual(result[0][1], result[1][0])

    def test_directional_derivative(self):
        f = lambda p: 3 * p[0] + 4 * p[1]
        self.assertAlmostEqual(directional_derivative(f, (2, 1), (3, 4)), 5.0, places=5)
        self.assertAlmostEqual(
            directional_derivative(f, (2, 1), (3, 4), normalize=False),
            25.0,
            places=5,
        )

    def test_vector_calculus(self):
        field = (lambda p: p[0] ** 2, lambda p: 2 * p[1], lambda p: p[2] ** 2)
        self.assertAlmostEqual(divergence(field, (2, 3, 5)), 16.0, places=5)
        rotation = (lambda p: -p[1] / 2, lambda p: p[0] / 2, lambda p: p[2])
        actual = curl(rotation, (2, 3, 5))
        for value, answer in zip(actual, (0.0, 0.0, 1.0)):
            self.assertAlmostEqual(value, answer, places=5)
        self.assertAlmostEqual(
            laplacian(lambda p: sum(x**2 for x in p), (2, 3, 5)), 6.0, places=3
        )

    def test_dimension_errors(self):
        with self.assertRaises(ValueError):
            gradient(lambda p: 1.0, [])
        with self.assertRaises(IndexError):
            partial_derivative(lambda p: 1.0, (1, 2), 2)
        with self.assertRaises(ValueError):
            directional_derivative(lambda p: 1.0, (1, 2), (1,))
        with self.assertRaises(ValueError):
            directional_derivative(lambda p: 1.0, (1, 2), (0, 0))
        with self.assertRaises(ValueError):
            divergence((lambda p: p[0],), (1, 2))
        with self.assertRaises(ValueError):
            curl((lambda p: p[0],), (1, 2, 3))


if __name__ == "__main__":
    unittest.main()
