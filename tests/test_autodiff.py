"""Analytic derivative oracles for forward-mode automatic differentiation."""

import math
import unittest

from calc.autodiff import (
    Dual,
    cos,
    derivative,
    directional_derivative,
    exp,
    gradient,
    jacobian,
    log,
    sin,
    sqrt,
    value_and_derivative,
)


class DualArithmeticTests(unittest.TestCase):
    def test_product_and_quotient(self):
        x = Dual(3, 1)
        y = (x * x + 2) / (x + 1)
        self.assertAlmostEqual(y.value, 11 / 4)
        self.assertAlmostEqual(y.tangent, 13 / 16)

    def test_polynomial_and_reflected_ops(self):
        self.assertAlmostEqual(derivative(lambda x: 3*x**4 - 7*x + 4, -2), -103)
        self.assertAlmostEqual(derivative(lambda x: 8 / x, 4), -0.5)
        self.assertAlmostEqual(derivative(lambda x: 2**x, 3), 8 * math.log(2))

    def test_variable_exponent(self):
        self.assertAlmostEqual(derivative(lambda x: x**x, 2), 4*(1 + math.log(2)))

    def test_constant_and_zero_tangent(self):
        self.assertEqual(value_and_derivative(lambda x: 5, 23), (5, 0))
        self.assertEqual(derivative(lambda x: x**2, 0), 0)

    def test_trigonometric_exponential_chain(self):
        x = 0.3
        expected = math.exp(x*x)*(math.cos(x)+2*x*math.sin(x))
        self.assertAlmostEqual(derivative(lambda t: exp(t*t)*sin(t), x), expected)
        self.assertAlmostEqual(derivative(lambda t: cos(t), x), -math.sin(x))

    def test_log_and_sqrt(self):
        self.assertAlmostEqual(derivative(lambda x: log(sqrt(x)), 4), 1 / 8)
        self.assertAlmostEqual(derivative(lambda x: sqrt(x), 9), 1 / 6)

    def test_invalid_domains(self):
        with self.assertRaises(ZeroDivisionError):
            derivative(lambda x: 1 / x, 0)
        with self.assertRaises(ValueError):
            derivative(lambda x: log(x), -1)
        with self.assertRaises(ValueError):
            derivative(lambda x: sqrt(x), 0)
        with self.assertRaises(ValueError):
            derivative(lambda x: x**0.5, -2)
        with self.assertRaises(ValueError):
            derivative(lambda x: x**x, -2)
        with self.assertRaises(ValueError):
            derivative(lambda x: x**0, 0)

    def test_invalid_inputs_and_unsupported_math(self):
        for bad in (float("nan"), float("inf"), True, "3"):
            with self.assertRaises((ValueError, TypeError)):
                Dual(bad, 1)
        with self.assertRaises(TypeError):
            derivative(lambda x: math.sin(x), 0.5)


class MultivariableTests(unittest.TestCase):
    def test_gradient_analytic(self):
        g = gradient(lambda v: v[0]**2 + sin(v[1])*v[2], (2, math.pi/2, 3))
        for actual, expected in zip(g, (4, 0, 1)):
            self.assertAlmostEqual(actual, expected)

    def test_jacobian_analytic(self):
        result = jacobian(lambda v: (v[0]*v[1], v[0]**2+v[1], exp(v[1])), (2, 3))
        expected = ((3, 2), (4, 1), (0, math.exp(3)))
        self.assertEqual(len(result), 3)
        for row, target in zip(result, expected):
            for actual, correct in zip(row, target):
                self.assertAlmostEqual(actual, correct)

    def test_directional_one_pass(self):
        answer = directional_derivative(lambda v: v[0]**2 * v[1], (2, 3), (-1, 2))
        self.assertAlmostEqual(answer, -12+8)

    def test_validation(self):
        with self.assertRaises(ValueError):
            gradient(lambda v: 0, ())
        with self.assertRaises(ValueError):
            directional_derivative(lambda v: v[0], (1, 2), (1,))
        with self.assertRaises(TypeError):
            jacobian(lambda v: v[0], (1,))


if __name__ == "__main__":
    unittest.main()
