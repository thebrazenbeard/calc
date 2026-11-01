"""Independent mathematical witnesses for numerical complex calculus."""

import cmath
import math
import unittest

from calc.complex_calc import (
    cauchy_riemann_residual,
    complex_derivative,
    contour_integral,
    winding_number,
)


def circle(radius=1.0, center=0j, orientation=1):
    def path(t):
        return center + radius * cmath.exp(1j * orientation * t)

    def tangent(t):
        return 1j * orientation * radius * cmath.exp(1j * orientation * t)

    return path, tangent


class ComplexDerivativeTests(unittest.TestCase):
    def test_polynomial_derivative_independent_of_direction(self):
        point = 1.25 - 2.0j
        expected = 2 * point
        for direction in (1 + 0j, 0 + 1j, 1 + 1j):
            with self.subTest(direction=direction):
                result = complex_derivative(lambda z: z*z, point, direction=direction)
                self.assertAlmostEqual(result.real, expected.real, delta=1e-8)
                self.assertAlmostEqual(result.imag, expected.imag, delta=1e-8)

    def test_exponential(self):
        point = 0.7 + 0.3j
        expected = cmath.exp(point)
        result = complex_derivative(cmath.exp, point)
        self.assertAlmostEqual(result.real, expected.real, delta=1e-8)
        self.assertAlmostEqual(result.imag, expected.imag, delta=1e-8)

    def test_conjugation_has_direction_dependent_derivative(self):
        f = lambda z: z.conjugate()
        self.assertAlmostEqual(complex_derivative(f, 1 + 2j, direction=1).real, 1, places=8)
        self.assertAlmostEqual(complex_derivative(f, 1 + 2j, direction=1j).real, -1, places=8)

    def test_cauchy_riemann_residuals(self):
        residual = cauchy_riemann_residual(lambda z: z**3, 0.3 + 0.8j)
        self.assertLess(max(abs(component) for component in residual), 1e-7)
        nonanalytic = cauchy_riemann_residual(lambda z: z.conjugate(), 1 + 2j)
        self.assertAlmostEqual(nonanalytic[0], 2, places=7)
        self.assertAlmostEqual(nonanalytic[1], 0, places=7)

    def test_derivative_errors(self):
        for kwargs in ({"step": 0}, {"step": -1}, {"step": 10**1000},
                       {"direction": 0}, {"direction": complex(float("nan"), 1)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                complex_derivative(lambda z: z, 1j, **kwargs)
        with self.assertRaises(ValueError):
            complex_derivative(lambda z: complex(float("inf"), 0), 0)


class ContourTests(unittest.TestCase):
    def test_residue_unit_circle_positive_and_negative(self):
        for orientation in (1, -1):
            with self.subTest(orientation=orientation):
                path, tangent = circle(orientation=orientation)
                integral = contour_integral(lambda z: 1/z, path, tangent, 0, 2*math.pi, panels=256)
                self.assertAlmostEqual(integral.real, 0, delta=1e-9)
                self.assertAlmostEqual(integral.imag, orientation * 2*math.pi, delta=1e-8)

    def test_holomorphic_polynomial_has_zero_integral(self):
        path, tangent = circle(radius=2, center=1 + 2j)
        integral = contour_integral(lambda z: z*z, path, tangent, 0, 2*math.pi)
        self.assertLess(abs(integral), 1e-8)

    def test_winding(self):
        for orientation in (1, -1):
            path, tangent = circle(radius=2, center=1 + 2j, orientation=orientation)
            self.assertEqual(winding_number(path, tangent, 1 + 2j), orientation)
            self.assertEqual(winding_number(path, tangent, -10j), 0)
        path = lambda t: cmath.exp(2j*t)
        tangent = lambda t: 2j*cmath.exp(2j*t)
        self.assertEqual(winding_number(path, tangent, 0j), 2)

    def test_open_path_and_pole(self):
        with self.assertRaises(ValueError):
            winding_number(lambda t: 1 + t, lambda _: 1, 0j)
        path, tangent = circle(center=0j)
        with self.assertRaises(ValueError):
            winding_number(path, tangent, 1 + 0j)

    def test_invalid_contour_inputs(self):
        path, tangent = circle()
        for kwargs in ({"panels": 3}, {"panels": 0}, {"panels": 2.5}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                contour_integral(lambda z: z, path, tangent, 0, 2*math.pi, **kwargs)
        with self.assertRaises(ValueError):
            contour_integral(lambda z: complex(float("nan"), 0), path, tangent, 0, 2*math.pi)
        with self.assertRaises(ValueError):
            contour_integral(lambda z: z, path, tangent, math.inf, 2*math.pi)


if __name__ == "__main__":
    unittest.main()
