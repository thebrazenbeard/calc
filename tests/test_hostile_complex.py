"""Independent, adversarial contracts for numerical complex calculus.

The expected values use analytic calculus, parameterization invariants and
representability checks, not the production quadrature implementation.
"""

import cmath
import math
import unittest

from calc.complex_calc import (
    cauchy_riemann_residual,
    complex_derivative,
    contour_integral,
    winding_number,
)


class ComplexHostileContracts(unittest.TestCase):
    def test_holomorphic_direction_independence(self):
        # For an entire f(z)=exp(z), f'(z)=exp(z) in every direction.
        z = 0.4 + 0.7j
        target = cmath.exp(z)
        for direction in (1, 1j, 3 + 4j, -1-1j):
            result = complex_derivative(cmath.exp, z, direction=direction)
            self.assertLess(abs(result - target), 1e-8)

    def test_antiholomorphic_direction_dependence(self):
        f = lambda z: z.conjugate()
        self.assertAlmostEqual(complex_derivative(f, 0.3+0.4j, direction=1).real, 1, delta=1e-10)
        self.assertAlmostEqual(complex_derivative(f, 0.3+0.4j, direction=1j).real, -1, delta=1e-10)
        residual = cauchy_riemann_residual(f, 0.3+0.4j)
        self.assertAlmostEqual(residual[0], 2, delta=1e-10)

    def test_contour_parameterization_for_line_integral(self):
        # ∫(2+z) dz over [0,1] = 5/2, regardless of affine time rescaling.
        first = contour_integral(lambda z: 2+z, lambda t: t, lambda t: 1, 0, 1, panels=8)
        second = contour_integral(lambda z: 2+z, lambda t: 2*t, lambda t: 2, 0, 0.5, panels=8)
        self.assertAlmostEqual(first.real, 2.5, delta=1e-12)
        self.assertAlmostEqual(second.real, 2.5, delta=1e-12)
        self.assertAlmostEqual(first.imag, 0, delta=1e-12)
        self.assertLess(abs(first-second), 1e-12)

    def test_circle_contour_residue_and_orientation(self):
        circle = lambda t: cmath.exp(1j*t)
        tangent = lambda t: 1j*cmath.exp(1j*t)
        integral = contour_integral(lambda z: 1/z, circle, tangent, 0, 2*math.pi, panels=128)
        backwards = contour_integral(lambda z: 1/z, circle, tangent, 2*math.pi, 0, panels=128)
        self.assertLess(abs(integral - 2j*math.pi), 1e-10)
        self.assertLess(abs(integral + backwards), 1e-10)

    def test_winding_of_open_path_cannot_depend_on_translation(self):
        # A genuine winding number is defined only for a closed loop.
        # Moving the same open segment by +10^12 cannot make it closed.
        with self.assertRaises(ValueError):
            winding_number(lambda t: 1000+0.01*t, lambda t: 0.01,
                           center=0, start=0, stop=1)
        with self.assertRaises(ValueError):
            winding_number(lambda t: 1e12+1000+0.01*t, lambda t: 0.01,
                           center=1e12, start=0, stop=1)

    def test_winding_of_closed_circle_cannot_depend_on_translation(self):
        circle = lambda t: cmath.exp(1j*t)
        derivative = lambda t: 1j*cmath.exp(1j*t)
        self.assertEqual(winding_number(circle, derivative, center=0), 1)
        self.assertEqual(winding_number(lambda t: 1e12+circle(t), derivative,
                                        center=1e12), 1)

    def test_unrepresentable_spatial_step_rejected(self):
        # 1 +/- 2**-1074 round to exactly 1.0; a zero numerator is a
        # floating-point artifact, not evidence that derivative(z**2)=0.
        with self.assertRaises(ValueError):
            complex_derivative(lambda z: z*z, 1, step=5e-324)
        with self.assertRaises(ValueError):
            cauchy_riemann_residual(lambda z: z.conjugate(), 1, step=5e-324)

    def test_unrepresentable_parameter_step_rejected(self):
        # start and stop are different finite floats, but (stop-start)/2
        # underflows to zero; quadrature must not silently return zero.
        tiny = 5e-324
        with self.assertRaises(ValueError):
            contour_integral(lambda z: 1, lambda t: complex(t),
                             lambda t: 1, 0, tiny, panels=2)


if __name__ == "__main__":
    unittest.main()

