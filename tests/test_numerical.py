"""Independent analytic checks for numerical calculus routines."""

import math
import unittest

from calc.numerical import adaptive_simpson, euler, rk4, simpson, trapezoid


class QuadratureTests(unittest.TestCase):
    def test_trapezoid_linear_exact(self):
        self.assertAlmostEqual(trapezoid(lambda x: 3*x + 2, -2, 4, 7), 30.0, places=12)

    def test_trapezoid_sine_convergence(self):
        self.assertAlmostEqual(trapezoid(math.sin, 0, math.pi, 1000), 2.0, delta=2e-6)

    def test_simpson_cubic_exact(self):
        # Integral x^3 from -2 to 3 = (3^4 - (-2)^4)/4 = 65/4.
        self.assertAlmostEqual(simpson(lambda x: x**3, -2, 3, 10), 65/4, places=12)

    def test_simpson_reverse_and_zero(self):
        self.assertAlmostEqual(simpson(lambda x: x*x, 1, 0, 10), -1/3, places=12)
        self.assertEqual(simpson(math.exp, 2, 2, 2), 0.0)

    def test_adaptive_simpson_exponential(self):
        self.assertAlmostEqual(adaptive_simpson(math.exp, 0, 1), math.e - 1, delta=1e-10)

    def test_adaptive_simpson_oscillatory(self):
        self.assertAlmostEqual(adaptive_simpson(lambda x: math.cos(8*x), 0, math.pi), 0.0, delta=1e-9)
        self.assertAlmostEqual(adaptive_simpson(lambda x: math.cos(32*x), 0, math.pi), 0.0, delta=1e-9)

    def test_bad_quadrature_arguments(self):
        for n in (0, -4, 1.2, True):
            with self.assertRaises((TypeError, ValueError)):
                trapezoid(math.sin, 0, 1, n)
        for n in (0, -2, 3, 2.5, False):
            with self.assertRaises((TypeError, ValueError)):
                simpson(math.sin, 0, 1, n)
        for tol in (0, -1, float('nan'), float('inf')):
            with self.assertRaises((TypeError, ValueError)):
                adaptive_simpson(math.sin, 0, 1, tol=tol)
        with self.assertRaises(ValueError):
            adaptive_simpson(lambda x: math.exp(x), 0, 1, tol=1e-16, max_depth=0)
        with self.assertRaises(ValueError):
            trapezoid(lambda x: math.inf, 0, 1)
        with self.assertRaises(ValueError):
            simpson(math.sin, 0, float('inf'))


class InitialValueTests(unittest.TestCase):
    def test_euler_exponential(self):
        trajectory = euler(lambda t, y: y, 0, 1, 1, 1000)
        self.assertEqual(len(trajectory), 1001)
        self.assertEqual(trajectory[0], (0.0, 1.0))
        self.assertEqual(trajectory[-1][0], 1.0)
        self.assertAlmostEqual(trajectory[-1][1], math.e, delta=0.002)

    def test_rk4_exponential(self):
        trajectory = rk4(lambda t, y: y, 0, 1, 1, 20)
        self.assertAlmostEqual(trajectory[-1][1], math.e, delta=2e-7)

    def test_vector_harmonic_oscillator(self):
        trajectory = rk4(lambda t, y: (y[1], -y[0]), 0, (0, 1), math.pi/2, 100)
        self.assertIsInstance(trajectory[-1][1], tuple)
        x, v = trajectory[-1][1]
        self.assertAlmostEqual(x, 1.0, delta=1e-8)
        self.assertAlmostEqual(v, 0.0, delta=1e-8)

    def test_reverse_and_zero_intervals(self):
        trajectory = rk4(lambda t, y: 2.0, 2, 5, 1, 4)
        self.assertEqual(trajectory[-1][0], 1.0)
        self.assertAlmostEqual(trajectory[-1][1], 3.0)
        self.assertEqual(euler(lambda t, y: 10, 4, 1, 4, 2),
                         [(4.0, 1.0), (4.0, 1.0), (4.0, 1.0)])

    def test_bad_ode_inputs(self):
        with self.assertRaises(ValueError):
            euler(lambda t, y: y, 0, 1, 1, 0)
        with self.assertRaises(ValueError):
            rk4(lambda t, y: y, 0, (), 1)
        with self.assertRaises(ValueError):
            rk4(lambda t, y: (1, 2), 0, (1,), 1)
        with self.assertRaises(ValueError):
            euler(lambda t, y: float('nan'), 0, 1, 1)
        with self.assertRaises(ValueError):
            rk4(lambda t, y: y, 0, 1, float('nan'))


if __name__ == '__main__':
    unittest.main()
