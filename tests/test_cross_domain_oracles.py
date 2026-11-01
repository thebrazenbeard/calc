"""Independent analytic and exact-rational oracles spanning calculus modules.

These tests deliberately compute expected values from hand-derived formulas,
not from the implementation's own recurrence or private helpers. They are
bounded executable checks, not formal proof for arbitrary inputs.
"""

from fractions import Fraction
import math
import random
import sys
import unittest

from calc.algebra import Polynomial
from calc.collocation import (
    barycentric_interpolate,
    chebyshev_lobatto_nodes,
    differentiate_at_nodes,
)
from calc.fractional import (
    caputo_power_derivative,
    grunwald_letnikov_derivative,
    riemann_liouville_power_derivative,
)
from calc.intervals import Interval
from calc.matrix_calc import frobenius_inner, matrix_gradient, matrix_jvp
from calc.numerical import adaptive_simpson, euler, rk4, simpson, trapezoid
from calc.optimization import bisection, golden_section_minimize, safeguarded_newton
from calc.spectral import eigen_residuals, jacobi_eigen, rayleigh_quotient


class ExactAndQuadratureOracles(unittest.TestCase):
    def test_random_polynomial_antiderivative_exact_endpoint_formula(self):
        rng = random.Random(20261101)
        for _ in range(60):
            coefficients = tuple(Fraction(rng.randint(-7, 7), rng.randint(1, 8))
                                 for _ in range(rng.randint(1, 6)))
            a = Fraction(rng.randint(-4, 0), rng.randint(1, 3))
            b = Fraction(rng.randint(1, 4), rng.randint(1, 3))
            # Integral of Σ c_k*x^k = Σ c_k/(k+1)*(b^(k+1)-a^(k+1)).
            independent = sum((c * (b ** (k + 1) - a ** (k + 1)) / (k + 1)
                               for k, c in enumerate(coefficients)), Fraction(0))
            primitive = Polynomial(coefficients).integrate(Fraction(17, 5))
            self.assertEqual(primitive(b) - primitive(a), independent)
            self.assertEqual(primitive.differentiate(), Polynomial(coefficients))

    def test_simpson_integrates_cubics_and_reverses_signed_bounds(self):
        # The primitive of 3*x^3 - 2*x^2 + 5*x - 7 is known independently.
        def cubic(x):
            return 3*x**3 - 2*x*x + 5*x - 7

        def primitive(x):
            return Fraction(3, 4)*x**4 - Fraction(2, 3)*x**3 + Fraction(5, 2)*x*x - 7*x

        a, b = Fraction(-3, 2), Fraction(9, 4)
        exact = float(primitive(b) - primitive(a))
        for panels in (2, 4, 10, 64):
            with self.subTest(panels=panels):
                self.assertAlmostEqual(simpson(cubic, float(a), float(b), panels), exact, delta=1e-11)
                self.assertAlmostEqual(simpson(cubic, float(b), float(a), panels), -exact, delta=1e-11)

    def test_trapezoid_affine_integral(self):
        a, b = -1.25, 2.75
        # ∫(11x−3) dx = 11*(b²−a²)/2−3*(b−a)
        exact = 11*(b*b - a*a)/2 - 3*(b-a)
        for panels in (1, 2, 7, 128):
            self.assertAlmostEqual(trapezoid(lambda x: 11*x - 3, a, b, panels), exact,
                                   delta=1e-12)

    def test_adaptive_simpson_sine_against_antiderivative(self):
        # Independent integral: [-cos(x)] from 0 to π equals 2.
        self.assertAlmostEqual(adaptive_simpson(math.sin, 0, math.pi, 1e-11, 16), 2,
                               delta=2e-10)


class DifferentialEquationOracles(unittest.TestCase):
    def test_rk4_exponential_global_fourth_order_convergence(self):
        # y'=y, y(0)=1 -> y(t)=exp(t). Error ratio under h -> h/2 ~16.
        errors = [abs(rk4(lambda t, y: y, 0, 1, 1, steps=n)[-1][1] - math.e)
                  for n in (8, 16, 32)]
        self.assertTrue(all(err > 0 for err in errors))
        self.assertGreater(errors[0] / errors[1], 14)
        self.assertGreater(errors[1] / errors[2], 14)
        self.assertLess(errors[-1], 3e-8)

    def test_euler_exponential_global_first_order_convergence(self):
        errors = [abs(euler(lambda t, y: y, 0, 1, 1, steps=n)[-1][1] - math.e)
                  for n in (32, 64, 128)]
        self.assertGreater(errors[0] / errors[1], 1.8)
        self.assertGreater(errors[1] / errors[2], 1.8)

    def test_rk4_rotation_against_trigonometric_solution(self):
        # y1'=-y2, y2'=y1, starting at (1,0) -> (cos t,sin t).
        final = rk4(lambda t, y: (-y[1], y[0]), 0, (1.0, 0.0), math.pi/2, 64)[-1][1]
        self.assertAlmostEqual(final[0], 0.0, delta=1e-7)
        self.assertAlmostEqual(final[1], 1.0, delta=1e-7)


class SpectralAndCollocationOracles(unittest.TestCase):
    def test_chebyshev_nodes_polynomial_derivative_with_closed_form(self):
        for degree in (5, 8, 16):
            nodes = chebyshev_lobatto_nodes(degree)
            polynomial = lambda x: x**5 - 3*x*x + 2*x - 7
            derivative = lambda x: 5*x**4 - 6*x + 2
            samples = tuple(polynomial(x) for x in nodes)
            got = differentiate_at_nodes(samples)
            for x, estimate in zip(nodes, got):
                self.assertAlmostEqual(estimate, derivative(x), delta=2e-9)
            for x in (-0.91, -0.13, 0.37, 0.83):
                self.assertAlmostEqual(barycentric_interpolate(samples, x), polynomial(x),
                                       delta=5e-12)

    def test_chebyshev_interpolation_of_exponential_improves_with_degree(self):
        errors = []
        for degree in (4, 8, 16):
            nodes = chebyshev_lobatto_nodes(degree)
            samples = tuple(math.exp(x) for x in nodes)
            errors.append(max(abs(barycentric_interpolate(samples, x)-math.exp(x))
                              for x in (-0.94, -0.77, -0.31, 0.21, 0.59, 0.92)))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[-1], 1e-12)

    def test_two_by_two_symmetric_eigenvalues_from_trace_and_determinant(self):
        # Matrix [[2,1],[1,2]] has roots of λ²−4λ+3: eigenvalues 1 and 3.
        matrix = ((2.0, 1.0), (1.0, 2.0))
        eigenvalues, vectors = jacobi_eigen(matrix)
        self.assertAlmostEqual(eigenvalues[0], 1.0, delta=1e-12)
        self.assertAlmostEqual(eigenvalues[1], 3.0, delta=1e-12)
        self.assertLess(max(eigen_residuals(matrix, eigenvalues, vectors)), 1e-12)
        self.assertAlmostEqual(rayleigh_quotient(matrix, (1, -1)), 1, delta=1e-12)
        self.assertAlmostEqual(rayleigh_quotient(matrix, (1, 1)), 3, delta=1e-12)

    def test_spectral_dynamic_range_preserves_representable_diagonal_eigenvalue(self):
        matrix = ((1e-300, 0.0), (0.0, 1e300))
        values, _ = jacobi_eigen(matrix)
        self.assertEqual(values, (1e-300, 1e300))
        self.assertEqual(rayleigh_quotient(matrix, (1, 0)), 1e-300)


class FractionalAndIntervalOracles(unittest.TestCase):
    def test_half_derivatives_from_gamma_closed_forms(self):
        # Γ(3/2)=√π/2, Γ(5/2)=3√π/4.
        self.assertAlmostEqual(caputo_power_derivative(.5, .5, 4), math.sqrt(math.pi)/2,
                               delta=1e-13)
        self.assertAlmostEqual(caputo_power_derivative(.5, 1.5, 4),
                               3*math.sqrt(math.pi), delta=1e-12)
        self.assertEqual(caputo_power_derivative(.5, 0, 4), 0.0)
        self.assertAlmostEqual(riemann_liouville_power_derivative(.5, 0, 4),
                               1/(2*math.sqrt(math.pi)), delta=1e-13)

    def test_grunwald_constant_approaches_rl_not_caputo(self):
        exact = 1/math.sqrt(math.pi)
        errors = [abs(grunwald_letnikov_derivative(lambda x: 1., 1., .5, steps=n)-exact)
                  for n in (256, 1024, 4096)]
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[-1], 8e-5)

    def test_seeded_binary64_interval_products_contain_exact_rational_extrema(self):
        rng = random.Random(960721)
        for _ in range(70):
            left = sorted((rng.randint(-200, 200)/7, rng.randint(-200, 200)/11))
            right = sorted((rng.randint(-200, 200)/13, rng.randint(-200, 200)/17))
            enclosure = Interval(*left)*Interval(*right)
            exact_values = tuple(Fraction.from_float(x)*Fraction.from_float(y)
                                 for x in left for y in right)
            self.assertLessEqual(Fraction.from_float(enclosure.lower), min(exact_values))
            self.assertGreaterEqual(Fraction.from_float(enclosure.upper), max(exact_values))

    def test_interval_large_cancellation_does_not_falsely_narrow(self):
        width = sys.float_info.max
        enclosure = Interval.point(width) * Interval.point(2.0)
        self.assertEqual(enclosure.lower, width)
        self.assertEqual(enclosure.upper, math.inf)
        with self.assertRaises(ValueError):
            _ = enclosure + Interval.point(1.0)
        with self.assertRaises(ZeroDivisionError):
            _ = Interval(-1.0, 2.0) / Interval(-1.0, 1.0)


class MatrixAndOptimizationOracles(unittest.TestCase):
    def test_frobenius_inner_exact_rational_dot(self):
        a = ((Fraction(1, 3), Fraction(-7, 5)), (Fraction(4), Fraction(2, 11)))
        b = ((Fraction(9, 2), Fraction(3)), (Fraction(-2, 7), Fraction(11)))
        # Hand expansion uses four scalar products and no matrix implementation.
        exact = Fraction(3, 2) - Fraction(21, 5) - Fraction(8, 7) + 2
        self.assertEqual(frobenius_inner(a, b), exact)

    def test_matrix_gradient_and_jvp_of_affine_maps(self):
        # f(X)=3 X00 -7 X01 +2 X10 +9 X11; Jacobian is constant.
        coeff = ((3.0, -7.0), (2.0, 9.0))
        x = ((.75, -.5), (1.5, 4.0))
        objective = lambda a: sum(coeff[i][j]*a[i][j] for i in range(2) for j in range(2))
        estimate = matrix_gradient(objective, x, step=1e-4)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(estimate[i][j], coeff[i][j], delta=3e-10)
        tangent = ((1.0, 2.0), (-3.0, 4.0))
        transform = lambda a: ((2*a[0][0]+3*a[1][1], a[0][1]-a[1][0]),)
        self.assertEqual(tuple(len(row) for row in matrix_jvp(transform, x, tangent)), (2,))
        jvp = matrix_jvp(transform, x, tangent, step=1e-4)
        self.assertAlmostEqual(jvp[0][0], 14.0, delta=2e-10)
        self.assertAlmostEqual(jvp[0][1], 5.0, delta=2e-10)

    def test_bracketed_roots_and_unimodal_minimum_against_closed_forms(self):
        root = 2**(1/3)
        self.assertAlmostEqual(bisection(lambda x: x*x*x-2, 0, 2), root, delta=2e-11)
        self.assertAlmostEqual(safeguarded_newton(lambda x: x*x-2,
                                                   lambda x: 2*x, 0, 2),
                               math.sqrt(2), delta=2e-11)
        argmin = golden_section_minimize(lambda x: (x-Fraction(2, 7))**2+3,
                                         -2.0, 2.0)
        self.assertAlmostEqual(argmin, 2/7, delta=2e-8)

    def test_smooth_sqrt_two_bracket_must_not_stop_before_small_residual(self):
        # A small x-bracket is not an exhausted iteration budget. With
        # f(x)=x²−2, the requested default residual is attainable after a
        # few more binary64 steps past a 1e-12-wide bracket.
        expected = math.sqrt(2)
        f = lambda x: x*x - 2
        slope = lambda x: 2*x
        self.assertAlmostEqual(bisection(f, 0, 2), expected, delta=2e-11)
        self.assertAlmostEqual(safeguarded_newton(f, slope, 0, 2), expected,
                               delta=2e-11)

    def test_invalid_domain_contracts(self):
        with self.assertRaises(ValueError):
            simpson(lambda x: x, 0, 1, 3)
        with self.assertRaises(ValueError):
            jacobi_eigen(((1.0, 2.0), (3.0, 4.0)))
        with self.assertRaises(ValueError):
            grunwald_letnikov_derivative(lambda x: x, 1, 1.0)
        with self.assertRaises(ValueError):
            barycentric_interpolate((0.0, 1.0), 2.0)


if __name__ == "__main__":
    unittest.main()
