"""Independent, adversarial mathematical witnesses for Wave 6 algorithms.

This suite deliberately derives expected values without invoking the routines
under review as its own oracle. Fixed seeds make failures reproducible.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import factorial, hypot, isclose
from random import Random
import unittest

from calc.discrete import (
    cumulative_sum, discrete_integral, forward_differences,
    power_sum, summation_by_parts,
)
from calc.forms import Form, basis_form, hodge_star, inner_product, wedge
from calc.nonsmooth import (
    ConvexPiecewiseLinear, absolute_subgradient, hinge_subgradient,
    moreau_envelope_absolute, prox_hinge, soft_threshold,
)
from calc.spectral import (
    eigen_residuals, jacobi_eigen, orthogonality_error, rayleigh_quotient,
)
from calc.stochastic_conventions import (
    brownian_square_integrals, ito_to_stratonovich_drift,
    stratonovich_heun, stratonovich_to_ito_drift,
)
from calc.evidence import get_evidence


class HostileExteriorTests(unittest.TestCase):
    def test_random_rational_hodge_pairing_is_exact(self):
        rng = Random(0xF0A6)
        for n in range(1, 6):
            for degree in range(n + 1):
                basis = tuple(combinations(range(n), degree))
                for _ in range(4):
                    alpha = Form(n, degree, {b: F(rng.randint(-5, 5), 3) for b in basis})
                    beta = Form(n, degree, {b: F(rng.randint(-5, 5), 7) for b in basis})
                    top = tuple(range(n))
                    self.assertEqual(wedge(alpha, hodge_star(beta)).as_dict().get(top, F(0)),
                                     inner_product(alpha, beta))
                    self.assertEqual(wedge(alpha, hodge_star(beta, orientation=-1))
                                     .as_dict().get(top, F(0)), -inner_product(alpha, beta))

    def test_wedge_associativity_and_graded_sign_for_sparse_mixtures(self):
        a = Form(5, 1, {(0,): F(2, 3), (3,): -4})
        b = Form(5, 2, {(1, 2): -F(3, 5), (2, 4): F(7, 3)})
        c = Form(5, 1, {(0,): 2, (4,): -5})
        self.assertEqual(wedge(wedge(a, b), c), wedge(a, wedge(b, c)))
        self.assertEqual(wedge(a, b), wedge(b, a))  # degree 1*2 is even
        self.assertEqual(wedge(a, c), -wedge(c, a))

    def test_hodge_star_involution_and_rational_linearity(self):
        for n in range(1, 6):
            for degree in range(n + 1):
                for basis in combinations(range(n), degree):
                    form = basis_form(n, basis, F(7, 11))
                    self.assertEqual(hodge_star(hodge_star(form)),
                                     (-1 if degree * (n - degree) % 2 else 1) * form)

    def test_form_dimension_cap_prevents_enormous_complement_iteration(self):
        with self.assertRaises(ValueError):
            Form(1_000_000, 0, {(): 1})


class HostileConvexTests(unittest.TestCase):
    def test_seeded_supporting_hyperplane_and_directional_derivative(self):
        rng = Random(0xC0B)
        for _ in range(90):
            affines = tuple((F(rng.randint(-6, 6), 3), F(rng.randint(-8, 8), 4))
                            for _ in range(rng.randint(1, 9)))
            f = ConvexPiecewiseLinear(affines)
            x = F(rng.randint(-8, 8), 2)
            y = F(rng.randint(-8, 8), 2)
            values = [m * x + b for m, b in affines]
            maximum = max(values)
            active = [m for (m, _), value in zip(affines, values) if value == maximum]
            self.assertEqual(f.value(x), maximum)
            self.assertEqual(f.subdifferential(x), (min(active), max(active)))
            for direction in (-3, -1, 0, 1, 4):
                self.assertEqual(f.directional_derivative(x, direction),
                                 max(m * direction for m in active))
            for slope in (min(active), max(active)):
                self.assertGreaterEqual(f.value(y), f.value(x) + slope * (y - x))

    def test_proximal_global_minimality_on_exact_candidate_grid(self):
        candidates = [F(i, 4) for i in range(-24, 25)]
        for weight in (F(1, 4), 1, F(5, 2)):
            for value in (F(-7, 2), -1, 0, F(3, 4), 3):
                for prox, objective, subgrad in (
                    (soft_threshold, lambda z: weight * abs(z) + (z-value)**2 / 2,
                     absolute_subgradient),
                    (prox_hinge, lambda z: weight * max(F(0), z) + (z-value)**2 / 2,
                     hinge_subgradient),
                ):
                    minimizer = prox(value, weight)
                    lo, hi = subgrad(minimizer)
                    self.assertLessEqual(weight * lo, value - minimizer)
                    self.assertGreaterEqual(weight * hi, value - minimizer)
                    self.assertTrue(all(objective(minimizer) <= objective(z)
                                        for z in candidates))

    def test_moreau_envelope_matches_independently_optimized_prox(self):
        for x in (F(-9, 4), -1, 0, F(1, 8), 2):
            for weight in (F(1, 3), 1, F(5, 2)):
                p = soft_threshold(x, weight)
                expected = abs(p) + (p - x)**2 / (2 * weight)
                self.assertEqual(moreau_envelope_absolute(x, weight), expected)


class HostileDiscreteTests(unittest.TestCase):
    def test_newton_finite_difference_of_polynomial_is_factorial(self):
        for degree in range(8):
            samples = [k ** degree for k in range(15)]
            self.assertTrue(all(value == factorial(degree)
                                for value in forward_differences(samples, order=degree)))
            if degree < 7:
                self.assertTrue(all(value == 0 for value in
                                    forward_differences(samples, order=degree + 1)))

    def test_seeded_faulhaber_telescope_including_max_n(self):
        rng = Random(0x5A12)
        for n in [0, 1, 2, 31, 10**12 - 1] + [rng.randrange(0, 250) for _ in range(25)]:
            for power in range(11):
                self.assertEqual(power_sum(n + 1, power) - power_sum(n, power),
                                 n ** power)

    def test_discrete_ftc_and_parts_for_seeded_fraction_sequences(self):
        rng = Random(0xBABA)
        for length in range(1, 10):
            for _ in range(6):
                values = tuple(F(rng.randint(-9, 9), 6) for _ in range(length))
                cumulative = cumulative_sum(values, initial=F(2, 3))
                self.assertEqual(tuple(cumulative[i+1] - cumulative[i]
                                       for i in range(length)), values)
                self.assertEqual(discrete_integral(values, 0, length),
                                 cumulative[-1] - cumulative[0])
                self.assertEqual(discrete_integral(values, length, 0),
                                 cumulative[0] - cumulative[-1])
                u = cumulative
                v = cumulative_sum(tuple(-value for value in values), initial=F(1, 2))
                lhs, boundary, correction = summation_by_parts(u, v)
                independent_lhs = sum((u[k]*(v[k+1]-v[k])
                                       for k in range(length)), F(0))
                independent_correction = sum((v[k+1]*(u[k+1]-u[k])
                                              for k in range(length)), F(0))
                self.assertEqual(lhs, independent_lhs)
                self.assertEqual(correction, independent_correction)
                self.assertEqual(lhs, boundary - correction)


class HostileSpectralTests(unittest.TestCase):
    def test_extreme_diagonal_dynamic_range_preserves_finite_tiny_eigenvalue(self):
        matrix = ((1e-300, 0.), (0., 1e300))
        values, vectors = jacobi_eigen(matrix)
        self.assertEqual(values, (1e-300, 1e300))
        self.assertEqual(rayleigh_quotient(matrix, (1., 0.)), 1e-300)
        self.assertEqual(rayleigh_quotient(matrix, (0., 1.)), 1e300)
        self.assertEqual(orthogonality_error(vectors), 0.)

    def test_seeded_symmetric_two_by_two_matches_closed_form(self):
        rng = Random(0xE19E)
        for _ in range(70):
            a, b, d = (rng.randint(-15, 15) for _ in range(3))
            matrix = ((float(a), float(b)), (float(b), float(d)))
            exact_analytic = sorted(((a + d - hypot(a-d, 2*b)) / 2,
                                     (a + d + hypot(a-d, 2*b)) / 2))
            actual, vectors = jacobi_eigen(matrix)
            self.assertTrue(all(isclose(x, y, rel_tol=1e-11, abs_tol=1e-11)
                                for x, y in zip(actual, exact_analytic)))
            self.assertLess(max(eigen_residuals(matrix, actual, vectors)), 1e-10)
            self.assertLess(orthogonality_error(vectors), 1e-10)

    def test_seeded_three_by_three_invariants(self):
        rng = Random(0xF3A7)
        for _ in range(45):
            a, b, c, d, e, f = [rng.randint(-6, 6) for _ in range(6)]
            matrix = ((a, b, c), (b, d, e), (c, e, f))
            values, vectors = jacobi_eigen(matrix)
            determinant = a*d*f + 2*b*c*e - a*e*e - d*c*c - f*b*b
            self.assertAlmostEqual(sum(values), a+d+f, delta=1e-8)
            self.assertAlmostEqual(values[0]*values[1]*values[2], determinant, delta=1e-7)
            self.assertLess(max(eigen_residuals(matrix, values, vectors)), 1e-9)
            self.assertLess(orthogonality_error(vectors), 1e-9)
            direction = (1.0, -2.0, 3.0)
            rayleigh = rayleigh_quotient(matrix, direction)
            self.assertLessEqual(values[0] - 1e-9, rayleigh)
            self.assertGreaterEqual(values[-1] + 1e-9, rayleigh)


class HostileStochasticTests(unittest.TestCase):
    def test_quadratic_variation_pathwise_identity_against_fraction_oracle(self):
        rng = Random(0x170)
        for _ in range(45):
            increments = tuple(F(rng.randrange(-8, 9), 8)
                               for _ in range(rng.randint(1, 16)))
            path = [F(0)]
            for delta in increments:
                path.append(path[-1] + delta)
            left = sum((2*path[i]*dw for i, dw in enumerate(increments)), F(0))
            midpoint = sum(((path[i]+path[i+1])*dw
                            for i, dw in enumerate(increments)), F(0))
            qv = sum((dw*dw for dw in increments), F(0))
            output = brownian_square_integrals(tuple(float(v) for v in increments))
            self.assertEqual(F(output.terminal_square), path[-1]**2)
            self.assertEqual(F(output.ito_left), left)
            self.assertEqual(F(output.stratonovich_midpoint), midpoint)
            self.assertEqual(F(output.quadratic_variation), qv)
            self.assertEqual(left + qv, midpoint)

    def test_conversion_roundtrip_matches_analytic_state_dependent_noise(self):
        diffusion = lambda t, y: y*y + 2
        derivative = lambda t, y: 2*y
        drift = lambda t, y: 3*y - 1
        for y in (-3., -0.5, 0., 0.75, 2.):
            expected = drift(0., y) + y*(y*y+2)
            converted = stratonovich_to_ito_drift(drift, diffusion, derivative, 0., y)
            self.assertAlmostEqual(converted, expected)
            restored = ito_to_stratonovich_drift(
                lambda t, state: converted, diffusion, derivative, 0., y)
            self.assertAlmostEqual(restored, drift(0., y))

    def test_conversion_recovery_when_intermediate_correction_overflows(self):
        self.assertEqual(
            stratonovich_to_ito_drift(
                lambda t, y: -1e308,
                lambda t, y: 1e308 + 4*y,
                lambda t, y: 4.,
                0., 0.),
            1e308,
        )
        self.assertEqual(
            ito_to_stratonovich_drift(
                lambda t, y: 1e308,
                lambda t, y: 1e308 + 4*y,
                lambda t, y: 4.,
                0., 0.),
            -1e308,
        )

    def test_additive_noise_heun_is_exact_for_constant_coefficients(self):
        path = stratonovich_heun(lambda t, y: 2., lambda t, y: -3.,
                                 0., 1., 0.25, 4,
                                 increments=(0.25, -0.5, 0.125, 0.125))
        self.assertEqual(len(path), 5)
        state = 1.
        for i, increment in enumerate((0.25, -0.5, 0.125, 0.125), start=1):
            state += 2*0.25 - 3*increment
            self.assertAlmostEqual(path[i][1], state)


class HostileEvidenceTests(unittest.TestCase):
    def test_promoted_topic_examples_are_real_modules_and_tests(self):
        examples = {
            "stratonovich": "calc.stochastic_conventions",
            "nonsmooth": "calc.nonsmooth",
            "subdifferential": "calc.nonsmooth",
            "exterior-algebra": "calc.forms",
        }
        for slug, module in examples.items():
            with self.subTest(topic=slug):
                record = get_evidence(slug)
                self.assertEqual(record["status"], "implemented")
                self.assertEqual(record["evidence_level"], "example_tested")
                self.assertFalse(record["formal_proof"])
                self.assertTrue(any(binding["module"] == module
                                    for binding in record["bindings"]))


if __name__ == "__main__":
    unittest.main()
