"""Stochastic numerics: deterministic witnesses and bounded statistical checks."""

import math
import unittest

from calc.stochastic import (
    brownian_increments,
    brownian_path,
    euler_maruyama,
    geometric_brownian_terminal,
    quadratic_variation,
)


class BrownianTests(unittest.TestCase):
    def test_seeded_increments_are_reproducible(self):
        a = brownian_increments(20, 0.125, seed=37)
        self.assertEqual(a, brownian_increments(20, 0.125, seed=37))
        self.assertNotEqual(a, brownian_increments(20, 0.125, seed=38))

    def test_brownian_scaling_same_standard_normals(self):
        full = brownian_increments(100, 1.0, seed=4)
        half = brownian_increments(100, 0.25, seed=4)
        for a, b in zip(full, half):
            self.assertAlmostEqual(a / 2, b, places=14)

    def test_seeded_moments_near_wiener_expectation(self):
        increments = brownian_increments(20_000, 0.5, seed=429)
        mean = sum(increments) / len(increments)
        variance = sum((x - mean) ** 2 for x in increments) / len(increments)
        self.assertAlmostEqual(mean, 0.0, delta=0.025)
        self.assertAlmostEqual(variance, 0.5, delta=0.02)

    def test_path_includes_initial_position(self):
        increments = brownian_increments(5, 0.2, seed=13)
        path = brownian_path(5, 0.2, x0=3.0, seed=13)
        self.assertEqual(len(path), 6)
        self.assertEqual(path[0], 3.0)
        self.assertAlmostEqual(path[-1], 3.0 + sum(increments))

    def test_known_quadratic_variation(self):
        self.assertEqual(quadratic_variation([0, 1, 0, 2]), 6.0)

    def test_bad_brownian_arguments(self):
        for n in (-1, 0, True, 2.0):
            with self.subTest(n=n), self.assertRaises((TypeError, ValueError)):
                brownian_increments(n, 0.25)
        for dt in (-1, 0, math.inf, math.nan, True):
            with self.subTest(dt=dt), self.assertRaises((TypeError, ValueError)):
                brownian_increments(1, dt)
        with self.assertRaises(ValueError):
            brownian_path(4, 0.2, x0=math.nan)
        with self.assertRaises(ValueError):
            quadratic_variation([0, math.inf])


class ItoSolverTests(unittest.TestCase):
    def test_zero_diffusion_reduces_to_euler(self):
        # y'=y, y(0)=1; 10 steps of explicit Euler have exact value 1.1**10.
        path = euler_maruyama(lambda t, y: y, lambda t, y: 0, 0, 1, 0.1, 10)
        self.assertEqual(len(path), 11)
        self.assertAlmostEqual(path[-1][0], 1)
        self.assertAlmostEqual(path[-1][1], 1.1**10, places=12)

    def test_known_noise_increments(self):
        increments = (0.1, -0.3, 0.2)
        path = euler_maruyama(lambda t, y: 0, lambda t, y: 2, 1, 5, 0.25, 3, increments=increments)
        self.assertEqual([round(p[1], 12) for p in path], [5, 5.2, 4.6, 5.0])
        self.assertEqual([p[0] for p in path], [1.0, 1.25, 1.5, 1.75])

    def test_seeded_solver_reproducible(self):
        run = lambda seed: euler_maruyama(lambda t, y: 0.03 * y, lambda t, y: .2 * y, 0, 2, .02, 50, seed=seed)
        self.assertEqual(run(123), run(123))
        self.assertNotEqual(run(123), run(234))

    def test_exact_gbm_with_ito_correction(self):
        self.assertAlmostEqual(geometric_brownian_terminal(2, 0, 1, 1, 0), 2 * math.exp(-0.5))
        self.assertAlmostEqual(geometric_brownian_terminal(3, .2, 0, 2, 7), 3 * math.exp(.4))
        self.assertEqual(geometric_brownian_terminal(7, .1, .2, 0, 0), 7.0)

    def test_bad_solver_arguments(self):
        for steps in (-1, 0, 1.5, True):
            with self.subTest(steps=steps), self.assertRaises((TypeError, ValueError)):
                euler_maruyama(lambda t, y: 0, lambda t, y: 1, 0, 1, .1, steps)
        with self.assertRaises(ValueError):
            euler_maruyama(lambda t, y: 0, lambda t, y: 1, 0, 1, .1, 2, increments=[0.1])
        with self.assertRaises(ValueError):
            euler_maruyama(lambda t, y: 0, lambda t, y: 1, 0, 1, .1, 2, increments=[math.nan, 0])
        with self.assertRaises(ValueError):
            euler_maruyama(lambda t, y: 0, lambda t, y: 1, 0, 1, .1, 2, increments=[0, 0], seed=4)
        with self.assertRaises(ValueError):
            euler_maruyama(lambda t, y: math.nan, lambda t, y: 1, 0, 1, .1, 2, increments=[0, 0])
        with self.assertRaises(ValueError):
            geometric_brownian_terminal(-1, .2, .3, 1, 0)
        with self.assertRaises(ValueError):
            geometric_brownian_terminal(1, 1000, 0, 2, 0)

    def test_euler_maruyama_approximates_same_noise_gbm(self):
        # Analytic comparison uses the *same Wiener displacement*, not a
        # second independent simulation or a copy of the numerical update.
        increments = brownian_increments(1000, 0.001, seed=289)
        path = euler_maruyama(
            lambda t, s: 0.05 * s, lambda t, s: 0.2 * s,
            0, 2, .001, 1000, increments=increments,
        )
        expected = geometric_brownian_terminal(2, 0.05, 0.2, 1, sum(increments))
        self.assertAlmostEqual(path[-1][1], expected, delta=0.03)


if __name__ == '__main__':
    unittest.main()
