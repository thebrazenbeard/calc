"""Finite-grid witnesses for scalar Itô and Stratonovich conventions."""

import math
import unittest

from calc.stochastic import brownian_increments, euler_maruyama, geometric_brownian_terminal
from calc.stochastic_conventions import (
    brownian_square_integrals,
    ito_to_stratonovich_drift,
    stratonovich_heun,
    stratonovich_to_ito_drift,
)


class DriftConventionTests(unittest.TestCase):
    def test_multiplicative_noise_drift_conversion(self):
        sigma, mu, y = 0.8, 0.3, 2.5
        a_strat = lambda t, x: mu * x
        b = lambda t, x: sigma * x
        bx = lambda t, x: sigma
        expected = (mu + sigma**2 / 2) * y
        converted = stratonovich_to_ito_drift(a_strat, b, bx, 1.0, y)
        self.assertAlmostEqual(converted, expected)
        self.assertAlmostEqual(
            ito_to_stratonovich_drift(lambda t, x: (mu + sigma**2 / 2) * x, b, bx, 1.0, y),
            mu * y,
        )

    def test_additive_noise_has_no_drift_correction(self):
        self.assertEqual(
            stratonovich_to_ito_drift(lambda t, y: 2, lambda t, y: 3, lambda t, y: 0, 0, 10),
            2.0,
        )

    def test_bad_drift_conversions(self):
        a = lambda t, y: 1.0
        b = lambda t, y: y
        bx = lambda t, y: 1.0
        with self.assertRaises(TypeError):
            stratonovich_to_ito_drift(None, b, bx, 0, 1)
        with self.assertRaises(ValueError):
            stratonovich_to_ito_drift(a, b, bx, math.inf, 1)
        with self.assertRaises(ValueError):
            stratonovich_to_ito_drift(a, b, lambda t, y: math.nan, 0, 1)
        with self.assertRaises(ValueError):
            stratonovich_to_ito_drift(a, lambda t, y: 1e308, lambda t, y: 1e308, 0, 1)


class StochasticHeunTests(unittest.TestCase):
    def test_zero_diffusion_and_constant_drift(self):
        path = stratonovich_heun(lambda t, y: -2, lambda t, y: 0, 0, 5, .25, 4,
                                 increments=[0.3, -0.4, 0.1, 0.2])
        self.assertEqual([t for t, _ in path], [0, .25, .5, .75, 1])
        self.assertAlmostEqual(path[-1][1], 3)

    def test_zero_diffusion_gives_deterministic_heun_method(self):
        path = stratonovich_heun(lambda t, y: y, lambda t, y: 0, 0, 1, .1, 10,
                                 increments=[0] * 10)
        self.assertAlmostEqual(path[-1][1], (1 + .1 + .1**2 / 2)**10, places=12)

    def test_additive_noise_and_constant_drift_is_exact_on_grid(self):
        increments = (.1, -.3, .2)
        path = stratonovich_heun(lambda t, y: 1.5, lambda t, y: -2, 0, 3, .2, 3,
                                 increments=increments)
        self.assertAlmostEqual(path[-1][1], 3 + 1.5 * .6 - 2 * sum(increments))

    def test_heun_and_ito_euler_share_gbm_brownian_path(self):
        # SDE dY=mu*Y dt+sigma*Y circ dW in Stratonovich form. The same
        # process has Itô drift (mu+sigma^2/2)*Y and exact solution
        # y0*exp(mu*T+sigma*W_T), with *the same sampled increments*.
        n, dt, y0, mu, sigma = 2_000, .0005, 2., .1, .3
        increments = brownian_increments(n, dt, seed=289)
        heun = stratonovich_heun(lambda t, y: mu*y, lambda t, y: sigma*y,
                                 0, y0, dt, n, increments=increments)
        ito = euler_maruyama(
            lambda t, y: stratonovich_to_ito_drift(
                lambda tt, yy: mu*yy, lambda tt, yy: sigma*yy,
                lambda tt, yy: sigma, t, y),
            lambda t, y: sigma*y, 0, y0, dt, n, increments=increments,
        )
        exact = y0 * math.exp(mu * (n * dt) + sigma * math.fsum(increments))
        self.assertAlmostEqual(heun[-1][1], exact, delta=.004)
        self.assertAlmostEqual(ito[-1][1], exact, delta=.04)
        self.assertAlmostEqual(
            geometric_brownian_terminal(y0, mu + sigma**2/2, sigma, n*dt, math.fsum(increments)),
            exact,
        )

    def test_invalid_time_noise_and_numerical_states(self):
        f = lambda t, y: 0
        for dt, n in [(0, 1), (-.1, 1), (math.nan, 1), (1, 0), (1, 1.2), (1, True), (1, 1_000_001)]:
            with self.subTest(dt=dt, steps=n), self.assertRaises((TypeError, ValueError)):
                stratonovich_heun(f, f, 0, 0, dt, n, increments=[0])
        for noise in ([0], [0, 0, 0], [math.inf, 0], [math.nan, 0]):
            with self.subTest(noise=noise), self.assertRaises(ValueError):
                stratonovich_heun(f, f, 0, 0, .1, 2, increments=noise)
        with self.assertRaises(ValueError):
            stratonovich_heun(f, f, 1e20, 0, 1e-9, 1, increments=[0])
        with self.assertRaises(ValueError):
            stratonovich_heun(f, f, 1e308, 0, 1e308, 1, increments=[0])
        with self.assertRaises(ValueError):
            stratonovich_heun(lambda t, y: math.nan, f, 0, 0, .1, 1, increments=[0])
        with self.assertRaises(ValueError):
            stratonovich_heun(f, lambda t, y: math.inf, 0, 0, .1, 1, increments=[0])
        with self.assertRaises(ValueError):
            stratonovich_heun(f, lambda t, y: 1e308, 0, 0, .1, 1, increments=[1e308])
        with self.assertRaises(TypeError):
            stratonovich_heun(3, f, 0, 0, 1, 1, increments=[0])


class BrownianSquareWitnessTests(unittest.TestCase):
    def test_finite_grid_ito_correction_identity(self):
        witness = brownian_square_integrals([.1, -.3, .2])
        self.assertAlmostEqual(witness.terminal_displacement, 0)
        self.assertAlmostEqual(witness.quadratic_variation, .14)
        self.assertAlmostEqual(witness.ito_left, -.14)
        self.assertAlmostEqual(witness.stratonovich_midpoint, 0)
        self.assertAlmostEqual(witness.terminal_square,
                               witness.ito_left + witness.quadratic_variation)
        self.assertAlmostEqual(witness.terminal_square, witness.stratonovich_midpoint)

    def test_quadratic_variation_is_only_approximately_elapsed_time(self):
        n, dt = 20_000, 1/20_000
        witness = brownian_square_integrals(brownian_increments(n, dt, seed=429))
        self.assertAlmostEqual(witness.quadratic_variation, n * dt, delta=.05)
        self.assertAlmostEqual(witness.stratonovich_midpoint, witness.terminal_square, delta=1e-10)
        self.assertAlmostEqual(witness.ito_left + witness.quadratic_variation,
                               witness.terminal_square, delta=1e-10)

    def test_invalid_witness_inputs(self):
        with self.assertRaises(ValueError):
            brownian_square_integrals([])
        with self.assertRaises(ValueError):
            brownian_square_integrals([1, math.inf])
        with self.assertRaises(ValueError):
            brownian_square_integrals([1e308, 1e308])
        with self.assertRaises(ValueError):
            brownian_square_integrals([1e154, -1e154])


if __name__ == '__main__':
    unittest.main()
