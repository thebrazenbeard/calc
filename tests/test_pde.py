"""Independent mathematical witnesses for one-dimensional PDE solvers."""

import math
import unittest

from calc.pde import heat_1d_dirichlet, wave_1d_dirichlet


class HeatEquationTests(unittest.TestCase):
    def test_zero_is_an_exact_solution(self):
        states = heat_1d_dirichlet([0] * 11, diffusivity=2, dx=.1, dt=.002, steps=20)
        self.assertEqual(len(states), 21)
        self.assertTrue(all(u == (0.0,) * 11 for u in states))

    def test_linear_steady_state_unequal_dirichlet_boundaries(self):
        initial = [1 + j / 10 for j in range(11)]
        states = heat_1d_dirichlet(initial, diffusivity=1, dx=.1, dt=.004, steps=200, left=1, right=2)
        for state in states:
            for actual, expected in zip(state, initial):
                self.assertAlmostEqual(actual, expected, delta=1e-12)

    def test_sine_mode_against_continuum_solution(self):
        # u(x,t) = exp(-pi^2*t) sin(pi*x) on [0,1] with zero boundaries.
        n, dx, dt, steps = 101, .01, .00002, 1000
        initial = [math.sin(math.pi * j / (n - 1)) for j in range(n)]
        states = heat_1d_dirichlet(initial, diffusivity=1, dx=dx, dt=dt, steps=steps)
        t = dt * steps
        self.assertAlmostEqual(states[-1][50], math.exp(-math.pi**2 * t), delta=2e-5)
        self.assertEqual(states[-1][0], 0)
        self.assertEqual(states[-1][-1], 0)

    def test_single_sine_step_matches_discrete_eigenvalue(self):
        n, dx, r = 11, .1, .2
        dt = r * dx**2
        initial = [math.sin(math.pi * j / (n-1)) for j in range(n)]
        actual = heat_1d_dirichlet(initial, diffusivity=1, dx=dx, dt=dt, steps=1)[1]
        amplification = 1 - 4 * r * math.sin(math.pi / (2 * (n - 1))) ** 2
        for j in range(1, n - 1):
            self.assertAlmostEqual(actual[j], amplification * initial[j], delta=1e-14)

    def test_heat_maximum_principle_stable_step(self):
        initial = [0, 0, 1, 0, 0]
        states = heat_1d_dirichlet(initial, diffusivity=1, dx=.25, dt=.03125, steps=20)
        self.assertTrue(all(0 <= x <= 1 for state in states for x in state))

    def test_heat_rejects_unstable_and_bad_inputs(self):
        kwargs = dict(diffusivity=1, dx=.1, dt=.006, steps=1)
        with self.assertRaises(ValueError):
            heat_1d_dirichlet([0, 1, 0], **kwargs)
        for bad in ([0, 1], [1, 1, 0], [0, math.inf, 0], "010"):
            with self.assertRaises((TypeError, ValueError)):
                heat_1d_dirichlet(bad, diffusivity=1, dx=.1, dt=.002, steps=1)
        with self.assertRaises(ValueError):
            heat_1d_dirichlet([0, 1, 0], diffusivity=-1, dx=.1, dt=.002, steps=1)
        with self.assertRaises(ValueError):
            heat_1d_dirichlet([0, 1, 0], diffusivity=1, dx=0, dt=.002, steps=1)
        with self.assertRaises(ValueError):
            heat_1d_dirichlet([0, 1, 0], diffusivity=1, dx=.1, dt=0, steps=1)
        with self.assertRaises(TypeError):
            heat_1d_dirichlet([0, 1, 0], diffusivity=1, dx=.1, dt=.002, steps=True)
        with self.assertRaises(TypeError):
            heat_1d_dirichlet((x for x in [0, 1, 0]), diffusivity=1, dx=.1,
                              dt=.002, steps=1)
        with self.assertRaises(ValueError):
            heat_1d_dirichlet([0, 1, 0], diffusivity=1, dx=.1,
                              dt=.002, steps=2_000_000)

    def test_heat_zero_diffusivity_preserves_profile(self):
        initial = [0, 3, -2, 7, 0]
        states = heat_1d_dirichlet(initial, diffusivity=0, dx=.25, dt=.9, steps=8)
        self.assertEqual(states[-1], tuple(float(v) for v in initial))

    def test_heat_zero_steps_and_exact_stability_boundary(self):
        initial = [0, 1, 0]
        self.assertEqual(heat_1d_dirichlet(initial, diffusivity=1, dx=.1, dt=.005,
                                            steps=0), ((0., 1., 0.),))
        states = heat_1d_dirichlet(initial, diffusivity=1, dx=.1, dt=.005, steps=1)
        # At r=1/2 the interior becomes the average of both neighbors.
        self.assertAlmostEqual(states[-1][1], 0., delta=1e-14)


class WaveEquationTests(unittest.TestCase):
    def test_zero_wave(self):
        states = wave_1d_dirichlet([0] * 7, [0] * 7, wave_speed=1, dx=.25, dt=.1, steps=10)
        self.assertTrue(all(state == (0.,) * 7 for state in states))

    def test_linear_steady_state(self):
        initial = [1 + .1 * j for j in range(11)]
        states = wave_1d_dirichlet(initial, [0] * 11, wave_speed=2, dx=.1, dt=.025, steps=50,
                                   left=1, right=2)
        for actual in states[-1]:
            self.assertTrue(math.isfinite(actual))
        for actual, expected in zip(states[-1], initial):
            self.assertAlmostEqual(actual, expected, delta=1e-11)

    def test_sine_mode_against_exact_wave(self):
        # For u_tt = c^2 u_xx, zero initial velocity, zero ends,
        # u(x,t)=sin(pi*x) cos(c*pi*t) on [0,1].
        n, dx, dt, steps, c = 101, .01, .003, 100, 2.
        initial = [math.sin(math.pi * j / (n - 1)) for j in range(n)]
        actual = wave_1d_dirichlet(initial, [0] * n, wave_speed=c, dx=dx, dt=dt, steps=steps)
        self.assertAlmostEqual(actual[-1][50], math.cos(c * math.pi * dt * steps), delta=1e-4)

    def test_wave_at_zero_speed_has_linear_time_displacement(self):
        initial = [0, 1, 0]
        velocity = [0, -2, 0]
        states = wave_1d_dirichlet(initial, velocity, wave_speed=0, dx=.2, dt=.1, steps=4)
        self.assertAlmostEqual(states[-1][1], .2, delta=1e-12)

    def test_wave_rejects_cfl_boundary_velocity_and_shapes(self):
        with self.assertRaises(ValueError):
            wave_1d_dirichlet([0, 1, 0], [0, 0, 0], wave_speed=1, dx=.1, dt=.11, steps=2)
        with self.assertRaises(ValueError):
            wave_1d_dirichlet([0, 1, 0], [1, 0, 0], wave_speed=1, dx=.1, dt=.05, steps=2)
        with self.assertRaises(ValueError):
            wave_1d_dirichlet([0, 1, 0], [0, 0], wave_speed=1, dx=.1, dt=.05, steps=2)
        with self.assertRaises(TypeError):
            wave_1d_dirichlet([0, 1, 0], [0, 0, 0], wave_speed=1, dx=.1, dt=.05, steps=False)
        with self.assertRaises(TypeError):
            wave_1d_dirichlet([0, 1, 0], (x for x in [0, 0, 0]), wave_speed=1,
                              dx=.1, dt=.05, steps=2)

    def test_wave_zero_steps_and_exact_courant_boundary(self):
        self.assertEqual(wave_1d_dirichlet([0, 1, 0], [0, 0, 0],
                                              wave_speed=1, dx=.1, dt=.1, steps=0),
                         ((0., 1., 0.),))
        states = wave_1d_dirichlet([0, 1, 0], [0, 0, 0],
                                    wave_speed=1, dx=.1, dt=.1, steps=2)
        # With c*dt/dx = 1, the centered stencil is exact at the grid nodes
        # for this one-step propagation pattern.
        self.assertAlmostEqual(states[1][1], 0., delta=1e-14)
        self.assertAlmostEqual(states[2][1], -1., delta=1e-14)


if __name__ == "__main__":
    unittest.main()
