"""Adversarial independent witnesses for the one-dimensional PDE solvers.

Tests derive expected results from the discretized PDE, stability conditions,
and continuum analytic sine modes. They do not mirror implementation helpers.
"""

import math
import unittest

from calc.pde import heat_1d_dirichlet, wave_1d_dirichlet


class HostileHeatTests(unittest.TestCase):
    def test_cfl_half_pulse_avoids_false_intermediate_overflow(self):
        # With r=1/2 the FTCS update is the mean of the two neighbours.
        # -2*1e308 overflows, but the mathematically correct result is zero.
        rows = heat_1d_dirichlet(
            [0.0, 1e308, 0.0], diffusivity=1, dx=1, dt=0.5, steps=1
        )
        self.assertEqual(rows[-1], (0.0, 0.0, 0.0))

    def test_maximum_principle_for_signed_profiles(self):
        initial = (3.0, -8.0, 12.0, -4.0, 2.0)
        rows = heat_1d_dirichlet(
            initial, diffusivity=1, dx=1, dt=0.4, steps=40, left=3, right=2
        )
        self.assertTrue(all(min(initial) - 1e-12 <= v <= max(initial) + 1e-12
                            for row in rows for v in row))
        self.assertTrue(all(row[0] == 3 and row[-1] == 2 for row in rows))

    def test_units_scaling_preserves_discrete_heat_history(self):
        profile = (0, 1, -3, 2, 0)
        a = heat_1d_dirichlet(profile, diffusivity=2, dx=1, dt=.2, steps=4)
        b = heat_1d_dirichlet(profile, diffusivity=2, dx=2, dt=.8, steps=4)
        self.assertEqual(a, b)

    def test_sine_mode_spatial_refinement_is_second_order(self):
        errors = []
        for n in (11, 21, 41):
            dx = 1 / (n - 1)
            dt = .4 * dx * dx
            steps = round(.04 / dt)
            profile = [math.sin(math.pi * j * dx) for j in range(n)]
            actual = heat_1d_dirichlet(profile, diffusivity=1, dx=dx,
                                       dt=dt, steps=steps)[-1][n // 2]
            exact = math.exp(-math.pi * math.pi * .04)
            errors.append(abs(actual - exact))
        self.assertLess(errors[1], errors[0] / 3)
        self.assertLess(errors[2], errors[1] / 3)

    def test_immutable_history_and_input(self):
        initial = [0, 1, 0]
        history = heat_1d_dirichlet(initial, diffusivity=1, dx=1, dt=.25, steps=2)
        self.assertEqual(initial, [0, 1, 0])
        self.assertEqual(history[0], (0., 1., 0.))
        self.assertNotEqual(history[0], history[1])
        with self.assertRaises(TypeError):
            history[0][1] = 42


class HostileWaveTests(unittest.TestCase):
    def test_courant_one_pulse_avoids_false_intermediate_overflow(self):
        # For Courant number one and zero initial velocity, center advances
        # 1e308 -> 0 -> -1e308, all representable finite values.
        rows = wave_1d_dirichlet([0.0, 1e308, 0.0], [0.0] * 3,
                                 wave_speed=1, dx=1, dt=1, steps=2)
        self.assertEqual(tuple(row[1] for row in rows), (1e308, 0.0, -1e308))

    def test_units_scaling_preserves_discrete_wave_history(self):
        initial = (0, 1, -3, 2, 0)
        velocity = (0, 0, 0, 0, 0)
        a = wave_1d_dirichlet(initial, velocity, wave_speed=2, dx=1,
                              dt=.25, steps=4)
        b = wave_1d_dirichlet(initial, velocity, wave_speed=2, dx=2,
                              dt=.5, steps=4)
        self.assertEqual(a, b)

    def test_initial_velocity_matches_sine_wave_solution(self):
        n = 101
        dx, dt, steps = .01, .005, 20
        velocity = [math.sin(math.pi * j * dx) for j in range(n)]
        rows = wave_1d_dirichlet([0.] * n, velocity, wave_speed=1,
                                 dx=dx, dt=dt, steps=steps)
        expected = math.sin(math.pi * .1) / math.pi
        self.assertAlmostEqual(rows[-1][50], expected, delta=1e-5)

    def test_zero_speed_and_fixed_boundary(self):
        initial = [2, 5, 8]
        velocity = [0, -3, 0]
        rows = wave_1d_dirichlet(initial, velocity, wave_speed=0,
                                 dx=1, dt=.25, steps=4, left=2, right=8)
        self.assertEqual(rows[-1], (2., 2., 8.))
        self.assertTrue(all(row[0] == 2. and row[-1] == 8. for row in rows))

    def test_invalid_signed_speed_and_zero_dt(self):
        with self.assertRaises(ValueError):
            wave_1d_dirichlet([0, 1, 0], [0, 0, 0], wave_speed=-1,
                              dx=1, dt=.1, steps=1)
        with self.assertRaises(ValueError):
            wave_1d_dirichlet([0, 1, 0], [0, 0, 0], wave_speed=1,
                              dx=1, dt=0, steps=1)


if __name__ == "__main__":
    unittest.main()
