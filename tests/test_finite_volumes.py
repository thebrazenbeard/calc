"""Independent finite-volume witnesses for periodic scalar transport."""

import math
import unittest
from fractions import Fraction

from calc.finite_volumes import (
    advance_upwind,
    cell_mass,
    periodic_total_variation,
    upwind_fluxes,
    upwind_step,
)


class FiniteVolumeTests(unittest.TestCase):
    def test_positive_velocity_donor_flux(self):
        self.assertEqual(upwind_fluxes((1, 2, 3, 4), 2), (2, 4, 6, 8))

    def test_negative_velocity_donor_flux(self):
        self.assertEqual(upwind_fluxes((1, 2, 3, 4), -2), (-4, -6, -8, -2))

    def test_zero_velocity_gives_zero_flux(self):
        self.assertEqual(upwind_fluxes((1, -2, 3), 0), (0, 0, 0))

    def test_positive_velocity_one_cell_translation_at_cfl_one(self):
        result = upwind_step((1, 2, 3, 4), dx=2, velocity=2, dt=1)
        self.assertEqual(result, (4, 1, 2, 3))

    def test_negative_velocity_one_cell_translation_at_cfl_one(self):
        result = upwind_step((1, 2, 3, 4), dx=2, velocity=-2, dt=1)
        self.assertEqual(result, (2, 3, 4, 1))

    def test_repeated_exact_translation_and_periodicity(self):
        state = (0.0, 1.0, 3.0, 2.0, -1.0)
        self.assertEqual(advance_upwind(state, dx=0.2, velocity=1.0, dt=0.2,
                                        steps=5), state)

    def test_zero_timestep_and_zero_velocity_are_identity(self):
        state = (2.0, 3.0, -7.0)
        self.assertEqual(upwind_step(state, dx=0.5, velocity=7, dt=0), state)
        self.assertEqual(upwind_step(state, dx=0.5, velocity=0, dt=1), state)
        self.assertEqual(advance_upwind(state, dx=0.5, velocity=1, dt=0.2,
                                        steps=0), state)

    def test_conservative_flux_difference_identity(self):
        state = (1.0, 2.0, -1.0, 4.0)
        flux = upwind_fluxes(state, -0.75)
        new = upwind_step(state, dx=0.5, velocity=-0.75, dt=0.25)
        for i, old in enumerate(state):
            expected = old - 0.25 / 0.5 * (flux[i] - flux[i - 1])
            self.assertAlmostEqual(new[i], expected)

    def test_constant_state_is_preserved(self):
        for velocity in (-2, 0, 2):
            state = (3.125,) * 13
            result = advance_upwind(state, dx=0.1, velocity=velocity,
                                    dt=0.04, steps=17)
            self.assertEqual(result, state)

    def test_mass_conservation_over_many_steps(self):
        state = (1.0, -2.0, 4.0, 0.25, -5.0, 3.5, 1.25)
        before = cell_mass(state, dx=0.125)
        for velocity in (-2.0, 2.0):
            after = advance_upwind(state, dx=0.125, velocity=velocity,
                                   dt=0.05, steps=300)
            self.assertAlmostEqual(cell_mass(after, dx=0.125), before, places=11)

    def test_maximum_principle_under_cfl(self):
        state = (-5.0, 10.0, -2.0, 7.0, 0.0, 3.0)
        low, high = min(state), max(state)
        for velocity in (-1.0, 1.0):
            result = advance_upwind(state, dx=0.2, velocity=velocity,
                                    dt=0.07, steps=100)
            self.assertTrue(all(low <= value <= high for value in result))

    def test_discrete_total_variation_does_not_increase(self):
        state = (0, 3, -4, 5, -2, 0, 8, -1)
        initial = periodic_total_variation(state)
        for velocity in (-1.0, 1.0):
            result = advance_upwind(state, dx=0.25, velocity=velocity,
                                    dt=0.125, steps=20)
            self.assertLessEqual(periodic_total_variation(result), initial + 1e-12)

    def test_total_variation_concrete_example(self):
        self.assertEqual(periodic_total_variation((0, 1, 1, 0)), 2)
        self.assertEqual(periodic_total_variation((2, 2)), 0)

    def test_first_order_convergence_on_periodic_sine(self):
        def error(n):
            dx = 1 / n
            mean_sine_scale = math.sin(math.pi / n) / (math.pi / n)
            samples = tuple(mean_sine_scale * math.sin(2 * math.pi * ((i + 0.5) / n))
                            for i in range(n))
            # With velocity one and nu=0.5, n/2 steps reach t=1/4.
            result = advance_upwind(samples, dx=dx, velocity=1.0,
                                    dt=0.5 * dx, steps=n // 2)
            exact = tuple(mean_sine_scale * math.sin(2 * math.pi * ((i + 0.5) / n - 0.25))
                          for i in range(n))
            return math.fsum(abs(a-b) for a, b in zip(result, exact)) / n

        coarse = error(40)
        fine = error(80)
        self.assertGreater(coarse, 0.005)
        self.assertLess(fine, 0.65 * coarse)

    def test_exact_binary64_cfl_handles_small_spacing(self):
        self.assertEqual(upwind_step((1.0, 2.0), dx=1e-300,
                                     velocity=1e-300, dt=1.0), (2.0, 1.0))

    def test_mass_recovers_finite_cancellation_after_term_overflow(self):
        self.assertEqual(cell_mass((1e308, -1e308), dx=2.0), 0.0)

    def test_fraction_inputs_are_permitted_but_float_results_documented(self):
        state = (Fraction(1, 2), Fraction(3, 2))
        result = upwind_step(state, dx=Fraction(1, 2),
                             velocity=Fraction(1, 4), dt=1)
        self.assertEqual(result, (1.0, 1.0))

    def test_cfl_violation_rejected(self):
        with self.assertRaisesRegex(ValueError, "CFL"):
            upwind_step((1, 2, 3), dx=0.1, velocity=1.0, dt=0.1000001)

    def test_invalid_grid_and_state_rejected(self):
        for state in ((1,), (), (1, float("nan")), (1, float("inf")), (1, True)):
            with self.subTest(state=state), self.assertRaises((TypeError, ValueError)):
                upwind_step(state, dx=1, velocity=1, dt=0.5)
        with self.assertRaises(TypeError):
            upwind_step(iter((1, 2)), dx=1, velocity=1, dt=0.5)
        for dx in (0, -1, float("inf"), float("nan"), True):
            with self.subTest(dx=dx), self.assertRaises((TypeError, ValueError)):
                upwind_step((1, 2), dx=dx, velocity=1, dt=0.5)

    def test_invalid_time_and_velocity_rejected(self):
        for bad in (float("inf"), float("nan"), True, "1"):
            with self.subTest(bad=bad), self.assertRaises((TypeError, ValueError)):
                upwind_step((1, 2), dx=1, velocity=bad, dt=0.5)
            with self.subTest(bad=bad), self.assertRaises((TypeError, ValueError)):
                upwind_step((1, 2), dx=1, velocity=1, dt=bad)
        with self.assertRaises(ValueError):
            upwind_step((1, 2), dx=1, velocity=1, dt=-0.5)

    def test_bound_steps_cells_and_workload(self):
        with self.assertRaises((TypeError, ValueError)):
            advance_upwind((1, 2), dx=1, velocity=1, dt=0.5, steps=True)
        for bad in (-1, 100_001):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                advance_upwind((1, 2), dx=1, velocity=1, dt=0.5, steps=bad)
        with self.assertRaises(ValueError):
            advance_upwind([0.0] * 1000, dx=1, velocity=1, dt=0.5,
                           steps=2001)

    def test_true_unrepresentable_magnitudes_fail_cleanly(self):
        with self.assertRaises(ValueError):
            cell_mass((1e308, 1e308), dx=2)
        with self.assertRaises(ValueError):
            upwind_fluxes((1e308, 2.0), 2)
        with self.assertRaises(ValueError):
            periodic_total_variation((1e308, -1e308))


if __name__ == "__main__":
    unittest.main()
