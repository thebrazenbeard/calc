"""Exact independent Bellman, Riccati, and open-loop LQR witnesses."""

import random
import unittest
from fractions import Fraction as Q
from itertools import product

from calc.control import (
    ScalarLQR, ScalarTrajectory, evaluate_trajectory, optimal_trajectory,
    solve_scalar_lqr,
)


class ScalarControlTests(unittest.TestCase):
    def test_one_stage_analytic_completion_of_square(self):
        solution = solve_scalar_lqr(1, 1, 1, 1, 1, 1)
        self.assertEqual(solution.cost_to_go, (Q(3, 2), Q(1)))
        self.assertEqual(solution.gains, (Q(1, 2),))
        self.assertEqual(solution.control(2, 0), -1)
        self.assertEqual(solution.value(2), 6)
        self.assertEqual(optimal_trajectory(solution, 2),
                         ScalarTrajectory((Q(2), Q(1)), (Q(-1),), Q(6)))

    def test_two_stage_closed_form_feedback(self):
        solution = solve_scalar_lqr(1, 1, 1, 1, 1, 2)
        self.assertEqual(solution.cost_to_go, (Q(8, 5), Q(3, 2), Q(1)))
        self.assertEqual(solution.gains, (Q(3, 5), Q(1, 2)))
        trajectory = optimal_trajectory(solution, 5)
        self.assertEqual(trajectory.states, (Q(5), Q(2), Q(1)))
        self.assertEqual(trajectory.controls, (Q(-3), Q(-1)))
        self.assertEqual(trajectory.total_cost, Q(40))

    def test_zero_horizon_is_terminal_cost_only(self):
        solution = solve_scalar_lqr(2, 1, 7, 3, Q(5, 2), 0)
        self.assertEqual(solution.cost_to_go, (Q(5, 2),))
        self.assertEqual(solution.gains, ())
        self.assertEqual(optimal_trajectory(solution, 3),
                         ScalarTrajectory((Q(3),), (), Q(45, 2)))
        self.assertEqual(evaluate_trajectory(solution, 3, ()),
                         optimal_trajectory(solution, 3))

    def test_no_input_control_uses_zero_optimal_effort(self):
        solution = solve_scalar_lqr(2, 0, 3, 7, 5, 3)
        self.assertEqual(solution.gains, (Q(0),) * 3)
        self.assertEqual(optimal_trajectory(solution, 1).states,
                         (Q(1), Q(2), Q(4), Q(8)))
        self.assertEqual(optimal_trajectory(solution, 1).controls, (Q(0),) * 3)

    def test_zero_state_weight_with_no_terminal_penalty(self):
        solution = solve_scalar_lqr(Q(3, 2), Q(2), 0, Q(1, 2), 0, 10)
        self.assertEqual(solution.cost_to_go, (Q(0),) * 11)
        self.assertEqual(solution.gains, (Q(0),) * 10)
        self.assertEqual(optimal_trajectory(solution, 0).total_cost, 0)

    def test_exact_scalar_scaling_of_value_and_policy(self):
        solution = solve_scalar_lqr(Q(5, 3), Q(-2, 3), Q(3, 5), Q(2), Q(7, 3), 4)
        for scalar in (Q(-5, 2), Q(0), Q(2, 3), Q(7)):
            self.assertEqual(solution.value(scalar, 2),
                             scalar * scalar * solution.value(1, 2))
            self.assertEqual(solution.control(scalar, 1),
                             scalar * solution.control(1, 1))

    def test_direct_open_loop_cost_oracle(self):
        solution = solve_scalar_lqr(Q(2), Q(-1), Q(1, 3), Q(4), Q(5), 3)
        controls = (Q(1, 2), Q(-2), Q(3, 4))
        x = Q(3)
        states = [x]
        cost = Q(0)
        for u in controls:
            cost += x * x / 3 + 4 * u * u
            x = 2 * x - u
            states.append(x)
        cost += 5 * x * x
        trajectory = evaluate_trajectory(solution, 3, controls)
        self.assertEqual(trajectory.states, tuple(states))
        self.assertEqual(trajectory.total_cost, cost)

    def test_exhaustive_small_rational_open_loop_grid(self):
        # Candidate grid deliberately includes controls that are not optimal;
        # exact rational comparisons certify all sampled candidates cost >= V.
        solution = solve_scalar_lqr(1, 1, 1, 1, 1, 2)
        baseline = optimal_trajectory(solution, 1).total_cost
        candidates = [Q(k, 4) for k in range(-8, 9)]
        for controls in product(candidates, repeat=2):
            self.assertGreaterEqual(
                evaluate_trajectory(solution, 1, controls).total_cost, baseline)

    def test_independent_bellman_completed_square_identity_seeded(self):
        rnd = random.Random(20261101)
        for _ in range(80):
            a = Q(rnd.randint(-3, 3), rnd.randint(1, 4))
            b = Q(rnd.randint(-3, 3), rnd.randint(1, 4))
            q = Q(rnd.randint(0, 4), rnd.randint(1, 4))
            r = Q(rnd.randint(1, 5), rnd.randint(1, 4))
            terminal = Q(rnd.randint(0, 4), rnd.randint(1, 4))
            horizon = rnd.randint(1, 5)
            solution = solve_scalar_lqr(a, b, q, r, terminal, horizon)
            for stage in range(horizon):
                x = Q(rnd.randint(-4, 4), 3)
                u = Q(rnd.randint(-4, 4), 2)
                next_p = solution.cost_to_go[stage + 1]
                # This identity follows by hand from completing the square:
                # cost(u)-V_k(x) = (r+b²P_(k+1))*(u+K_k*x)^2.
                lhs = q * x * x + r * u * u + next_p * (a * x + b * u) ** 2
                lhs -= solution.cost_to_go[stage] * x * x
                rhs = (r + b * b * next_p) * (u + solution.gains[stage] * x) ** 2
                self.assertEqual(lhs, rhs)
            trajectory = optimal_trajectory(solution, Q(2, 3))
            self.assertEqual(trajectory.total_cost, solution.value(Q(2, 3)))

    def test_all_cost_coefficients_nonnegative(self):
        solution = solve_scalar_lqr(Q(-7, 5), Q(3, 2), Q(1, 3), Q(2, 5), 4, 20)
        self.assertTrue(all(v >= 0 for v in solution.cost_to_go))

    def test_stability_not_assumed_or_forced(self):
        solution = solve_scalar_lqr(3, 0, 1, 1, 1, 3)
        self.assertEqual(optimal_trajectory(solution, 1).states[-1], Q(27))

    def test_negative_weights_and_control_convexity_rejected(self):
        for kwargs in (
            dict(q=-1, r=1, terminal_q=0),
            dict(q=1, r=0, terminal_q=0),
            dict(q=1, r=-2, terminal_q=0),
            dict(q=1, r=1, terminal_q=-1),
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                solve_scalar_lqr(1, 1, horizon=1, **kwargs)

    def test_floats_booleans_text_and_huge_inputs_rejected(self):
        for invalid in (True, 1.0, 0.5, None, "1", complex(1, 0), 1 << 5000):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                solve_scalar_lqr(invalid, 1, 1, 1, 1, 1)
        with self.assertRaises(TypeError):
            solve_scalar_lqr(1, 1, 1, 1, 1, True)

    def test_horizon_bounds(self):
        for invalid in (-1, 65, Q(2), "3"):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                solve_scalar_lqr(1, 1, 1, 1, 1, invalid)

    def test_state_and_stage_bounds(self):
        solution = solve_scalar_lqr(1, 1, 1, 1, 1, 2)
        for invalid in (True, 1.5, 1 << 5000, "2"):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                solution.value(invalid)
        for invalid in (-1, 3, True, Q(1)):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                solution.value(1, invalid)
        with self.assertRaises(ValueError):
            solution.control(1, 2)

    def test_control_sequence_length_and_type_bounded(self):
        solution = solve_scalar_lqr(1, 1, 1, 1, 1, 3)
        for invalid in ((), (1, 2), (1, 2, 3, 4), (1, 2, 0.5), "123"):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                evaluate_trajectory(solution, 0, invalid)
        with self.assertRaises(TypeError):
            evaluate_trajectory(solution, 0, None)
        with self.assertRaises(ValueError):
            evaluate_trajectory(solution, 0, (Q(0) for _ in iter(int, 1)))

    def test_trajectory_bit_growth_budget(self):
        solution = solve_scalar_lqr(2, 1, 1, 1, 1, 64)
        with self.assertRaises(ValueError):
            evaluate_trajectory(solution, Q(1, 2 ** 4095), (Q(1) for _ in range(64)))


if __name__ == "__main__":
    unittest.main()
