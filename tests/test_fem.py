"""Independent weak-form, exactness, convergence and negative FEM witnesses."""

import unittest
from fractions import Fraction as Q

from calc.fem import assemble_poisson_1d, solve_poisson_1d


class FiniteElementTests(unittest.TestCase):
    def test_two_elements_exact_stiffness_and_load(self):
        system = assemble_poisson_1d((0, Q(1, 2), 1))
        self.assertEqual(system.diagonal, (Q(4),))
        self.assertEqual(system.upper, ())
        self.assertEqual(system.load, (Q(1, 2),))
        self.assertEqual(system.forcing, (Q(1), Q(1)))

    def test_nonuniform_three_element_assembly(self):
        system = assemble_poisson_1d((0, Q(1, 4), Q(1, 2), 1))
        self.assertEqual(system.diagonal, (Q(8), Q(6)))
        self.assertEqual(system.upper, (Q(-4),))
        self.assertEqual(system.load, (Q(1, 4), Q(3, 8)))

    def test_poisson_quadratic_nodal_exactness_nonuniform(self):
        nodes = (Q(-1), Q(-3, 5), Q(-1, 7), Q(2, 5), Q(2))
        solution = solve_poisson_1d(nodes, 3)
        a, b = nodes[0], nodes[-1]
        expected = tuple(Q(3, 2) * (x - a) * (b - x) for x in nodes)
        self.assertEqual(solution.nodal_values, expected)
        self.assertEqual(solution.system.residual(solution.nodal_values[1:-1]),
                         (Q(0),) * (len(nodes) - 2))

    def test_piecewise_constant_load_assembly(self):
        system = assemble_poisson_1d((0, Q(1, 2), 1), (2, 0))
        self.assertEqual(system.forcing, (Q(2), Q(0)))
        self.assertEqual(system.load, (Q(1, 2),))
        self.assertEqual(system.solve().nodal_values, (Q(0), Q(1, 8), Q(0)))

    def test_piecewise_load_nodal_solution_residual(self):
        system = assemble_poisson_1d((0, Q(1, 5), Q(1, 2), 1), (2, -1, 4))
        solution = system.solve()
        self.assertEqual(system.residual(solution.nodal_values[1:-1]), (Q(0), Q(0)))
        self.assertTrue(all(type(v) is Q for v in solution.nodal_values))

    def test_galerkin_equation_against_independent_test_vector(self):
        system = assemble_poisson_1d((0, Q(1, 4), Q(1, 2), Q(3, 4), 1),
                                      (2, 1, 0, 3))
        sol = system.solve().nodal_values[1:-1]
        witness = (Q(2), Q(-3), Q(5))
        self.assertEqual(sum(v * r for v, r in zip(witness, system.residual(sol))), 0)
        # Independent assembly of bilinear energy ∑ (Δu)(Δv)/h.
        full_u, full_v = (Q(0), *sol, Q(0)), (Q(0), *witness, Q(0))
        left = sum((full_u[i + 1] - full_u[i]) *
                   (full_v[i + 1] - full_v[i]) /
                   (system.nodes[i + 1] - system.nodes[i])
                   for i in range(len(system.nodes) - 1))
        right = sum(w * b for w, b in zip(witness, system.load))
        self.assertEqual(left, right)

    def test_continuous_piecewise_affine_interpolation(self):
        sol = solve_poisson_1d((0, Q(1, 2), 1))
        self.assertEqual(sol.evaluate(0), 0)
        self.assertEqual(sol.evaluate(1), 0)
        self.assertEqual(sol.evaluate(Q(1, 2)), Q(1, 8))
        self.assertEqual(sol.evaluate(Q(1, 4)), Q(1, 16))
        self.assertEqual(sol.evaluate(Q(3, 4)), Q(1, 16))

    def test_interpolation_error_quadratic_convergence(self):
        errors = []
        for intervals in (2, 4, 8, 16):
            nodes = tuple(Q(i, intervals) for i in range(intervals + 1))
            sol = solve_poisson_1d(nodes)
            # Concave exact solution x(1-x)/2 is above every affine chord.
            mid = Q(1, 2 * intervals)
            exact = mid * (1 - mid) / 2
            error = exact - sol.evaluate(mid)
            self.assertEqual(error, Q(1, 8 * intervals * intervals))
            errors.append(error)
        self.assertEqual([errors[i] / errors[i + 1] for i in range(3)], [4] * 3)

    def test_energy_seminorm_exact(self):
        sol = solve_poisson_1d((0, Q(1, 2), 1))
        self.assertEqual(sol.energy_seminorm_squared(), Q(1, 16))

    def test_scaling_and_sign_symmetry(self):
        nodes = (Q(-2), Q(-1), Q(0), Q(2))
        plus = solve_poisson_1d(nodes, 7).nodal_values
        minus = solve_poisson_1d(nodes, -7).nodal_values
        base = solve_poisson_1d(nodes, 1).nodal_values
        self.assertEqual(plus, tuple(7 * v for v in base))
        self.assertEqual(minus, tuple(-v for v in plus))

    def test_zero_load_and_single_element(self):
        z = solve_poisson_1d((0, 1, 2, 3), 0)
        self.assertEqual(z.nodal_values, (Q(0),) * 4)
        self.assertEqual(z.energy_seminorm_squared(), 0)
        only = assemble_poisson_1d((3, 5), 42)
        self.assertEqual(only.degrees_of_freedom, 0)
        self.assertEqual(only.residual(()), ())
        self.assertEqual(only.solve().nodal_values, (Q(0), Q(0)))

    def test_rejects_invalid_meshes(self):
        bad = ((), (0,), (0, 0, 1), (2, 1), (0, 1, 1),
               (0, float('inf')), (0, .5, 1), (False, 1), (0, '1'))
        for nodes in bad:
            with self.subTest(nodes=nodes), self.assertRaises((ValueError, TypeError)):
                assemble_poisson_1d(nodes)

    def test_rejects_invalid_forcing(self):
        mesh = (0, Q(1, 2), 1)
        for f in ((1,), (1, 2, 3), (1, 1.5), True, .5, "1", (1, None)):
            with self.subTest(forcing=f), self.assertRaises((ValueError, TypeError)):
                assemble_poisson_1d(mesh, f)

    def test_invalid_residual_and_evaluation(self):
        system = assemble_poisson_1d((0, Q(1, 2), 1))
        for values in ((), (1, 2), (True,), (.3,)):
            with self.subTest(values=values), self.assertRaises((ValueError, TypeError)):
                system.residual(values)
        sol = system.solve()
        for point in (-1, Q(11, 10), .5, True, '1'):
            with self.subTest(point=point), self.assertRaises((ValueError, TypeError)):
                sol.evaluate(point)

    def test_rejects_mesh_node_and_rational_size_over_budget(self):
        with self.assertRaises(ValueError):
            assemble_poisson_1d(range(514))
        with self.assertRaises(ValueError):
            assemble_poisson_1d((0, 1 << 5000))
        with self.assertRaises(ValueError):
            assemble_poisson_1d((0, 1), 1 << 5000)

    def test_immutable_solution_and_system(self):
        sol = solve_poisson_1d((0, Q(1, 2), 1))
        with self.assertRaises(AttributeError):
            sol.nodal_values = (Q(0),)
        with self.assertRaises(AttributeError):
            sol.system.diagonal = ()


if __name__ == '__main__':
    unittest.main()
