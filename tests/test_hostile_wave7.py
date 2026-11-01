"""Independent adversarial mathematical witnesses for Wave 7.

Oracles come from algebraic identities, weak forms, and exact closed forms;
they intentionally do not import internal helper functions from the targets.
"""

from fractions import Fraction
import itertools
import random
import unittest

from calc.clifford import CliffordAlgebra
from calc.distributions import DiracMeasure, PiecewiseConstant
from calc.fem import solve_poisson_1d
from calc.finite_volumes import advance_upwind, cell_mass, periodic_total_variation, upwind_step
from calc.lie import (adjoint_matrix, commutator, jacobiator, linear_combination,
                      sl2_basis, so3_basis, structure_constants)


class HostileCliffordTests(unittest.TestCase):
    def test_independent_generator_reordering_oracle(self):
        # In the free associative algebra, move right-side generators to
        # canonical order by anticommutation, then replace repeated e_i e_i
        # with the specified metric sign. This derives each expected blade
        # product without calling the implementation's bitmask helper.
        for metric in ((1, 1, 1, 1), (-1, -1, 1, -1), (1, -1, 1, -1)):
            algebra = CliffordAlgebra(metric)
            for left in range(16):
                for right in range(16):
                    generators = [i for i in range(4) if left & (1 << i)]
                    generators.extend(i for i in range(4) if right & (1 << i))
                    swaps = sum(generators[i] > generators[j]
                                for i in range(len(generators))
                                for j in range(i + 1, len(generators)))
                    sign = -1 if swaps % 2 else 1
                    for i in range(4):
                        if generators.count(i) == 2:
                            sign *= metric[i]
                    expected_mask = left ^ right
                    self.assertEqual(algebra.basis(left) * algebra.basis(right),
                                     algebra.multivector({expected_mask: sign}))

    def test_reversion_quadratic_indefinite_sign_and_sparse_exactness(self):
        algebra = CliffordAlgebra((1, -1, 1, -1))
        coefficients = {0: Fraction(1, 2), 3: Fraction(2, 5),
                        6: Fraction(-7, 3), 8: Fraction(4, 7)}
        x = algebra.multivector(coefficients)
        # For a diagonal orthogonal metric, the scalar part of ~X X is
        # sum_B coefficient[B]^2 * product_{i in B} g_ii.
        oracle = Fraction(0)
        for mask, coefficient in coefficients.items():
            metric_sign = 1
            for i, metric in enumerate(algebra.signature):
                if mask & (1 << i):
                    metric_sign *= metric
            oracle += coefficient * coefficient * metric_sign
        self.assertEqual(x.reversion_quadratic(), oracle)
        self.assertEqual((x * x.algebra.scalar(3)).as_dict(),
                         {mask: 3 * c for mask, c in coefficients.items()})


class HostileLieTests(unittest.TestCase):
    def test_sl2_structure_constants_against_exact_commutation_table(self):
        h, e, f = sl2_basis()
        basis = (h, e, f)
        constants = structure_constants(basis)
        expected = {
            (0, 1): (0, 2, 0), (1, 0): (0, -2, 0),
            (0, 2): (0, 0, -2), (2, 0): (0, 0, 2),
            (1, 2): (1, 0, 0), (2, 1): (-1, 0, 0),
        }
        for i, j in itertools.product(range(3), repeat=2):
            self.assertEqual(constants[i][j], tuple(map(Fraction, expected.get((i, j), (0, 0, 0)))))
            self.assertEqual(commutator(basis[i], basis[j]),
                             linear_combination(constants[i][j], basis))

    def test_adjoint_is_lie_representation(self):
        basis = sl2_basis()
        x = linear_combination((Fraction(2, 3), -1, 4), basis)
        y = linear_combination((1, Fraction(5, 7), -3), basis)
        ax = adjoint_matrix(x, basis)
        ay = adjoint_matrix(y, basis)
        axy = adjoint_matrix(commutator(x, y), basis)
        self.assertEqual(axy, commutator(ax, ay))
        self.assertEqual(jacobiator(x, y, basis[0]), ((0, 0), (0, 0)))

    def test_so3_structure_bracket_is_cyclic(self):
        a, b, c = so3_basis()
        self.assertEqual(commutator(a, b), c)
        self.assertEqual(commutator(b, c), a)
        self.assertEqual(commutator(c, a), b)
        self.assertEqual(structure_constants((a, b, c))[0][1], (0, 0, 1))


class HostileFEMTests(unittest.TestCase):
    def test_nonuniform_poisson_nodes_match_closed_form_green_solution(self):
        # -u'' = q, u(a)=u(b)=0 => u(x)=q*(x-a)*(b-x)/2.
        # P1 FEM with exact constant forcing reproduces its nodal values on
        # every partition in 1D, even when the interior interpolant differs.
        rng = random.Random(710231)
        for _ in range(48):
            widths = [Fraction(rng.randrange(1, 9), rng.randrange(1, 12))
                      for _ in range(rng.randrange(2, 11))]
            start = Fraction(rng.randrange(-7, 6), 5)
            nodes = [start]
            for width in widths:
                nodes.append(nodes[-1] + width)
            q = Fraction(rng.randrange(-9, 10), rng.randrange(1, 10))
            solution = solve_poisson_1d(nodes, q)
            a, b = nodes[0], nodes[-1]
            expected = tuple(q * (x-a) * (b-x) / 2 for x in nodes)
            self.assertEqual(solution.nodal_values, expected)
            self.assertEqual(solution.system.residual(expected[1:-1]),
                             (Fraction(0),) * (len(nodes) - 2))

    def test_galerkin_energy_equals_piecewise_load_pairing(self):
        nodes = (Fraction(-2), Fraction(-1, 2), Fraction(1, 3), Fraction(2))
        loads = (Fraction(2, 3), Fraction(-5, 7), Fraction(3, 4))
        solution = solve_poisson_1d(nodes, loads)
        weak_load = sum((loads[k] * (nodes[k+1] - nodes[k]) *
                         (solution.nodal_values[k] + solution.nodal_values[k+1]) / 2
                         for k in range(len(loads))), Fraction(0))
        self.assertEqual(solution.energy_seminorm_squared(), weak_load)

    def test_piecewise_constant_force_matches_green_integral_at_nodes(self):
        # Integrate -u''=f twice independently in physical coordinates.
        # With u(a)=u(b)=0 the initial derivative is
        # s0 = (1/(b-a)) integral_a^b (b-t) f(t) dt.
        rng = random.Random(94122)
        for _ in range(35):
            left = Fraction(rng.randrange(-5, 3), 3)
            nodes = [left]
            for _k in range(rng.randrange(2, 8)):
                nodes.append(nodes[-1] + Fraction(rng.randrange(1, 6),
                                                   rng.randrange(1, 9)))
            loads = [Fraction(rng.randrange(-8, 9), rng.randrange(1, 7))
                     for _k in range(len(nodes)-1)]
            right = nodes[-1]
            s0 = sum((f * (right * (x1-x0) - (x1*x1-x0*x0)/2)
                      for f, x0, x1 in zip(loads, nodes, nodes[1:])), Fraction(0))
            s0 /= right - left
            oracle = []
            for x in nodes:
                u = s0 * (x-left)
                for f, x0, x1 in zip(loads, nodes, nodes[1:]):
                    if x <= x0:
                        break
                    stop = min(x, x1)
                    u -= f * (x*(stop-x0) - (stop*stop-x0*x0)/2)
                oracle.append(u)
            self.assertEqual(solve_poisson_1d(nodes, loads).nodal_values, tuple(oracle))


class HostileFiniteVolumeTests(unittest.TestCase):
    def test_dyadic_periodic_upwind_exact_convex_oracle(self):
        # With nu=1/4 the coefficients are exactly representable binary64.
        # Compare repeated float solver steps to an independently built
        # rational convex recurrence and its exact conserved mean.
        rng = random.Random(415)
        for velocity in (-2, 2):
            for _ in range(36):
                start = tuple(Fraction(rng.randrange(-7, 8)) for _ in range(8))
                oracle = start
                for _step in range(4):
                    oracle = tuple(
                        Fraction(3, 4) * oracle[i] + Fraction(1, 4) *
                        oracle[(i-1) % len(oracle) if velocity > 0 else (i+1) % len(oracle)]
                        for i in range(len(oracle))
                    )
                result = advance_upwind(start, dx=2, velocity=velocity,
                                        dt=Fraction(1, 4), steps=4)
                self.assertEqual(result, tuple(map(float, oracle)))
                self.assertEqual(cell_mass(result, dx=2), float(2 * sum(start)))
                self.assertLessEqual(periodic_total_variation(result),
                                     periodic_total_variation(start))

    def test_maximum_principle_across_signed_speeds(self):
        values = (-7.0, 3.5, -2.5, 1.0, 9.0)
        for velocity in (-3, -1, 0, 1, 3):
            result = upwind_step(values, dx=6, velocity=velocity, dt=1.0)
            self.assertGreaterEqual(min(result), min(values))
            self.assertLessEqual(max(result), max(values))


class HostileDistributionTests(unittest.TestCase):
    def test_jump_measure_matches_independent_one_sided_values(self):
        rng = random.Random(8990)
        for _ in range(48):
            raw = []
            right = rng.randrange(-15, -9)
            for _ in range(rng.randrange(1, 8)):
                left = right + rng.randrange(0, 3)
                right = left + rng.randrange(1, 4)
                raw.append((left, right, rng.randrange(-7, 8)))
            function = PiecewiseConstant(raw)
            atoms = dict(function.weak_derivative().atoms)
            endpoints = sorted({point for l, r, _ in raw for point in (l, r)})
            for point in endpoints:
                delta = Fraction(1, 4)
                jump = function.value_at(Fraction(point) + delta) - function.value_at(Fraction(point) - delta)
                self.assertEqual(atoms.get(Fraction(point), 0), jump)
            self.assertEqual(function.weak_derivative().total_mass(), 0)

    def test_closed_form_distribution_pairings(self):
        f = PiecewiseConstant(((-3, -1, 2), (0, 2, -1)))
        self.assertEqual(f.weak_derivative().atoms,
                         ((Fraction(-3), Fraction(2)), (Fraction(-1), Fraction(-2)),
                          (Fraction(0), Fraction(-1)), (Fraction(2), Fraction(1))))
        # p(x)=x^2+2x+3, p'=2x+2, integrate exactly without using f.pair.
        p = lambda x: x*x + 2*x + 3
        oracle = 2 * p(-3) - 2 * p(-1) - p(0) + p(2)
        self.assertEqual(f.weak_derivative().pair_polynomial((3, 2, 1)), oracle)
        self.assertEqual(f.integration_by_parts_residual((3, 2, 1)), 0)

    def test_signed_atomic_measure_coalesces_exactly(self):
        measure = DiracMeasure(((0, 3), (Fraction(0), -1), (2, -4), (2, 4)))
        self.assertEqual(measure.atoms, ((Fraction(0), Fraction(2)),))
        self.assertEqual(measure.total_variation(), 2)


if __name__ == "__main__":
    unittest.main()
