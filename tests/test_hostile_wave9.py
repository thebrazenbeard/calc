"""Independent Wave 9 mathematical oracles and bounded adversarial inputs.

The expected answers below come from finite measure identities, free words,
matrix multiplication, and Bellman completion of squares rather than from
calling one implementation helper to predict another.
"""

from fractions import Fraction as F
import itertools
import random
import unittest

from calc.measure import FiniteMeasure
from calc.free_algebra import FreeAlgebra
from calc.operator_algebra import (
    adjoint, product, trace, is_positive_semidefinite_2x2,
    star_algebra,
)
from calc.control import solve_scalar_lqr, evaluate_trajectory, optimal_trajectory
from calc.differential_algebra import JetPolynomial


class HostileFiniteMeasureTests(unittest.TestCase):
    def test_radon_nikodym_exact_change_of_measure_with_null_atoms(self):
        rng = random.Random(99135)
        for _ in range(60):
            a, b, c = [F(rng.randrange(0, 8), rng.randrange(1, 7)) for _ in range(3)]
            ref = FiniteMeasure({"a": a, "b": b, "c": c})
            densities = {atom: F(rng.randrange(0, 9), rng.randrange(1, 5))
                         for atom in ref.atoms}
            target = FiniteMeasure({atom: mass * densities[atom]
                                    for atom, mass in zip(ref.atoms, ref.weights)})
            derivative = ref.density_of(target)
            self.assertEqual(ref.integrate(derivative), target.total_mass)
            for atom, mass in zip(ref.atoms, ref.weights):
                self.assertEqual(derivative[atom] * mass,
                                 dict(zip(target.atoms, target.weights))[atom])
                if mass == 0:
                    self.assertEqual(derivative[atom], 0)

    def test_conditional_expectation_preserves_cell_integrals_and_tower(self):
        rng = random.Random(91022)
        for _ in range(75):
            weights = {str(i): F(rng.randrange(0, 9), rng.randrange(1, 7))
                       for i in range(5)}
            values = {str(i): F(rng.randrange(-13, 14), rng.randrange(1, 7))
                      for i in range(5)}
            mu = FiniteMeasure(weights)
            partition = (("0", "1"), ("2", "3"), ("4",))
            conditioned = mu.conditional_expectation(values, partition)
            self.assertEqual(mu.integrate(values), mu.integrate(conditioned))
            for cell in partition:
                selector = {atom: (1 if atom in cell else 0) for atom in mu.atoms}
                self.assertEqual(
                    mu.integrate({a: selector[a] * values[a] for a in mu.atoms}),
                    mu.integrate({a: selector[a] * conditioned[a] for a in mu.atoms}),
                )
            whole = mu.conditional_expectation(conditioned, (mu.atoms,))
            self.assertEqual(mu.integrate(whole), mu.integrate(values))

    def test_singularity_and_zero_mass_conventions(self):
        base = FiniteMeasure({"a": 0, "b": F(1, 3)})
        with self.assertRaises(ValueError):
            base.density_of(FiniteMeasure({"a": 1, "b": 0}))
        self.assertEqual(base.conditional_expectation({"a": 7, "b": -3},
                                                     (("a",), ("b",))),
                         {"a": F(0), "b": F(-3)})
        with self.assertRaises(ValueError):
            FiniteMeasure({"a": 0}).normalize()
        with self.assertRaises(ValueError):
            base.conditional_expectation({"a": 1, "b": 2}, (("a",), ("a", "b")))


class HostileFreeAlgebraTests(unittest.TestCase):
    def setUp(self):
        self.a = FreeAlgebra(("x", "y"))
        self.x, self.y = self.a.variable("x"), self.a.variable("y")

    def test_noncommuting_words_and_exact_custom_derivation(self):
        x, y = self.x, self.y
        self.assertNotEqual(x * y, y * x)
        self.assertEqual((x * y * x).derivative("x"), y * x + x * y)
        images = {"x": x * y, "y": x}
        p = x * y * x
        expected = x * y * y * x + x * x * x + x * y * x * y
        self.assertEqual(p.derivation(images), expected)

    def test_seeded_free_leibniz_and_jacobi(self):
        rng = random.Random(12209)
        a = self.a
        words = [(), ("x",), ("y",), ("x", "y"), ("y", "x"),
                 ("x", "x"), ("y", "y", "x")]
        for _ in range(55):
            polys = [a.zero() for _ in range(3)]
            for i in range(3):
                for word in rng.sample(words, 3):
                    polys[i] += rng.randrange(-3, 4) * a.monomial(word)
            p, q, r = polys
            for var in a.generators:
                self.assertEqual((p * q).derivative(var),
                                 p.derivative(var) * q + p * q.derivative(var))
            self.assertEqual(
                p.commutator(q.commutator(r))
                + q.commutator(r.commutator(p))
                + r.commutator(p.commutator(q)), a.zero())

    def test_noncommutative_input_resource_guards(self):
        with self.assertRaises(ValueError):
            FreeAlgebra(("x", "x"))
        with self.assertRaises((ValueError, TypeError)):
            self.a.monomial(("z",))
        with self.assertRaises(ValueError):
            _ = self.x ** -1


class HostileFiniteOperatorTests(unittest.TestCase):
    def test_adjoint_product_reversal_trace_cyclicity_seeded(self):
        rng = random.Random(8889)
        for _ in range(45):
            a = tuple(tuple(F(rng.randrange(-5, 6), rng.randrange(1, 5))
                            for _ in range(3)) for _ in range(3))
            b = tuple(tuple(F(rng.randrange(-5, 6), rng.randrange(1, 5))
                            for _ in range(3)) for _ in range(3))
            self.assertEqual(adjoint(product(a, b)), product(adjoint(b), adjoint(a)))
            self.assertEqual(trace(product(a, b)), trace(product(b, a)))

    def test_psd_2x2_as_quadratic_form(self):
        # PSD iff all real quadratic forms are nonnegative. The forward
        # implication is checked on an independently enumerated rational grid.
        rng = random.Random(5007)
        from calc.operator_algebra import quadratic_form
        for _ in range(40):
            v = [F(rng.randrange(-6, 7), rng.randrange(1, 5)) for _ in range(2)]
            a = ((v[0] * v[0], v[0] * v[1]),
                 (v[0] * v[1], v[1] * v[1]))
            self.assertTrue(is_positive_semidefinite_2x2(a))
            for x, y in itertools.product(range(-3, 4), repeat=2):
                self.assertGreaterEqual(quadratic_form(a, (x, y)), 0)
        self.assertFalse(is_positive_semidefinite_2x2(((0, 1), (1, 0))))

    def test_matrix_star_algebra_closed_and_bad_span(self):
        d1 = ((1, 0), (0, 0))
        d2 = ((0, 0), (0, 1))
        good = star_algebra((d1, d2))
        self.assertEqual(good.unit, (F(1), F(1)))
        self.assertEqual(good.multiplication[0][1], (F(0), F(0)))
        with self.assertRaises(ValueError):
            star_algebra((d1,))  # unital closure requires I outside span


class HostileFiniteControlTests(unittest.TestCase):
    def test_completed_square_one_step_independent_closed_form(self):
        rng = random.Random(2305)
        for _ in range(65):
            a = F(rng.randrange(-5, 6), rng.randrange(1, 6))
            b = F(rng.randrange(-5, 6), rng.randrange(1, 6))
            q = F(rng.randrange(0, 7), rng.randrange(1, 6))
            r = F(rng.randrange(1, 8), rng.randrange(1, 6))
            qt = F(rng.randrange(0, 7), rng.randrange(1, 6))
            x = F(rng.randrange(-6, 7), rng.randrange(1, 6))
            sol = solve_scalar_lqr(a, b, q, r, qt, 1)
            denom = r + qt * b * b
            independent_u = -a * b * qt * x / denom
            independent_v = q*x*x + qt*a*a*x*x - (qt*a*b*x)**2 / denom
            self.assertEqual(sol.control(x, 0), independent_u)
            self.assertEqual(sol.value(x), independent_v)
            self.assertEqual(optimal_trajectory(sol, x).total_cost, independent_v)
            for candidate in (-2, -1, 0, 1, 2, independent_u + F(1, 17)):
                candidate_cost = evaluate_trajectory(sol, x, (candidate,)).total_cost
                self.assertGreaterEqual(candidate_cost, independent_v)

    def test_two_stage_bellman_identity_and_zero_horizon(self):
        for n in range(4):
            sol = solve_scalar_lqr(F(3, 2), F(2, 3), F(1, 5), F(2, 7), F(5, 4), n)
            self.assertEqual(len(sol.cost_to_go), n+1)
            for x in (F(-2, 3), F(0), F(7, 5)):
                path = optimal_trajectory(sol, x)
                self.assertEqual(path.total_cost, sol.value(x))
                for i in range(n):
                    state, chosen = path.states[i], path.controls[i]
                    self.assertEqual(
                        sol.value(state, i), sol.q*state*state + sol.r*chosen*chosen
                        + sol.value(path.states[i+1], i+1))


class HostileDifferentialAlgebraTests(unittest.TestCase):
    def test_exact_jet_leibniz_and_formal_chain_rule(self):
        u = JetPolynomial.jet(0, max_jet=3)
        up = JetPolynomial.jet(1, max_jet=3)
        upp = JetPolynomial.jet(2, max_jet=3)
        self.assertEqual((u**2 * up**3).total_derivative(),
                         2*u*up**4 + 3*u**2*up**2*upp)
        self.assertEqual((u**5).total_derivative(), 5*u**4*up)
        self.assertEqual(JetPolynomial.constant(7, max_jet=3).total_derivative(),
                         JetPolynomial.constant(0, max_jet=3))

    def test_seeded_jet_derivation_against_independent_monomial_oracle(self):
        rng = random.Random(66819)
        for _ in range(65):
            coefficients: dict[tuple[int, ...], F] = {}
            for _ in range(4):
                powers = (rng.randrange(0, 4), rng.randrange(0, 4), 0, 0)
                coefficients[powers] = F(rng.randrange(-6, 7), rng.randrange(1, 8))
            p = JetPolynomial(3, coefficients)
            values = tuple(F(rng.randrange(-7, 8), rng.randrange(1, 6)) for _ in range(4))
            # Dp = sum_i u_(i+1) ∂_i p, independently evaluated as a
            # rational polynomial in four freely chosen jet values.
            expected = F(0)
            for exponents, coefficient in p.coefficients:
                for i in range(3):
                    if exponents[i]:
                        expected += (coefficient * exponents[i] * values[i+1]
                                     * _independent_monomial(values, exponents, i))
            self.assertEqual(p.total_derivative().evaluate(values), expected)

    def test_seeded_jet_leibniz_and_truncation_rejection(self):
        rng = random.Random(2029)
        for _ in range(55):
            a = JetPolynomial(3, {(rng.randrange(0, 3), rng.randrange(0, 3), 0, 0):
                                  F(rng.randrange(1, 4), rng.randrange(1, 4))})
            b = JetPolynomial(3, {(rng.randrange(0, 3), rng.randrange(0, 3), 0, 0):
                                  F(rng.randrange(1, 4), rng.randrange(1, 4))})
            self.assertEqual((a*b).total_derivative(),
                             a.total_derivative()*b + a*b.total_derivative())
        top = JetPolynomial.jet(3, max_jet=3)
        with self.assertRaises(ValueError):
            top.total_derivative()
        with self.assertRaises(TypeError):
            JetPolynomial(3, {(1,0,0,0): 1.2})


def _independent_monomial(values: tuple[F, ...], powers: tuple[int, ...],
                          differentiated_index: int) -> F:
    product = F(1)
    for i, power in enumerate(powers):
        product *= values[i] ** (power - int(i == differentiated_index))
    return product


if __name__ == "__main__":
    unittest.main()
