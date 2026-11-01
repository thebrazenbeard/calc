"""Independent adversarial checks for the Wave 8 mathematical primitives.

The oracles below are analytic identities or independently specified finite
structures. They do not derive expected results from production algorithms.
"""

from fractions import Fraction
from itertools import product, repeat
import math
import random
import unittest

from calc.categories import (
    FiniteFunctor, NaturalTransformation, poset_category, monoid_category,
)
from calc.collocation import (
    barycentric_interpolate, chebyshev_differentiation_matrix,
    chebyshev_lobatto_nodes, differentiate_at_nodes,
)
from calc.homology import ChainComplex, rational_rank
from calc.tensors import Tensor, identity_tensor
from calc.umbral import (
    falling_factorial, rising_factorial, stirling_first, stirling_second,
    monomial_to_falling, falling_to_monomial, forward_difference,
    newton_forward_interpolate,
)


def _eval(coefficients, x):
    return sum((coefficient * x**degree for degree, coefficient in
                enumerate(coefficients)), Fraction(0))


class HostileSpectralCollocation(unittest.TestCase):
    def test_polynomial_identity_interpolation_and_nodal_derivative(self):
        # A degree-seven polynomial is reproduced by degree-eight Lobatto data.
        n = 8
        nodes = chebyshev_lobatto_nodes(n)
        p = lambda x: 1.25 - 3*x + 2*x**3 - x**5 + .125*x**7
        dp = lambda x: -3 + 6*x*x - 5*x**4 + .875*x**6
        nodal = tuple(p(x) for x in nodes)
        for sample in (-.913, -.531, -.07, .0, .123, .78, .987):
            self.assertAlmostEqual(barycentric_interpolate(nodal, sample), p(sample), delta=1e-12)
        derivatives = differentiate_at_nodes(nodal)
        for actual, x in zip(derivatives, nodes):
            self.assertAlmostEqual(actual, dp(x), delta=2e-10)

    def test_matrix_differentiates_linear_and_annihilates_constants(self):
        for n in (1, 2, 5, 17, 33):
            nodes = chebyshev_lobatto_nodes(n)
            matrix = chebyshev_differentiation_matrix(n)
            for row in matrix:
                self.assertLess(abs(math.fsum(row)), 1e-10)
                self.assertAlmostEqual(math.fsum(a * x for a, x in zip(row, nodes)),
                                       1., delta=2e-10)

    def test_nodes_reversed_and_exact_at_saved_nodes(self):
        nodes = chebyshev_lobatto_nodes(16)
        samples = tuple(math.cos(3*x) for x in nodes)
        self.assertEqual(nodes[0], 1.)
        self.assertEqual(nodes[-1], -1.)
        self.assertTrue(all(a > b for a, b in zip(nodes, nodes[1:])))
        for node, sample in zip(nodes, samples):
            self.assertEqual(barycentric_interpolate(samples, node), sample)


class HostileUmbralCalculus(unittest.TestCase):
    def test_stirling_inverse_triangular_matrices(self):
        # The signed first-kind and second-kind Stirling matrices are inverses.
        for n in range(13):
            for k in range(n + 1):
                actual = sum(stirling_second(n, j) * stirling_first(j, k)
                             for j in range(k, n + 1))
                self.assertEqual(actual, int(n == k), (n, k))

    def test_seeded_rational_polynomial_basis_roundtrip(self):
        random_source = random.Random(202611011)
        for degree in range(12):
            for _ in range(6):
                p = tuple(Fraction(random_source.randint(-17, 17),
                                   random_source.randint(1, 17))
                          for _ in range(degree + 1))
                restored = falling_to_monomial(monomial_to_falling(p))
                normalized = tuple(p[:degree + 1])
                while len(normalized) > 1 and not normalized[-1]:
                    normalized = normalized[:-1]
                self.assertEqual(restored, normalized)

    def test_falling_factorial_forward_step_and_negative_grid(self):
        for n in range(1, 16):
            for x in (Fraction(-8, 3), Fraction(1, 2), Fraction(9, 5)):
                self.assertEqual(falling_factorial(x + 1, n) - falling_factorial(x, n),
                                 n * falling_factorial(x, n - 1))
                self.assertEqual(rising_factorial(x, n),
                                 (-1)**n * falling_factorial(-x, n))
        polynomial = (Fraction(3, 7), -2, Fraction(4, 5), 0, -3)
        start, step = Fraction(7, 4), Fraction(-3, 8)
        samples = tuple(_eval(polynomial, start + index * step) for index in range(5))
        self.assertEqual(newton_forward_interpolate(samples, start=start, step=step), polynomial)
        delta = forward_difference(polynomial, step=step, order=3)
        for x in (Fraction(-7, 3), Fraction(1, 4)):
            # Independent binomial shift formula Δ_h^3 p(x).
            oracle = sum(((-1)**(3-j) * math.comb(3, j) *
                          _eval(polynomial, x + j * step) for j in range(4)), Fraction(0))
            self.assertEqual(_eval(delta, x), oracle)


class HostileFiniteCategories(unittest.TestCase):
    def test_diamond_poset_two_paths_have_unique_composite(self):
        objs = ("bottom", "left", "right", "top")
        relation = {(a, a) for a in objs} | {
            ("bottom", "left"), ("bottom", "right"),
            ("left", "top"), ("right", "top"), ("bottom", "top")}
        category = poset_category(objs, relation)
        bl = category.hom("bottom", "left")[0]
        br = category.hom("bottom", "right")[0]
        lt = category.hom("left", "top")[0]
        rt = category.hom("right", "top")[0]
        self.assertEqual(category.compose(lt, bl), category.compose(rt, br))
        self.assertEqual(category.hom("left", "right"), ())

    def test_natural_transformation_of_monotone_endofunctors(self):
        objs = ("a", "b", "c")
        relation = {(a, b) for i, a in enumerate(objs) for b in objs[i:]}
        category = poset_category(objs, relation)

        def functor(map_obj):
            arrow_map = {arrow.name: category.hom(
                map_obj[arrow.source], map_obj[arrow.target])[0]
                for arrow in category.arrows}
            return FiniteFunctor(category, category, map_obj, arrow_map)

        f = functor({"a": "a", "b": "a", "c": "b"})
        g = functor({"a": "b", "b": "b", "c": "c"})
        components = {obj: category.hom(f.on_objects[obj], g.on_objects[obj])[0]
                      for obj in objs}
        eta = NaturalTransformation(f, g, components)
        for arrow in category.arrows:
            left = category.compose(g.on_arrows[arrow.name], eta.components[arrow.source])
            right = category.compose(eta.components[arrow.target], f.on_arrows[arrow.name])
            self.assertEqual(left, right)

    def test_bad_monoid_associativity_rejected(self):
        # A unital but nonassociative multiplication table.
        elements = ("e", "a", "b")
        table = {(g, f): (f if g == "e" else g if f == "e" else
                          {("a", "a"): "a", ("a", "b"): "a",
                           ("b", "a"): "b", ("b", "b"): "a"}[(g, f)])
                 for g in elements for f in elements}
        with self.assertRaisesRegex(ValueError, "associativity"):
            monoid_category(elements, "e", table)

    def test_group_of_order_three_and_nonnatural_automorphisms(self):
        elements = ("e", "a", "b")
        multiplication = {(g, f): elements[(i + j) % 3]
                          for i, g in enumerate(elements)
                          for j, f in enumerate(elements)}
        category = monoid_category(elements, "e", multiplication)
        original = FiniteFunctor(category, category, {"•": "•"},
                                 {name: name for name in elements})
        inversion = FiniteFunctor(category, category, {"•": "•"},
                                  {"e": "e", "a": "b", "b": "a"})
        for eta in elements:
            with self.subTest(component=eta), self.assertRaisesRegex(ValueError, "naturality"):
                NaturalTransformation(original, inversion, {"•": eta})


class HostileRationalHomology(unittest.TestCase):
    def test_oriented_square_cycle_and_filled_disc(self):
        # Four oriented edges around a cycle; rank ∂1 = |V|-1 = 3.
        edge_boundary = ((-1, 0, 0, 1),
                         (1, -1, 0, 0),
                         (0, 1, -1, 0),
                         (0, 0, 1, -1))
        perimeter = ChainComplex((4, 4), {1: edge_boundary})
        self.assertEqual(perimeter.betti_numbers(), (1, 1))
        disc = ChainComplex((4, 4, 1), {1: edge_boundary,
                                        2: ((1,), (1,), (1,), (1,))})
        self.assertEqual(disc.betti_numbers(), (1, 0, 0))
        self.assertEqual(disc.euler_characteristic(), 1)
        self.assertEqual(disc.homology_euler_characteristic(), 1)

    def test_rational_basis_scaling_preserves_betti(self):
        rng = random.Random(202611012)
        for _ in range(60):
            a, b, c = (Fraction(rng.choice((-1, 1)) * rng.randint(1, 12),
                                rng.randint(1, 12)) for _ in range(3))
            # ∂1=(a,b,c), ∂2 spans vectors (b,-a,0),(c,0,-a), so ∂²=0.
            chain = ChainComplex((1, 3, 2), {
                1: ((a, b, c),),
                2: ((b, c), (-a, 0), (0, -a)),
            })
            self.assertEqual(chain.betti_numbers(), (0, 0, 0))
            self.assertEqual(chain.euler_characteristic(), 0)

    def test_exact_rank_oracle_fraction_cancellation(self):
        # Third row equals first row minus half the second, identically.
        first = (Fraction(1, 7), Fraction(2, 5), Fraction(1, 2))
        second = (Fraction(3, 2), Fraction(-5, 3), Fraction(4, 9))
        third = tuple(a - b / 2 for a, b in zip(first, second))
        self.assertEqual(rational_rank((first, second, third)), 2)
        with self.assertRaisesRegex(ValueError, "chain identity"):
            ChainComplex((1, 1, 1), {1: ((1,),), 2: ((1,),)})


class HostileTensorAlgebra(unittest.TestCase):
    def test_kronecker_delta_contraction_equals_matrix(self):
        matrix = Tensor((2, 2), (Fraction(1, 3), -5, Fraction(7, 9), 2))
        identity = identity_tensor(2)
        self.assertEqual(matrix.tensor_product(identity).contract(1, 2), matrix)
        self.assertEqual(identity.tensor_product(matrix).contract(1, 2), matrix)

    def test_axis_permutation_inverse_over_nonsymmetric_rectangular_shape(self):
        values = tuple(Fraction(3*i-11, i+2) for i in range(24))
        tensor = Tensor((2, 3, 4), values)
        permuted = tensor.permute_axes((1, 2, 0))
        self.assertEqual(permuted.shape, (3, 4, 2))
        self.assertEqual(permuted.permute_axes((2, 0, 1)), tensor)
        for i, j, k in product(range(2), range(3), range(4)):
            self.assertEqual(tensor[i, j, k], permuted[j, k, i])

    def test_symmetry_decomposition_reconstructs_original(self):
        tensor = Tensor((3, 3), tuple(Fraction(i*i-5*i+2, i+1) for i in range(9)))
        symmetric = tensor.symmetrize(0, 1)
        antisymmetric = tensor.antisymmetrize(0, 1)
        self.assertEqual(symmetric + antisymmetric, tensor)
        self.assertEqual(symmetric.transpose(), symmetric)
        self.assertEqual(antisymmetric.transpose(), -antisymmetric)
        self.assertEqual(antisymmetric.contract(0, 1), Tensor.scalar(0))

    def test_exact_matrix_product_associativity(self):
        a = Tensor((2, 2), (Fraction(1, 3), 2, -3, Fraction(5, 4)))
        b = Tensor((2, 2), (1, Fraction(-2, 7), 4, 0))
        c = Tensor((2, 2), (2, -1, Fraction(5, 2), 7))
        self.assertEqual((a @ b) @ c, a @ (b @ c))
        self.assertEqual((a @ identity_tensor(2)), a)


class HostileResourceBoundaries(unittest.TestCase):
    def test_infinite_inputs_are_bounded(self):
        # Every lazy input should be consumed no farther than a documented cap.
        with self.assertRaises(ValueError):
            Tensor((2, 2), repeat(1))
        with self.assertRaises(ValueError):
            newton_forward_interpolate(repeat(1))
        with self.assertRaises(ValueError):
            barycentric_interpolate(repeat(1.), 0.)
        with self.assertRaises(ValueError):
            poset_category(repeat("x"), ())

    def test_rational_and_dimensional_limits_enforced(self):
        with self.assertRaises(ValueError):
            ChainComplex((33,))
        with self.assertRaises(ValueError):
            Tensor((33,), (1,) * 33)
        with self.assertRaises(ValueError):
            falling_factorial(1, 33)
        with self.assertRaises(ValueError):
            chebyshev_lobatto_nodes(129)


if __name__ == "__main__":
    unittest.main()
