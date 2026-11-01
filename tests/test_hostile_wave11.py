"""Independent bounded adversarial oracles for Wave 11 mathematics.

These public-API checks derive expected values from finite enumeration,
elementary coefficient identities, and independently written group actions.
They are not a formal verification of unrestricted algebraic structures.
"""

from fractions import Fraction
from itertools import product
import random
import unittest

from calc.affine_geometry import AffineSystem, partial_derivative
from calc.asymptotics import RationalAsymptotic, compare_growth
from calc.grobner import MultiPolynomial
from calc.modules import ResidueFreeModule
from calc.representations import (
    RationalRepresentation, character_pairing, cyclic_group, direct_sum, regular_representation,
    trivial_representation,
)
from calc.rings import ResidueRing, chinese_remainder


class HostileFiniteRingTests(unittest.TestCase):
    def test_exhaustive_principal_ideals_and_annihilators(self):
        for modulus in range(2, 24):
            ring = ResidueRing(modulus)
            for generator in range(modulus):
                ideal = ring.ideal(generator)
                expected = {(multiplier * generator) % modulus
                            for multiplier in range(modulus)}
                with self.subTest(modulus=modulus, generator=generator):
                    self.assertEqual(set(ideal.elements()), expected)
                    self.assertEqual(ideal.order, len(expected))
                    self.assertEqual(
                        set(ring.annihilator(generator).elements()),
                        {x for x in range(modulus) if (x * generator) % modulus == 0},
                    )

    def test_ideal_sum_intersection_product_direct_set_oracle(self):
        for modulus in range(2, 15):
            ring = ResidueRing(modulus)
            for a, b in product(range(modulus), repeat=2):
                left, right = ring.ideal(a), ring.ideal(b)
                aa, bb = set(left.elements()), set(right.elements())
                self.assertEqual(set(left.sum(right).elements()),
                                 {(u + v) % modulus for u in aa for v in bb})
                self.assertEqual(set(left.intersection(right).elements()), aa & bb)
                self.assertEqual(set(left.product(right).elements()),
                                 {(u * v) % modulus for u in aa for v in bb})

    def test_crt_finite_independent_congruences(self):
        moduli = (4, 9, 5, 7)
        for value in (0, 1, 17, 319, 1259):
            solution, modulus = chinese_remainder((value % m for m in moduli), moduli)
            self.assertEqual(modulus, 1260)
            self.assertEqual(solution, value % modulus)
        with self.assertRaises(ValueError):
            chinese_remainder((1, 2), (4, 6))


class HostileFiniteRepresentationTests(unittest.TestCase):
    def test_regular_action_character_and_matrix_oracle(self):
        for order in range(1, 10):
            group = cyclic_group(order)
            representation = regular_representation(group)
            character = representation.character()
            self.assertEqual(character[group.identity], order)
            for exponent, element in enumerate(group.elements):
                matrix = representation.matrix(element)
                for row, column in product(range(order), repeat=2):
                    self.assertEqual(matrix[row][column], int(row == (column + exponent) % order))
                self.assertEqual(character[element], order if exponent == 0 else 0)
            self.assertEqual(character_pairing(representation, trivial_representation(group)), 1)

    def test_direct_sum_character_and_intertwiners(self):
        from calc.representations import intertwines
        for order in range(2, 7):
            group = cyclic_group(order)
            trivial = trivial_representation(group)
            repeated = direct_sum(trivial, trivial)
            self.assertEqual(repeated.dimension, 2)
            self.assertEqual(set(repeated.character().values()), {Fraction(2)})
            self.assertTrue(intertwines(trivial, repeated, ((1,), (3,))))
            self.assertTrue(intertwines(repeated, trivial, ((1, 5),)))

    def test_nontrivial_character_orthogonality_and_invalid_generator(self):
        c2 = cyclic_group(2)
        sign = RationalRepresentation(c2, {"e": ((1,),), "g": ((-1,),)})
        trivial = trivial_representation(c2)
        self.assertEqual(character_pairing(sign, trivial), 0)
        self.assertEqual(character_pairing(sign, sign), 1)
        c3 = cyclic_group(3)
        with self.assertRaisesRegex(ValueError, "homomorphism"):
            RationalRepresentation(c3, {
                "e": ((1,),), "g": ((-1,),), "g2": ((1,),)
            })


class HostileAffineGeometryTests(unittest.TestCase):
    def test_circle_grid_and_exact_jacobian(self):
        circle = MultiPolynomial({(2, 0): 1, (0, 2): 1, (0, 0): -1}, nvars=2)
        system = AffineSystem((circle,))
        grid = (-1, 0, 1)
        self.assertEqual(set(system.points_on_grid((grid, grid))),
                         {(Fraction(1), Fraction(0)), (Fraction(-1), Fraction(0)),
                          (Fraction(0), Fraction(1)), (Fraction(0), Fraction(-1))})
        self.assertEqual(system.jacobian_at((0, 1)), ((Fraction(0), Fraction(2)),))
        self.assertEqual(partial_derivative(circle, 0),
                         MultiPolynomial({(1, 0): 2}, nvars=2))

    def test_ideal_membership_differs_from_point_sample(self):
        x2_minus_1 = MultiPolynomial({(2,): 1, (0,): -1}, nvars=1)
        system = AffineSystem((x2_minus_1,))
        x4_minus_1 = MultiPolynomial({(4,): 1, (0,): -1}, nvars=1)
        x_minus_1 = MultiPolynomial({(1,): 1, (0,): -1}, nvars=1)
        self.assertTrue(system.ideal_contains(x4_minus_1))
        self.assertFalse(system.ideal_contains(x_minus_1))
        self.assertTrue(system.contains((1,)))
        self.assertTrue(system.contains((-1,)))
        irrational_roots = AffineSystem((MultiPolynomial({(2,): 1, (0,): -2}, nvars=1),))
        self.assertEqual(irrational_roots.points_on_grid(((-2, -1, 0, 1, 2),)), ())


class HostileAsymptoticsTests(unittest.TestCase):
    def test_random_exact_degree_classification(self):
        rng = random.Random(20261101)
        for _ in range(110):
            first_degree, second_degree = rng.randrange(6), rng.randrange(6)
            a, b = rng.choice((-7, -3, 1, 2, 5)), rng.choice((-5, -2, 1, 3, 7))
            first = RationalAsymptotic((0,) * first_degree + (a,))
            second = RationalAsymptotic((0,) * second_degree + (b,))
            infinity = compare_growth(first, second)
            near_zero = compare_growth(first, second, at="zero_plus")
            difference = first_degree - second_degree
            self.assertEqual(infinity.big_o, difference <= 0)
            self.assertEqual(infinity.little_o, difference < 0)
            self.assertEqual(infinity.theta, difference == 0)
            self.assertEqual(near_zero.big_o, difference >= 0)
            self.assertEqual(near_zero.little_o, difference > 0)
            self.assertEqual(near_zero.theta, difference == 0)
            if difference == 0:
                self.assertEqual(infinity.ratio_limit, Fraction(a, b))
                self.assertEqual(near_zero.ratio_limit, Fraction(a, b))


class HostileFiniteModuleTests(unittest.TestCase):
    def test_generated_submodule_against_independent_scalar_enumeration(self):
        rng = random.Random(20261102)
        for modulus in (2, 3, 4, 5):
            for rank in (1, 2, 3):
                ambient = ResidueFreeModule(modulus, rank)
                for _ in range(8):
                    generators = [tuple(rng.randrange(modulus) for _ in range(rank))
                                  for _ in range(rng.randrange(4))]
                    actual = ambient.generated_submodule(*generators)
                    independent = set()
                    # Every integer coefficient tuple suffices modulo n.
                    for coefficients in product(range(modulus), repeat=len(generators)):
                        vector = tuple(sum(k * gen[i] for k, gen in zip(coefficients, generators)) % modulus
                                       for i in range(rank))
                        independent.add(vector)
                    self.assertEqual(set(actual.elements), independent)
                    self.assertEqual(actual.index * actual.order, ambient.order)

    def test_quotient_cosets_and_representative_independence(self):
        for modulus, rank, generators in (
            (4, 2, ((2, 0),)), (3, 2, ((1, 1),)),
            (5, 1, ((0,),)), (2, 3, ((1, 1, 0), (0, 1, 1))),
        ):
            ambient = ResidueFreeModule(modulus, rank)
            submodule = ambient.generated_submodule(*generators)
            quotient = submodule.quotient()
            self.assertEqual(quotient.order * submodule.order, ambient.order)
            flat = tuple(value for coset in quotient.cosets for value in coset)
            self.assertEqual(len(flat), len(set(flat)))
            self.assertEqual(set(flat), set(ambient.elements()))
            for index, coset in enumerate(quotient.cosets):
                for representative in coset:
                    self.assertEqual(quotient.class_of(representative), index)
                    for scalar in range(modulus):
                        self.assertEqual(quotient.scalar_multiply(scalar, index),
                                         quotient.class_of(ambient.scalar_multiply(scalar, representative)))
                self.assertEqual(quotient.add(index, quotient.zero_class), index)
                self.assertEqual(quotient.add(index, quotient.neg(index)), quotient.zero_class)

    def test_invalid_large_module_and_nonintegral_coordinates(self):
        with self.assertRaises(ValueError):
            ResidueFreeModule(2, 13)
        module = ResidueFreeModule(4, 2)
        with self.assertRaises(TypeError):
            module.normalize((True, 0))
        with self.assertRaises(TypeError):
            module.scalar_multiply(0.5, (0, 1))
        with self.assertRaises(ValueError):
            module.normalize((1,))


if __name__ == "__main__":
    unittest.main()
