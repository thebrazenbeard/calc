"""Independent algebraic witnesses for bounded rational Clifford algebras."""

from fractions import Fraction
import itertools
import unittest

from calc.clifford import CliffordAlgebra, Multivector


class CliffordTests(unittest.TestCase):
    def setUp(self):
        self.alg = CliffordAlgebra((1, 1, -1))
        self.e = tuple(self.alg.basis(1 << i) for i in range(3))

    def test_units_signature(self):
        one = self.alg.scalar(1)
        for e, metric in zip(self.e, (1, 1, -1)):
            self.assertEqual(e * e, metric * one)

    def test_anticommutation(self):
        z = self.alg.zero()
        for i in range(3):
            for j in range(3):
                if i != j:
                    self.assertEqual(self.e[i] * self.e[j] + self.e[j] * self.e[i], z)

    def test_associativity_exhaustive_basis(self):
        blades = [self.alg.basis(mask) for mask in range(8)]
        for a, b, c in itertools.product(blades, repeat=3):
            self.assertEqual((a * b) * c, a * (b * c))

    def test_bivector_square_and_reversion(self):
        e01 = self.e[0].wedge(self.e[1])
        self.assertEqual(e01 * e01, self.alg.scalar(-1))
        self.assertEqual(e01.reverse(), -e01)
        self.assertEqual(self.e[0].reverse(), self.e[0])

    def test_wedge_graded_commutativity(self):
        a = self.e[0] + self.e[2]
        b = self.e[1].wedge(self.e[2])
        self.assertEqual(a.wedge(b), b.wedge(a))  # 1*2 even
        self.assertEqual(self.e[0].wedge(self.e[1]),
                         -self.e[1].wedge(self.e[0]))
        self.assertEqual(self.e[0].wedge(self.e[0]), self.alg.zero())

    def test_left_contraction(self):
        e0, e1, e2 = self.e
        self.assertEqual(e0.left_contract(e0.wedge(e1)), e1)
        self.assertEqual(e1.left_contract(e0.wedge(e1)), -e0)
        self.assertEqual(e2.left_contract(e2), self.alg.scalar(-1))
        self.assertEqual(e0.wedge(e1).left_contract(e0), self.alg.zero())

    def test_fraction_exactness_and_grade(self):
        a = self.alg.multivector({0: Fraction(1, 7), 1: Fraction(2, 3), 6: 1})
        self.assertEqual(a.scalar_part(), Fraction(1, 7))
        self.assertEqual(a.grade(1), Fraction(2, 3) * self.e[0])
        self.assertEqual(a.grade(2), self.alg.basis(6))
        self.assertEqual(a.grade(3), self.alg.zero())
        self.assertEqual(a.as_dict()[0], Fraction(1, 7))
        self.assertEqual((a * self.alg.scalar(Fraction(3, 2))).as_dict()[0], Fraction(3, 14))

    def test_involution_identities(self):
        a = self.alg.multivector({0: 2, 1: -3, 3: 4, 7: 5})
        b = self.alg.multivector({0: 3, 2: 1, 5: -1})
        self.assertEqual((a*b).reverse(), b.reverse() * a.reverse())
        self.assertEqual((a*b).grade_involution(),
                         a.grade_involution() * b.grade_involution())
        self.assertEqual(a.reverse().reverse(), a)
        self.assertEqual(a.clifford_conjugate(), a.reverse().grade_involution())

    def test_positive_euclidean_norm_and_indefinite_quadratic(self):
        euclidean = CliffordAlgebra((1, 1))
        a = euclidean.multivector({0: 2, 1: 3, 2: Fraction(1, 2), 3: 1})
        self.assertEqual(a.reversion_quadratic(), Fraction(57, 4))
        self.assertEqual(self.e[2].reversion_quadratic(), Fraction(-1))

    def test_zero_and_identity(self):
        a = self.alg.multivector({0: 4, 2: 7})
        self.assertEqual(a + (-a), self.alg.zero())
        self.assertEqual(a * self.alg.scalar(1), a)
        self.assertEqual(self.alg.scalar(1) * a, a)
        self.assertEqual(a * self.alg.zero(), self.alg.zero())
        self.assertEqual(3 * a, a * 3)

    def test_algebra_mismatch(self):
        other = CliffordAlgebra((1, -1))
        with self.assertRaises(ValueError):
            self.alg.scalar(1) + other.scalar(1)
        with self.assertRaises(ValueError):
            self.alg.scalar(1) * other.scalar(1)

    def test_input_validation(self):
        for metric in ((), (1, 0), (1, 2), (True, 1), (1.0, -1), (1,) * 9):
            with self.subTest(metric=metric), self.assertRaises((TypeError, ValueError)):
                CliffordAlgebra(metric)
        with self.assertRaises(TypeError):
            self.alg.multivector({1: 0.5})
        with self.assertRaises(TypeError):
            self.alg.scalar(True)
        with self.assertRaises(ValueError):
            self.alg.basis(8)
        with self.assertRaises(TypeError):
            self.alg.basis(True)
        with self.assertRaises(ValueError):
            self.alg.multivector({-1: 1})
        with self.assertRaises(ValueError):
            self.alg.scalar(Fraction(1, 1 << 5000))
        with self.assertRaises(ValueError):
            self.e[0].grade(-1)

    def test_immutable_storage_and_dict_copy(self):
        a = self.alg.multivector({2: 3, 1: 4, 7: 0})
        self.assertIsInstance(a, Multivector)
        self.assertEqual(tuple(a.as_dict()), (1, 2))
        copy = a.as_dict()
        copy[1] = Fraction(100)
        self.assertEqual(a.as_dict()[1], Fraction(4))

    def test_wedge_distributivity_and_mixed_grades(self):
        a = self.alg.scalar(2) + self.e[0]
        b = self.e[1] + self.e[2]
        c = self.alg.scalar(Fraction(1, 5)) + self.e[2]
        self.assertEqual(a.wedge(b + c), a.wedge(b) + a.wedge(c))
        self.assertEqual(a.wedge(b).wedge(c), a.wedge(b.wedge(c)))

    def test_left_contraction_grade_mixed(self):
        a = self.alg.scalar(3) + self.e[0]
        b = self.alg.scalar(2) + self.e[0].wedge(self.e[1])
        self.assertEqual(a.left_contract(b),
                         3 * b + self.e[1])

    def test_dimension_bound_and_generator_validation(self):
        big = CliffordAlgebra((1,) * 8)
        self.assertEqual(big.basis(255).grade(8), big.basis(255))
        self.assertEqual(big.positive, 8)
        self.assertEqual(big.negative, 0)
        with self.assertRaises(ValueError):
            big.generator(8)
        with self.assertRaises(TypeError):
            big.generator(True)

    def test_arithmetic_result_resource_bound(self):
        max_accepted = self.alg.scalar(1 << 4095)
        with self.assertRaises(ValueError):
            _ = max_accepted * 2


if __name__ == "__main__":
    unittest.main()
