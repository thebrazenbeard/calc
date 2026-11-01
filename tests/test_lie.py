"""Algebraic witnesses for the bounded exact matrix Lie algebra module."""

from fractions import Fraction
import unittest

from calc.lie import (
    adjoint_matrix,
    as_matrix,
    basis_coordinates,
    commutator,
    jacobiator,
    linear_combination,
    matrix_product,
    sl2_basis,
    so3_basis,
    structure_constants,
)


class MatrixLieTests(unittest.TestCase):
    def setUp(self):
        self.h, self.e, self.f = sl2_basis()

    def test_square_matrix_and_product_exact(self):
        a = as_matrix(((Fraction(1, 3), 2), (3, Fraction(1, 7))))
        self.assertEqual(a[0][0], Fraction(1, 3))
        self.assertEqual(matrix_product(a, ((2, 0), (0, 3))),
                         as_matrix(((Fraction(2, 3), 6), (6, Fraction(3, 7)))))

    def test_commutator_is_antisymmetric(self):
        a, b = self.h, self.e
        zero = as_matrix(((0, 0), (0, 0)))
        self.assertEqual(commutator(a, a), zero)
        self.assertEqual(commutator(a, b),
                         linear_combination((0, 2, 0), (self.h, self.e, self.f)))
        self.assertEqual(commutator(a, b),
                         linear_combination((-1,), (commutator(b, a),)))

    def test_sl2_commutator_relations(self):
        self.assertEqual(commutator(self.h, self.e),
                         linear_combination((2,), (self.e,)))
        self.assertEqual(commutator(self.h, self.f),
                         linear_combination((-2,), (self.f,)))
        self.assertEqual(commutator(self.e, self.f), self.h)

    def test_jacobi_identity(self):
        matrices = (
            as_matrix(((1, Fraction(1, 3)), (2, -3))),
            as_matrix(((2, 1), (Fraction(-1, 5), 4))),
            as_matrix(((Fraction(7, 4), 0), (3, Fraction(2, 7)))),
        )
        self.assertEqual(jacobiator(*matrices), as_matrix(((0, 0), (0, 0))))
        self.assertEqual(jacobiator(self.h, self.e, self.f),
                         as_matrix(((0, 0), (0, 0))))

    def test_jacobi_three_by_three(self):
        x, y, z = so3_basis()
        self.assertEqual(jacobiator(x, y, z), as_matrix(((0, 0, 0),) * 3))

    def test_rational_basis_coordinates(self):
        target = linear_combination((Fraction(2, 3), -7, Fraction(9, 5)),
                                    (self.h, self.e, self.f))
        self.assertEqual(basis_coordinates(target, (self.h, self.e, self.f)),
                         (Fraction(2, 3), Fraction(-7), Fraction(9, 5)))

    def test_sl2_structure_constants(self):
        c = structure_constants((self.h, self.e, self.f))
        self.assertEqual(len(c), 3)
        self.assertEqual(c[0][1], (0, 2, 0))
        self.assertEqual(c[0][2], (0, 0, -2))
        self.assertEqual(c[1][2], (1, 0, 0))
        for i in range(3):
            for j in range(3):
                self.assertEqual(c[i][j], tuple(-v for v in c[j][i]))

    def test_so3_bracket_and_structure(self):
        x, y, z = so3_basis()
        self.assertEqual(commutator(x, y), z)
        self.assertEqual(commutator(y, z), x)
        self.assertEqual(commutator(z, x), y)
        c = structure_constants((x, y, z))
        self.assertEqual(c[0][1], (0, 0, 1))
        self.assertEqual(c[1][2], (1, 0, 0))
        self.assertEqual(c[2][0], (0, 1, 0))

    def test_adjoint_matrix_column_convention(self):
        ad = adjoint_matrix(self.h, (self.h, self.e, self.f))
        self.assertEqual(ad, as_matrix(((0, 0, 0), (0, 2, 0), (0, 0, -2))))
        target = linear_combination((0, Fraction(5, 2), 3), (self.h, self.e, self.f))
        self.assertEqual(basis_coordinates(commutator(self.h, target),
                                           (self.h, self.e, self.f)), (0, 5, -6))

    def test_abelian_one_dimensional_algebra(self):
        identity = as_matrix(((1, 0), (0, 1)))
        self.assertEqual(structure_constants((identity,)), (((Fraction(0),),),))
        self.assertEqual(adjoint_matrix(identity, (identity,)), ((Fraction(0),),))

    def test_dependent_basis_rejected(self):
        with self.assertRaisesRegex(ValueError, "dependent"):
            basis_coordinates(self.h, (self.h,
                                       linear_combination((2,), (self.h,))))
        with self.assertRaisesRegex(ValueError, "dependent"):
            structure_constants((self.h, self.h))

    def test_matrix_outside_span_rejected(self):
        with self.assertRaisesRegex(ValueError, "not in the rational span"):
            basis_coordinates(self.e, (self.h,))

    def test_nonclosed_span_rejected(self):
        with self.assertRaisesRegex(ValueError, "not in the rational span"):
            structure_constants((self.e, self.f))

    def test_adjoint_requires_member(self):
        with self.assertRaisesRegex(ValueError, "not in the rational span"):
            adjoint_matrix(self.f, (self.h, self.e))

    def test_nonsquare_or_mismatched_dimensions_rejected(self):
        with self.assertRaisesRegex(ValueError, "square"):
            as_matrix(((1, 2), (3,)))
        with self.assertRaisesRegex(ValueError, "equal dimensions"):
            commutator(((1,),), ((1, 0), (0, 1)))
        with self.assertRaisesRegex(ValueError, "common dimension"):
            structure_constants((((1,),), self.h))

    def test_invalid_scalar_types_rejected(self):
        for value in (True, 0.5, complex(1, 0), "1", None):
            with self.subTest(value=value), self.assertRaises(TypeError):
                as_matrix(((value,),))
        with self.assertRaises(TypeError):
            linear_combination((0.5,), (self.h,))

    def test_resource_bounds_rejected(self):
        with self.assertRaisesRegex(ValueError, "between 1 and 8"):
            as_matrix(tuple(tuple(0 for _ in range(9)) for _ in range(9)))
        with self.assertRaisesRegex(ValueError, "bit-size"):
            as_matrix(((1 << 2049,),))
        with self.assertRaisesRegex(ValueError, "16 elements"):
            structure_constants(tuple(self.h for _ in range(17)))

    def test_invalid_coefficients(self):
        with self.assertRaisesRegex(ValueError, "count"):
            linear_combination((1, 2), (self.h,))


if __name__ == "__main__":
    unittest.main()
