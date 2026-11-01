"""Independent exact-rational finite operator-algebra witnesses."""

from fractions import Fraction
import unittest

from calc.operator_algebra import (
    adjoint, anticommutator, commutator, identity, is_normal,
    is_orthogonal_projection, is_positive_semidefinite_2x2,
    is_self_adjoint, product, quadratic_form, star_algebra, trace,
)


E11 = ((1, 0), (0, 0))
E12 = ((0, 1), (0, 0))
E21 = ((0, 0), (1, 0))
E22 = ((0, 0), (0, 1))


class OperatorAlgebraTests(unittest.TestCase):
    def test_identity_and_adjoint(self):
        a = ((Fraction(1, 3), 2), (-5, Fraction(7, 4)))
        self.assertEqual(identity(2), ((1, 0), (0, 1)))
        self.assertEqual(adjoint(adjoint(a)), a)
        self.assertEqual(adjoint(a), ((Fraction(1, 3), -5), (2, Fraction(7, 4))))

    def test_adjoint_reverses_products(self):
        a, b = ((2, 3), (7, -1)), ((0, 9), (4, 8))
        self.assertEqual(adjoint(product(a, b)), product(adjoint(b), adjoint(a)))

    def test_commutator_antisymmetry_and_trace_zero(self):
        a, b = ((1, 2), (3, 4)), ((0, 1), (5, 6))
        ab, ba = commutator(a, b), commutator(b, a)
        self.assertEqual(ab, tuple(tuple(-x for x in row) for row in ba))
        self.assertEqual(trace(ab), 0)

    def test_anticommutator_of_symmetric_is_symmetric(self):
        a, b = ((1, 2), (2, 3)), ((3, 5), (5, 9))
        self.assertEqual(adjoint(anticommutator(a, b)), anticommutator(a, b))
        self.assertEqual(anticommutator(E11, E11), ((2, 0), (0, 0)))

    def test_trace_cyclicity_independent(self):
        a, b = ((Fraction(1, 2), 3), (-2, 7)), ((8, -1), (9, 6))
        self.assertEqual(trace(product(a, b)), trace(product(b, a)))
        self.assertEqual(trace(product(a, b)), Fraction(75))

    def test_quadratic_form_exact(self):
        matrix = ((2, Fraction(1, 2)), (Fraction(1, 2), 3))
        x = (Fraction(1, 3), Fraction(-2, 5))
        expected = Fraction(2, 9) - Fraction(2, 15) + Fraction(12, 25)
        self.assertEqual(quadratic_form(matrix, x), expected)

    def test_orthogonal_projection(self):
        p = ((Fraction(1, 2), Fraction(1, 2)),
             (Fraction(1, 2), Fraction(1, 2)))
        self.assertTrue(is_orthogonal_projection(p))
        self.assertTrue(is_self_adjoint(p))
        self.assertEqual(product(p, p), p)
        self.assertFalse(is_orthogonal_projection(((1, 1), (0, 0))))

    def test_normal_and_nonnormal(self):
        self.assertTrue(is_normal(((0, -1), (1, 0))))
        self.assertTrue(is_normal(((2, 0), (0, 4))))
        self.assertFalse(is_normal(E12))

    def test_psd_2x2_boundary(self):
        self.assertTrue(is_positive_semidefinite_2x2(((1, 2), (2, 4))))
        self.assertTrue(is_positive_semidefinite_2x2(((0, 0), (0, 0))))
        self.assertTrue(is_positive_semidefinite_2x2(((2, Fraction(1, 2)),
                                                      (Fraction(1, 2), 2))))
        self.assertFalse(is_positive_semidefinite_2x2(((1, 2), (2, 1))))
        self.assertFalse(is_positive_semidefinite_2x2(((-1, 0), (0, 9))))

    def test_psd_2x2_rejects_bad_domain(self):
        with self.assertRaises(ValueError):
            is_positive_semidefinite_2x2(((1, 2), (3, 4)))
        with self.assertRaises(ValueError):
            is_positive_semidefinite_2x2(((1,),))

    def test_diagonal_star_algebra(self):
        structure = star_algebra((E11, E22))
        self.assertEqual(structure.unit, (1, 1))
        self.assertEqual(structure.involution, ((1, 0), (0, 1)))
        self.assertEqual(structure.multiplication[0][0], (1, 0))
        self.assertEqual(structure.multiplication[0][1], (0, 0))
        self.assertEqual(structure.multiplication[1][1], (0, 1))

    def test_full_matrix_algebra_structure(self):
        basis = (E11, E12, E21, E22)
        structure = star_algebra(basis)
        self.assertEqual(structure.unit, (1, 0, 0, 1))
        self.assertEqual(structure.multiplication[1][2], (1, 0, 0, 0))
        self.assertEqual(structure.multiplication[2][1], (0, 0, 0, 1))
        self.assertEqual(structure.involution[1], (0, 0, 1, 0))

    def test_nonunital_star_algebra_allowed_when_requested(self):
        structure = star_algebra((E11,), require_unit=False)
        self.assertIsNone(structure.unit)
        self.assertEqual(structure.multiplication, (((Fraction(1),),),))
        with self.assertRaises(ValueError):
            star_algebra((E11,))

    def test_rejects_nonclosed_product_span(self):
        a = ((0, 0, 0), (0, 1, 0), (0, 0, 2))
        basis = (identity(3), a)
        with self.assertRaisesRegex(ValueError, "span"):
            star_algebra(basis)

    def test_rejects_nonclosed_adjoint_span(self):
        basis = (identity(2), E12)
        with self.assertRaisesRegex(ValueError, "span"):
            star_algebra(basis)

    def test_rejects_dependent_basis(self):
        with self.assertRaises(ValueError):
            star_algebra((E11, E11))

    def test_dimension_mismatch_and_shapes(self):
        with self.assertRaises(ValueError):
            product(E11, identity(3))
        with self.assertRaises(ValueError):
            commutator(E11, identity(3))
        with self.assertRaises(ValueError):
            quadratic_form(E11, (1,))
        with self.assertRaises(ValueError):
            trace(((1, 2, 3), (4, 5, 6)))

    def test_rejects_inexact_and_bool(self):
        with self.assertRaises(TypeError):
            adjoint(((1.0,),))
        with self.assertRaises(TypeError):
            trace(((True,),))
        with self.assertRaises(TypeError):
            quadratic_form(((1,),), (False,))
        with self.assertRaises(TypeError):
            star_algebra((E11,), require_unit=1)

    def test_dimension_and_scalar_bounds(self):
        with self.assertRaises(ValueError):
            identity(9)
        with self.assertRaises(ValueError):
            identity(True)
        with self.assertRaises(ValueError):
            trace(((1 << 4096,),))
        with self.assertRaises(ValueError):
            quadratic_form(((1,),), (1 << 4096,))


if __name__ == "__main__":
    unittest.main()
