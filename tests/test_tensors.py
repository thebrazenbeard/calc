"""Exact independent component-index witnesses for bounded tensor arrays."""

from fractions import Fraction
from itertools import permutations
import unittest

from calc.tensors import Tensor, identity_tensor


class ExactTensorTests(unittest.TestCase):
    def test_rank_zero_scalar_and_empty_axis_tuple(self):
        scalar = Tensor.scalar(Fraction(7, 9))
        self.assertEqual(scalar.shape, ())
        self.assertEqual(scalar.rank, 0)
        self.assertEqual(scalar[()], Fraction(7, 9))
        self.assertEqual(scalar.tensor_product(Tensor((2,), (2, 3))),
                         Tensor((2,), (Fraction(14, 9), Fraction(7, 3))))

    def test_row_major_three_axis_indexing(self):
        tensor = Tensor((2, 3, 2), range(12))
        self.assertEqual(tensor[0, 0, 0], 0)
        self.assertEqual(tensor[0, 2, 1], 5)
        self.assertEqual(tensor[1, 0, 0], 6)
        self.assertEqual(tensor[1, 2, 1], 11)
        with self.assertRaises(IndexError):
            _ = tensor[1, 2]
        with self.assertRaises(IndexError):
            _ = tensor[2, 0, 0]
        with self.assertRaises(TypeError):
            _ = tensor[1, "0", 0]

    def test_finite_float_rejection_and_exact_mixed_rationals(self):
        tensor = Tensor((3,), (3, Fraction(-2, 7), Fraction(1, 4)))
        self.assertEqual(tensor.values, (Fraction(3), Fraction(-2, 7), Fraction(1, 4)))
        for invalid in (1.0, 0.5, True, complex(1, 0)):
            with self.subTest(invalid=invalid), self.assertRaises(TypeError):
                Tensor((1,), (invalid,))

    def test_exact_linear_operations(self):
        a = Tensor((2,), (Fraction(1, 2), 3))
        b = Tensor((2,), (Fraction(-2, 3), 1))
        self.assertEqual(a + b, Tensor((2,), (Fraction(-1, 6), 4)))
        self.assertEqual(a - b, a + (-b))
        self.assertEqual((a * 6) / 3, a * 2)
        self.assertEqual(a + Tensor.zeros((2,)), a)
        self.assertEqual((a * Fraction(2, 3))[0], Fraction(1, 3))
        with self.assertRaises(ZeroDivisionError):
            _ = a / 0
        with self.assertRaises(ValueError):
            _ = a + Tensor((1,), (3,))

    def test_outer_product_by_manual_components(self):
        u = Tensor((2,), (Fraction(1, 2), -3))
        v = Tensor((3,), (2, Fraction(-1, 4), 5))
        outer = u.tensor_product(v)
        self.assertEqual(outer.shape, (2, 3))
        for i in range(2):
            for j in range(3):
                self.assertEqual(outer[i, j], u[i] * v[j])
        self.assertEqual(outer.values, (1, Fraction(-1, 8), Fraction(5, 2),
                                        -6, Fraction(3, 4), -15))

    def test_outer_product_associativity_and_bilinearity(self):
        a = Tensor((2,), (1, Fraction(2, 3)))
        b = Tensor((2,), (Fraction(-1, 5), 4))
        c = Tensor((2,), (3, 7))
        self.assertEqual(a.tensor_product(b).tensor_product(c),
                         a.tensor_product(b.tensor_product(c)))
        self.assertEqual((a + b).tensor_product(c),
                         a.tensor_product(c) + b.tensor_product(c))
        self.assertEqual(c.tensor_product(a + b),
                         c.tensor_product(a) + c.tensor_product(b))

    def test_permutation_all_six_orderings_independent_coordinate_oracle(self):
        t = Tensor((2, 3, 2), range(12))
        for order in permutations((0, 1, 2)):
            with self.subTest(order=order):
                rearranged = t.permute_axes(order)
                self.assertEqual(rearranged.shape, tuple(t.shape[axis] for axis in order))
                for i in range(t.shape[0]):
                    for j in range(t.shape[1]):
                        for k in range(t.shape[2]):
                            coords = (i, j, k)
                            self.assertEqual(rearranged[tuple(coords[axis] for axis in order)],
                                             t[coords])

    def test_reverse_transpose_and_involution(self):
        t = Tensor((2, 3, 4), range(24))
        self.assertEqual(t.transpose().shape, (4, 3, 2))
        self.assertEqual(t.transpose().transpose(), t)
        self.assertEqual(Tensor.scalar(2).transpose(), Tensor.scalar(2))

    def test_arbitrary_axis_contraction_manual_oracle(self):
        t = Tensor((2, 3, 2, 4), range(48))
        traced = t.contract(0, 2)
        self.assertEqual(traced.shape, (3, 4))
        for j in range(3):
            for ell in range(4):
                self.assertEqual(traced[j, ell], t[0, j, 0, ell] + t[1, j, 1, ell])

    def test_trace_returns_rank_zero_and_exact_fraction(self):
        matrix = Tensor((2, 2), (Fraction(1, 3), 6, 7, Fraction(-1, 5)))
        result = matrix.contract(0, 1)
        self.assertEqual(result.shape, ())
        self.assertEqual(result[()], Fraction(2, 15))
        self.assertEqual(Tensor.zeros((4, 4)).contract(0, 1), Tensor.scalar(0))

    def test_matrix_product_equals_tensor_product_contraction(self):
        a = Tensor((2, 3), (1, 2, Fraction(-1, 3), 4, 5, 6))
        b = Tensor((3, 2), (2, 7, Fraction(1, 2), 0, -3, 1))
        product = a @ b
        self.assertEqual(product, a.tensor_product(b).contract(1, 2))
        self.assertEqual(product[0, 0], 1 * 2 + 2 * Fraction(1, 2) + Fraction(-1, 3) * -3)
        self.assertEqual(product[1, 1], 4 * 7 + 5 * 0 + 6 * 1)

    def test_kronecker_delta_contraction_and_matrix_identity(self):
        for n in range(1, 5):
            identity = identity_tensor(n)
            self.assertEqual(identity.contract(0, 1), Tensor.scalar(n))
            a = Tensor((n, n), [Fraction(i * n + j + 1, 3) for i in range(n) for j in range(n)])
            self.assertEqual(identity @ a, a)
            self.assertEqual(a @ identity, a)
            for i in range(n):
                for j in range(n):
                    self.assertEqual(sum(identity[i, k] * identity[k, j]
                                         for k in range(n)), int(i == j))

    def test_matrix_multiply_associativity(self):
        a = Tensor((2, 3), (1, Fraction(1, 2), -3, 0, 4, 2))
        b = Tensor((3, 2), (3, 1, Fraction(2, 3), 0, 5, 2))
        c = Tensor((2, 2), (1, 7, 2, -1))
        self.assertEqual((a @ b) @ c, a @ (b @ c))

    def test_symmetric_and_antisymmetric_decomposition(self):
        matrix = Tensor((3, 3), range(9))
        sym, anti = matrix.symmetrize(0, 1), matrix.antisymmetrize(0, 1)
        self.assertEqual(sym + anti, matrix)
        self.assertEqual(sym.transpose(), sym)
        self.assertEqual(anti.transpose(), -anti)
        for i in range(3):
            self.assertEqual(anti[i, i], 0)

    def test_symmetrize_interior_axes_three_tensor(self):
        t = Tensor((2, 3, 3), range(18))
        s = t.symmetrize(1, 2)
        a = t.antisymmetrize(1, 2)
        self.assertEqual(s + a, t)
        for i in range(2):
            for j in range(3):
                for k in range(3):
                    self.assertEqual(s[i, j, k], s[i, k, j])
                    self.assertEqual(a[i, j, k], -a[i, k, j])

    def test_reshape_keeps_row_major_values(self):
        a = Tensor((2, 3), range(6))
        self.assertEqual(a.reshape((3, 2)).values, a.values)
        self.assertEqual(a.reshape((6,)).reshape((2, 3)), a)
        with self.assertRaises(ValueError):
            a.reshape((5,))

    def test_invalid_axis_selection_and_shapes(self):
        t = Tensor((2, 3), range(6))
        for invalid in ((0, 0), (0, 2), (-1, 1), (True, 1)):
            with self.subTest(invalid=invalid):
                with self.assertRaises((TypeError, ValueError)):
                    t.contract(*invalid)
                with self.assertRaises((TypeError, ValueError)):
                    t.symmetrize(*invalid)
        with self.assertRaises(ValueError):
            t.contract(0, 1)
        with self.assertRaises(ValueError):
            t.permute_axes((1, 1))
        with self.assertRaises(ValueError):
            t.permute_axes((0,))
        with self.assertRaises(ValueError):
            Tensor((2, 2), range(3))
        with self.assertRaises(ValueError):
            Tensor((2, 2), range(5))
        with self.assertRaises(ValueError):
            Tensor((0,), [])
        with self.assertRaises(TypeError):
            Tensor((True,), [0])

    def test_rank_size_and_coefficient_resource_fences(self):
        with self.assertRaises(ValueError):
            Tensor((2,) * 7, [])
        with self.assertRaises(ValueError):
            Tensor((33,), range(33))
        with self.assertRaises(ValueError):
            Tensor((32, 32, 8), [])
        with self.assertRaises(ValueError):
            Tensor((1,), (1 << 4097,))
        with self.assertRaises(ValueError):
            Tensor((1,), (Fraction(1, 1 << 4097),))
        a = Tensor((32, 32), range(1024))
        with self.assertRaises(ValueError):
            a.tensor_product(a)

    def test_unbounded_generators_are_capped(self):
        def infinite():
            while True:
                yield 1
        with self.assertRaises(ValueError):
            Tensor(infinite(), ())
        with self.assertRaises(ValueError):
            Tensor((1,), infinite())

    def test_rejected_implicit_metric_or_broadcast(self):
        with self.assertRaises(ValueError):
            _ = Tensor((2, 1), (1, 2)) + Tensor((2,), (1, 2))
        with self.assertRaises(TypeError):
            _ = Tensor((1,), (1,)) * Tensor((1,), (2,))
        with self.assertRaises(TypeError):
            _ = Tensor((1,), (1,)) @ Tensor((1,), (2,))


if __name__ == "__main__":
    unittest.main()
