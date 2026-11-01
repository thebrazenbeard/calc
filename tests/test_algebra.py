"""Independent examples and edge cases for exact algebra routines."""

from fractions import Fraction
from itertools import permutations
from math import prod
from random import Random
import unittest

from calc.algebra import (
    Polynomial,
    chinese_remainder,
    determinant,
    extended_gcd,
    inverse_matrix,
    matmul,
    matrix_rank,
    mod_inverse,
    solve_linear_system,
    transpose,
    trace,
)


class PolynomialTests(unittest.TestCase):
    def test_normalization_and_degree(self):
        self.assertEqual(Polynomial([1, 0, 0]).coefficients, (Fraction(1),))
        self.assertEqual(Polynomial([]).coefficients, (Fraction(0),))
        self.assertEqual(Polynomial([0]).degree, -1)

    def test_exact_arithmetic(self):
        p = Polynomial([1, 2, 3])
        q = Polynomial([3, -1])
        self.assertEqual((p + q).coefficients, (4, 1, 3))
        self.assertEqual((p - q).coefficients, (-2, 3, 3))
        self.assertEqual((p * q).coefficients, (3, 5, 7, -3))
        self.assertEqual((-q).coefficients, (-3, 1))
        self.assertEqual(p(2), 17)

    def test_exact_derivative_and_antiderivative(self):
        p = Polynomial([7, 2, 3])
        self.assertEqual(p.differentiate().coefficients, (2, 6))
        self.assertEqual(p.differentiate(2).coefficients, (6,))
        self.assertEqual(p.differentiate(3).coefficients, (0,))
        self.assertEqual(p.integrate(5).coefficients, (5, 7, 1, 1))
        self.assertEqual(p.integrate().differentiate(), p)

    def test_rational_coefficients_remain_exact(self):
        p = Polynomial([Fraction(1, 3), Fraction(1, 2)])
        self.assertEqual(p(Fraction(2, 3)), Fraction(2, 3))
        self.assertEqual(p.integrate().coefficients, (0, Fraction(1, 3), Fraction(1, 4)))

    def test_invalid_polynomial_inputs(self):
        with self.assertRaises(TypeError):
            Polynomial([0.1])
        with self.assertRaises(ValueError):
            Polynomial([1]).differentiate(-1)
        with self.assertRaises(TypeError):
            Polynomial([1]).integrate(0.5)

    def test_polynomial_division_and_gcd(self):
        numerator = Polynomial([-1, 0, 1])  # x**2 - 1
        divisor = Polynomial([-1, 1])  # x - 1
        self.assertEqual(divmod(numerator, divisor), (Polynomial([1, 1]), Polynomial()))
        self.assertEqual(divmod(Polynomial([1, 0, 1]), Polynomial([1, 1])),
                         (Polynomial([-1, 1]), Polynomial([2])))
        self.assertEqual(numerator.gcd(divisor), divisor)
        self.assertEqual(Polynomial([0, 0]).gcd(Polynomial([-2, 2])), divisor)
        self.assertEqual(Polynomial().gcd(Polynomial()), Polynomial())
        with self.assertRaises(ZeroDivisionError):
            divmod(numerator, Polynomial())

    def test_nonnegative_integer_power(self):
        self.assertEqual((Polynomial([1, 1]) ** 3).coefficients, (1, 3, 3, 1))
        self.assertEqual(Polynomial() ** 0, Polynomial([1]))
        with self.assertRaises(ValueError):
            Polynomial([1]) ** -1


class MatrixTests(unittest.TestCase):
    def test_transpose_and_multiply(self):
        self.assertEqual(transpose([[1, 2, 3], [4, 5, 6]]),
                         ((1, 4), (2, 5), (3, 6)))
        self.assertEqual(matmul([[1, 2], [3, 4]], [[2, 0], [1, 2]]),
                         ((4, 4), (10, 8)))

    def test_determinant_exact_including_swap(self):
        self.assertEqual(determinant([[2, 1, 3], [1, -1, 0], [4, 2, 1]]), 15)
        self.assertEqual(determinant([[0, 1], [5, 2]]), -5)
        self.assertEqual(determinant([[1, 2], [2, 4]]), 0)
        self.assertEqual(determinant([[Fraction(1, 2)]]), Fraction(1, 2))

    def test_linear_system_unique_fractional_solution(self):
        self.assertEqual(solve_linear_system([[2, 1], [1, -1]], [5, 1]), (2, 1))
        self.assertEqual(solve_linear_system([[2, 1], [1, 2]], [1, 0]),
                         (Fraction(2, 3), Fraction(-1, 3)))
        self.assertEqual(solve_linear_system([[0, 1], [1, 0]], [4, 7]), (7, 4))

    def test_invalid_shape_and_singularity(self):
        with self.assertRaises(ValueError):
            transpose([[1], [2, 3]])
        with self.assertRaises(ValueError):
            matmul([[1, 2]], [[1, 2]])
        with self.assertRaises(ValueError):
            determinant([[1, 2, 3]])
        with self.assertRaises(ValueError):
            solve_linear_system([[1, 2], [2, 4]], [3, 6])
        with self.assertRaises(ValueError):
            solve_linear_system([[1, 2], [2, 4]], [3, 7])
        with self.assertRaises(ValueError):
            solve_linear_system([[1, 2], [3, 4]], [1])
        with self.assertRaises(ValueError):
            determinant([])

    def test_exact_inputs_not_floats(self):
        with self.assertRaises(TypeError):
            determinant([[0.1]])

    def test_small_determinants_against_permutation_formula(self):
        # Independent Leibniz-formula oracle, not another elimination routine.
        random = Random(3901)
        for size in (1, 2, 3, 4):
            for _ in range(16):
                matrix = [[random.randrange(-4, 5) for _ in range(size)] for _ in range(size)]
                expected = 0
                for perm in permutations(range(size)):
                    inversions = sum(perm[i] > perm[j]
                                     for i in range(size) for j in range(i + 1, size))
                    expected += (-1 if inversions % 2 else 1) * prod(
                        matrix[i][perm[i]] for i in range(size)
                    )
                self.assertEqual(determinant(matrix), expected, matrix)

    def test_solutions_against_substitution(self):
        random = Random(222)
        for _ in range(40):
            # Strict diagonal dominance guarantees a nonsingular matrix.
            matrix = [[random.randrange(-2, 3) for _ in range(3)] for _ in range(3)]
            for i in range(3):
                matrix[i][i] = 10
            solution = tuple(random.randrange(-4, 5) for _ in range(3))
            rhs = [sum(row[j] * solution[j] for j in range(3)) for row in matrix]
            self.assertEqual(solve_linear_system(matrix, rhs), solution)

    def test_rank_inverse_and_trace(self):
        self.assertEqual(matrix_rank([[1, 2, 3], [2, 4, 6], [1, 0, 1]]), 2)
        self.assertEqual(matrix_rank([[0, 0], [0, 0]]), 0)
        self.assertEqual(matrix_rank([[1, 0, 0], [0, 1, 0]]), 2)
        matrix = [[1, 2], [3, 5]]
        self.assertEqual(inverse_matrix(matrix), ((-5, 2), (3, -1)))
        self.assertEqual(matmul(matrix, inverse_matrix(matrix)), ((1, 0), (0, 1)))
        self.assertEqual(trace(matrix), 6)
        with self.assertRaises(ValueError):
            inverse_matrix([[1, 2], [2, 4]])
        with self.assertRaises(ValueError):
            trace([[1, 2, 3]])


class ModularArithmeticTests(unittest.TestCase):
    def test_extended_gcd_identity(self):
        for a, b, expected_gcd in [(-30, 12, 6), (0, 0, 0), (17, 5, 1)]:
            g, x, y = extended_gcd(a, b)
            self.assertEqual(g, expected_gcd)
            self.assertEqual(a*x + b*y, g)

    def test_inverse(self):
        self.assertEqual(mod_inverse(3, 11), 4)
        self.assertEqual(mod_inverse(-3, 11), 7)
        with self.assertRaises(ValueError):
            mod_inverse(6, 9)
        with self.assertRaises(ValueError):
            mod_inverse(2, 1)

    def test_generalized_chinese_remainder(self):
        self.assertEqual(chinese_remainder([(2, 3), (3, 5)]), (8, 15))
        self.assertEqual(chinese_remainder([(2, 4), (6, 8)]), (6, 8))
        self.assertEqual(chinese_remainder([(-1, 7), (0, 2)]), (6, 14))
        self.assertEqual(chinese_remainder([]), (0, 1))
        with self.assertRaises(ValueError):
            chinese_remainder([(1, 2), (0, 4)])
        with self.assertRaises(ValueError):
            chinese_remainder([(1, 0)])

    def test_small_modular_inverse_exhaustive(self):
        for modulus in range(2, 24):
            for residue in range(-30, 31):
                candidates = [i for i in range(modulus) if residue * i % modulus == 1]
                if candidates:
                    self.assertEqual(mod_inverse(residue, modulus), candidates[0])
                else:
                    with self.assertRaises(ValueError):
                        mod_inverse(residue, modulus)


if __name__ == "__main__":
    unittest.main()
