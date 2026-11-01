"""Independent analytic and negative witnesses for matrix calculus."""

from fractions import Fraction
import math
import unittest

from calc.matrix_calc import (
    frobenius_inner, matrix_gradient, matrix_jvp,
    squared_frobenius_gradient, trace_axb_gradient,
)


class MatrixCalculusTests(unittest.TestCase):
    def test_exact_frobenius_inner(self):
        a = ((Fraction(1, 2), 2), (3, 4))
        b = ((2, 3), (4, 5))
        self.assertEqual(frobenius_inner(a, b), Fraction(39))
        self.assertIsInstance(frobenius_inner(a, b), Fraction)

    def test_float_frobenius_inner(self):
        self.assertAlmostEqual(frobenius_inner(((0.1, 0.2),), ((1.0, 2.0),)), 0.5)

    def test_frobenius_dimensions_and_type(self):
        for a, b in [(((1, 2),), ((1,),)), (((1, 2), (3,)), ((1, 2), (3, 4))),
                     (((True,),), ((1,),)), (((math.nan,),), ((2,),))]:
            with self.subTest(a=a, b=b), self.assertRaises((ValueError, TypeError)):
                frobenius_inner(a, b)

    def test_squared_gradient_exact(self):
        self.assertEqual(squared_frobenius_gradient(((1, 2), (-3, 4))),
                         ((Fraction(2), Fraction(4)), (Fraction(-6), Fraction(8))))

    def test_trace_axb_gradient_known(self):
        # F(X) = tr(A X B): gradient = (B A)^T, independent analytic formula.
        self.assertEqual(trace_axb_gradient(((1, 2), (3, 4)), ((5, 6), (7, 8))),
                         ((Fraction(23), Fraction(31)), (Fraction(34), Fraction(46))))

    def test_trace_axb_gradient_rectangular(self):
        a = ((1, 2, 3), (4, 5, 6))  # 2x3
        b = ((1, 2), (0, 1), (2, 0), (3, -1))  # 4x2
        result = trace_axb_gradient(a, b)
        self.assertEqual((len(result), len(result[0])), (3, 4))
        for i in range(3):
            for j in range(4):
                self.assertEqual(result[i][j], sum(a[k][i] * b[j][k] for k in range(2)))

    def test_trace_axb_incompatible_shapes(self):
        with self.assertRaises(ValueError):
            trace_axb_gradient(((1, 2), (3, 4)), ((5, 6, 7),))

    def test_frobenius_rejects_overflowed_float_product(self):
        with self.assertRaises(ValueError):
            frobenius_inner(((1e308,),), ((1e308,),))

    def test_trace_formula_vs_independent_perturbation(self):
        a = ((1, 2, 3), (4, 5, 6))
        b = ((1, 2), (0, 1), (2, 0), (3, -1))
        x = ((0.5, -1.0, 0.0, 2.0), (3.0, 2.0, -2.0, 1.0), (1.0, 0.0, 4.0, -1.0))
        def scalar_trace(m):
            # Full three-matrix contraction, independently differentiated by
            # matrix_gradient rather than reusing analytic gradient entries.
            return sum(sum(sum(a[k][i]*m[i][j] for i in range(3))*b[j][k]
                           for j in range(4)) for k in range(2))
        numeric = matrix_gradient(scalar_trace, x)
        analytic = trace_axb_gradient(a, b)
        for i in range(3):
            for j in range(4):
                self.assertAlmostEqual(numeric[i][j], analytic[i][j], places=5)

    def test_gradient_independent_quadratic_oracle(self):
        x = ((1.0, 2.0), (-3.0, 4.0))
        def quadratic(m):
            return sum(value * value for row in m for value in row)
        result = matrix_gradient(quadratic, x)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(result[i][j], 2*x[i][j], places=5)

    def test_gradient_step_errors(self):
        for step in (0, -1, math.inf, math.nan):
            with self.subTest(step=step), self.assertRaises(ValueError):
                matrix_gradient(lambda x: x[0][0], ((1,),), step=step)
        with self.assertRaises(ValueError):
            matrix_gradient(lambda x: x[0][0], ((1e20,),), step=1e-20)
        with self.assertRaises(ValueError):
            matrix_gradient(lambda x: x[0][0], ((float.fromhex('0x1.fffffffffffffp+1023'),),))

    def test_jvp_linear_exact_oracle(self):
        x = ((2.0, 1.0), (0.0, -3.0))
        tangent = ((1.0, -2.0), (3.0, 4.0))
        def transform(m):
            return tuple(tuple(2*v+1 for v in row) for row in m)
        result = matrix_jvp(transform, x, tangent)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(result[i][j], 2*tangent[i][j], places=5)

    def test_jvp_shape_and_output_errors(self):
        with self.assertRaises(ValueError):
            matrix_jvp(lambda x: x, ((1.0, 2.0),), ((1.0,),))
        with self.assertRaises(ValueError):
            matrix_jvp(lambda x: ((math.inf,),), ((1.0,),), ((1.0,),))


if __name__ == "__main__":
    unittest.main()
