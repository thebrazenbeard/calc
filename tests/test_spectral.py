"""Independent spectral witnesses for real symmetric matrices."""

import math
import unittest

from calc.spectral import (
    JacobiConvergenceError, jacobi_eigen, rayleigh_quotient,
    eigen_residuals, orthogonality_error,
)


class JacobiEigenTests(unittest.TestCase):
    def test_known_two_by_two_eigenpairs(self):
        matrix = ((2.0, 1.0), (1.0, 2.0))
        values, vectors = jacobi_eigen(matrix)
        self.assertAlmostEqual(values[0], 1.0, delta=1e-12)
        self.assertAlmostEqual(values[1], 3.0, delta=1e-12)
        self.assertLess(max(eigen_residuals(matrix, values, vectors)), 1e-12)
        self.assertLess(orthogonality_error(vectors), 1e-12)

    def test_three_by_three_repeated_eigenvalue(self):
        # 3I + J has eigenvalues 3, 3, 6 because rank(J)=1.
        matrix = ((4., 1., 1.), (1., 4., 1.), (1., 1., 4.))
        values, vectors = jacobi_eigen(matrix)
        for actual, truth in zip(values, (3., 3., 6.)):
            self.assertAlmostEqual(actual, truth, delta=1e-11)
        self.assertLess(max(eigen_residuals(matrix, values, vectors)), 1e-11)
        self.assertLess(orthogonality_error(vectors), 1e-11)

    def test_distinct_three_by_three_analytic_eigenvalues(self):
        # Tridiagonal Toeplitz eigenvalues: 2, 2 ± sqrt(2).
        matrix = ((2., 1., 0.), (1., 2., 1.), (0., 1., 2.))
        values, vectors = jacobi_eigen(matrix)
        for actual, truth in zip(values, (2 - math.sqrt(2), 2., 2 + math.sqrt(2))):
            self.assertAlmostEqual(actual, truth, delta=1e-12)
        self.assertLess(max(eigen_residuals(matrix, values, vectors)), 1e-12)

    def test_reconstruction_and_bounds(self):
        matrix = (
            (4., -1., .5, 0.),
            (-1., 3., 0., .25),
            (.5, 0., 2., .75),
            (0., .25, .75, 5.),
        )
        values, vectors = jacobi_eigen(matrix)
        self.assertLess(max(eigen_residuals(matrix, values, vectors)), 1e-11)
        self.assertLess(orthogonality_error(vectors), 1e-11)
        for i in range(4):
            for j in range(4):
                # The spectral decomposition A = V diag(lambda) V^T.
                reconstructed = sum(vectors[i][k] * values[k] * vectors[j][k]
                                    for k in range(4))
                self.assertAlmostEqual(reconstructed, matrix[i][j], delta=1e-10)
        # Independent Gershgorin disk bounds for every real eigenvalue.
        lo = min(matrix[i][i] - sum(abs(matrix[i][j]) for j in range(4) if i != j)
                 for i in range(4))
        hi = max(matrix[i][i] + sum(abs(matrix[i][j]) for j in range(4) if i != j)
                 for i in range(4))
        self.assertGreaterEqual(values[0], lo - 1e-12)
        self.assertLessEqual(values[-1], hi + 1e-12)
        self.assertAlmostEqual(sum(values), sum(matrix[i][i] for i in range(4)),
                               delta=1e-11)

        arbitrary_vector = (1., -2., 3., 4.)
        quotient = rayleigh_quotient(matrix, arbitrary_vector)
        self.assertGreaterEqual(quotient, values[0] - 1e-11)
        self.assertLessEqual(quotient, values[-1] + 1e-11)

    def test_zero_diagonal_singleton_and_sign(self):
        self.assertEqual(jacobi_eigen(((0., 0.), (0., 0.)))[0], (0., 0.))
        self.assertEqual(jacobi_eigen(((7.,),))[0], (7.,))
        self.assertEqual(jacobi_eigen(((-2., 0.), (0., 3.)))[0], (-2., 3.))

    def test_tiny_and_huge_scales(self):
        tiny = ((2e-300, 1e-300), (1e-300, 2e-300))
        values, vectors = jacobi_eigen(tiny)
        self.assertAlmostEqual(values[0] / 1e-300, 1.0, delta=1e-12)
        self.assertAlmostEqual(values[1] / 1e-300, 3.0, delta=1e-12)
        self.assertLess(max(eigen_residuals(tiny, values, vectors)), 1e-290)
        huge = ((1e308, 0.), (0., -1e308))
        self.assertEqual(jacobi_eigen(huge)[0], (-1e308, 1e308))

    def test_overflowed_true_eigenvalue_is_rejected(self):
        matrix = ((1e308, 1e308), (1e308, 1e308))
        with self.assertRaises(ArithmeticError):
            jacobi_eigen(matrix)

    def test_orthogonal_similarity_independent_spectrum(self):
        theta = 0.4
        c, s = math.cos(theta), math.sin(theta)
        # Rotate diagonal diag(1,4); eigenvalues invariant under Q D Q^T.
        matrix = ((c*c + 4*s*s, -3*c*s), (-3*c*s, s*s + 4*c*c))
        values, vectors = jacobi_eigen(matrix)
        self.assertAlmostEqual(values[0], 1., delta=1e-12)
        self.assertAlmostEqual(values[1], 4., delta=1e-12)
        self.assertLess(max(eigen_residuals(matrix, values, vectors)), 1e-12)

    def test_budget_and_invalid_matrices(self):
        matrix = ((2., 1., .3), (1., 3., .7), (.3, .7, 4.))
        with self.assertRaises(JacobiConvergenceError):
            jacobi_eigen(matrix, max_sweeps=1, tolerance=1e-30)
        for bad in (
            (), ((1., 2.),), ((1., 2.), (3., 4.)),
            ((1., float('nan')), (float('nan'), 1.)),
            ((True,),), ((1., 2., 3.), (2., 1., 4.), (3., 4.)),
        ):
            with self.subTest(bad=bad), self.assertRaises((ValueError, TypeError)):
                jacobi_eigen(bad)
        with self.assertRaises(ValueError):
            jacobi_eigen(matrix, tolerance=0)
        with self.assertRaises(ValueError):
            jacobi_eigen(matrix, max_sweeps=0)
        with self.assertRaises(TypeError):
            jacobi_eigen(matrix, max_sweeps=True)


class SpectralDiagnosticTests(unittest.TestCase):
    def test_rayleigh_on_exact_eigenvector(self):
        matrix = ((4., 1.), (1., 4.))
        self.assertAlmostEqual(rayleigh_quotient(matrix, (1., 1.)), 5.)
        self.assertAlmostEqual(rayleigh_quotient(matrix, (1., -1.)), 3.)

    def test_rayleigh_extreme_scaling(self):
        self.assertAlmostEqual(rayleigh_quotient(((2e-300,),), (1e308,)) / 1e-300, 2.)
        self.assertEqual(rayleigh_quotient(((0.,),), (1.,)), 0.)

    def test_diagnostics_detect_corrupt_vectors(self):
        matrix = ((2., 0.), (0., 3.))
        eigenvalues = (2., 3.)
        identity = ((1., 0.), (0., 1.))
        self.assertLess(max(eigen_residuals(matrix, eigenvalues, identity)), 1e-14)
        self.assertEqual(orthogonality_error(identity), 0.)
        wrong = ((1., 1.), (0., 1.))
        self.assertGreater(max(eigen_residuals(matrix, eigenvalues, wrong)), 0.1)
        self.assertGreater(orthogonality_error(wrong), 0.5)

    def test_diagnostic_errors(self):
        with self.assertRaises(ValueError):
            rayleigh_quotient(((2.,),), (0.,))
        with self.assertRaises(ValueError):
            rayleigh_quotient(((2.,),), (1., 2.))
        with self.assertRaises(ValueError):
            eigen_residuals(((2.,),), (2., 3.), ((1.,),))
        with self.assertRaises(ValueError):
            orthogonality_error(((1., 2.),))
        with self.assertRaises(ValueError):
            rayleigh_quotient(((1., 2.), (3., 4.)), (1., 2.))
        with self.assertRaises(ValueError):
            eigen_residuals(((2.,),), (float('inf'),), ((1.,),))
        with self.assertRaises(ValueError):
            orthogonality_error(((float('nan'),),))


if __name__ == "__main__":
    unittest.main()
