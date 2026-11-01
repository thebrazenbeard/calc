"""Independent contract properties spanning algebra and analysis.

These oracles use mathematical invariants or closed forms, not the library's
own intermediate routines. They are intentionally deterministic.
"""

from fractions import Fraction
import math
import unittest

from calc.advanced import wedge
from calc.algebra import Polynomial, chinese_remainder, determinant, matmul
from calc.differential import curl, gradient
from calc.numerical import rk4, simpson


class AlgebraContractTests(unittest.TestCase):
    def test_polynomial_distributive_law_exact(self):
        a = Polynomial([Fraction(1, 2), -3, 7])
        b = Polynomial([4, Fraction(2, 3), 1])
        c = Polynomial([-9, Fraction(1, 7)])
        self.assertEqual(a * (b + c), a * b + a * c)
        self.assertEqual((a * b).differentiate(),
                         a.differentiate() * b + a * b.differentiate())

    def test_determinant_multiplicativity_exact(self):
        a = ((1, 2, 3), (-1, 0, 2), (2, 3, -1))
        b = ((2, 1, 0), (3, 2, -2), (0, 1, 1))
        self.assertEqual(determinant(matmul(a, b)), determinant(a) * determinant(b))

    def test_chinese_remainder_matches_brute_force_small_moduli(self):
        for m in range(1, 9):
            for n in range(1, 9):
                for a in range(m):
                    for b in range(n):
                        solutions = [x for x in range(math.lcm(m, n))
                                     if x % m == a and x % n == b]
                        if not solutions:
                            with self.assertRaises(ValueError):
                                chinese_remainder(((a, m), (b, n)))
                        else:
                            residue, modulus = chinese_remainder(((a, m), (b, n)))
                            self.assertEqual(modulus, math.lcm(m, n))
                            self.assertEqual(residue, solutions[0])


class CalculusContractTests(unittest.TestCase):
    def test_vector_calculus_gradient_has_zero_curl(self):
        # For a smooth scalar potential phi, curl(grad(phi)) = 0.
        phi = lambda p: p[0] * p[1] + p[2] ** 2 + p[1] ** 3
        field = tuple(lambda p, i=i: gradient(phi, p)[i] for i in range(3))
        for component in curl(field, (1.0, -0.5, 2.0)):
            self.assertAlmostEqual(component, 0.0, delta=3e-4)

    def test_simpson_exact_for_independent_cubic(self):
        f = lambda x: -2 * x**3 + 5 * x**2 - 9 * x + 4
        primitive = lambda x: -x**4 / 2 + 5*x**3 / 3 - 9*x*x / 2 + 4*x
        self.assertAlmostEqual(simpson(f, -3, 2, 4), primitive(2) - primitive(-3), places=11)

    def test_rk4_error_reduces_with_finer_mesh(self):
        # y'=y, y(0)=1 -> y(1)=e. RK4 global error O(h^4).
        error_coarse = abs(rk4(lambda t, y: y, 0, 1, 1, 5)[-1][1] - math.e)
        error_fine = abs(rk4(lambda t, y: y, 0, 1, 1, 10)[-1][1] - math.e)
        self.assertGreater(error_coarse, error_fine * 10)

    def test_exterior_graded_anticommutativity_and_associativity(self):
        alpha = {(0,): Fraction(2, 3), (2,): 1}
        beta = {(1,): 3, (3,): -1}
        gamma = {(4,): Fraction(-5, 2)}
        ab = wedge(alpha, beta)
        ba = wedge(beta, alpha)
        self.assertEqual({indices: -coeff for indices, coeff in ab.items()}, ba)
        self.assertEqual(wedge(ab, gamma), wedge(alpha, wedge(beta, gamma)))


if __name__ == "__main__":
    unittest.main()
