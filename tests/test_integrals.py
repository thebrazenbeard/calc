"""Independent exact polynomial-integration identities and failure bounds."""

from fractions import Fraction
from random import Random
import unittest

from calc.algebra import Polynomial
from calc.integrals import (
    definite_integral,
    polynomial_antiderivative,
    substitution_integral,
    integration_by_parts,
)


class ExactPolynomialIntegralTests(unittest.TestCase):
    def test_primitive_constant_and_fundamental_theorem(self):
        p = Polynomial((Fraction(2, 3), -3, 0, Fraction(5, 4)))
        primitive = polynomial_antiderivative(p, constant=Fraction(-7, 11))
        self.assertEqual(primitive(0), Fraction(-7, 11))
        self.assertEqual(primitive.differentiate(), p)
        self.assertEqual(primitive.coefficients[2], Fraction(-3, 2))

    def test_hand_computed_exact_definite_integrals(self):
        self.assertEqual(definite_integral(Polynomial((0, 1)), 0, 2), 2)
        self.assertEqual(definite_integral(Polynomial((0, 0, 1)), -1, 2), 3)
        self.assertEqual(definite_integral(Polynomial((1, 0, 3)), 0, 1), 2)
        self.assertEqual(definite_integral(Polynomial((0, 0, 0, 1)), -2, 2), 0)

    def test_orientation_and_zero_interval(self):
        p = Polynomial((Fraction(3, 7), -2, Fraction(5, 6)))
        a, b = Fraction(-4, 9), Fraction(5, 11)
        self.assertEqual(definite_integral(p, a, b), -definite_integral(p, b, a))
        self.assertEqual(definite_integral(p, a, a), Fraction(0))
        self.assertEqual(definite_integral(Polynomial((0,)), a, b), Fraction(0))

    def test_seeded_monomial_closed_form_oracles(self):
        rng = Random(0x51A7)
        for _ in range(80):
            degree = rng.randrange(0, 11)
            coefficient = Fraction(rng.randint(-15, 15), rng.randint(1, 9))
            a = Fraction(rng.randint(-4, 4), rng.randint(1, 6))
            b = Fraction(rng.randint(-4, 4), rng.randint(1, 6))
            polynomial = Polynomial((0,) * degree + (coefficient,))
            expected = coefficient * (b ** (degree + 1) - a ** (degree + 1)) / (degree + 1)
            self.assertEqual(definite_integral(polynomial, a, b), expected)

    def test_linearity(self):
        p = Polynomial((1, Fraction(2, 3), -1))
        q = Polynomial((Fraction(-1, 5), 0, 3, Fraction(1, 7)))
        a, b = -2, Fraction(3, 5)
        self.assertEqual(definite_integral(p + q, a, b),
                         definite_integral(p, a, b) + definite_integral(q, a, b))

    def test_substitution_closed_form(self):
        # ∫0¹ (x²+1)² (2x) dx = [((x²+1)³)/3]₀¹ = 7/3.
        self.assertEqual(
            substitution_integral(Polynomial((0, 0, 1)), Polynomial((1, 0, 1)), 0, 1),
            Fraction(7, 3),
        )
        self.assertEqual(
            substitution_integral(Polynomial((0, 0, 1)), Polynomial((1, 0, 1)), 1, 0),
            Fraction(-7, 3),
        )

    def test_substitution_affine_change_of_variables(self):
        # ∫0¹ (3x+2)^2 * 3 dx = [y³/3]₂⁵ = 39.
        self.assertEqual(
            substitution_integral(Polynomial((0, 0, 1)), Polynomial((2, 3)), 0, 1),
            39,
        )
        self.assertEqual(substitution_integral(Polynomial((1,)), Polynomial((7,)), 0, 1), 0)

    def test_integration_by_parts_exact(self):
        # ∫₀² (x²)' x³ dx = ∫₀² 2x⁴ dx = 64/5.
        u, v = Polynomial((0, 0, 1)), Polynomial((0, 0, 0, 1))
        self.assertEqual(integration_by_parts(u, v, 0, 2), Fraction(64, 5))
        self.assertEqual(integration_by_parts(u, v, 2, 0), Fraction(-64, 5))
        self.assertEqual(integration_by_parts(Polynomial((5,)), v, 0, 2), 0)

    def test_integration_by_parts_matches_direct_polynomial(self):
        rng = Random(201)
        for _ in range(60):
            u = Polynomial(rng.randint(-4, 4) for __ in range(rng.randrange(1, 6)))
            v = Polynomial(rng.randint(-4, 4) for __ in range(rng.randrange(1, 6)))
            a = Fraction(rng.randint(-3, 2), 3)
            b = Fraction(rng.randint(1, 4), 5)
            expected = definite_integral(u.differentiate() * v, a, b)
            self.assertEqual(integration_by_parts(u, v, a, b), expected)

    def test_all_results_are_exact_rationals(self):
        p = Polynomial((Fraction(1, 7), Fraction(2, 11)))
        self.assertIsInstance(definite_integral(p, Fraction(1, 8), Fraction(7, 9)), Fraction)
        self.assertIsInstance(substitution_integral(p, Polynomial((0, Fraction(1, 6))), 0, 1), Fraction)
        self.assertIsInstance(integration_by_parts(p, p, 0, 1), Fraction)

    def test_reject_inexact_inputs_and_nonpolynomials(self):
        p = Polynomial((1,))
        for value in (True, False, 1.0, "1", complex(1)):
            with self.subTest(value=repr(value)):
                with self.assertRaises(TypeError):
                    definite_integral(p, value, 1)
                with self.assertRaises(TypeError):
                    polynomial_antiderivative(p, constant=value)
        with self.assertRaises(TypeError):
            definite_integral((1, 2), 0, 1)
        with self.assertRaises(TypeError):
            substitution_integral(p, (1, 2), 0, 1)

    def test_resource_bound_degree(self):
        p = Polynomial((0,) * 65 + (1,))
        with self.assertRaises(ValueError):
            definite_integral(p, 0, 1)
        with self.assertRaises(ValueError):
            polynomial_antiderivative(p)

    def test_resource_bound_input_coefficient_and_endpoint(self):
        p = Polynomial((1 << 4100,))
        with self.assertRaises(ValueError):
            definite_integral(p, 0, 1)
        with self.assertRaises(ValueError):
            definite_integral(Polynomial((1,)), 0, 1 << 4100)

    def test_resource_bound_intermediate_growth(self):
        p = Polynomial((0,) * 64 + (1,))
        with self.assertRaises(ValueError):
            definite_integral(p, 0, 1 << 2048)

    def test_resource_bound_compositions(self):
        p = Polynomial((0,) * 33 + (1,))
        q = Polynomial((0, 0, 1))
        with self.assertRaises(ValueError):
            substitution_integral(p, q, 0, 1)
        with self.assertRaises(ValueError):
            integration_by_parts(Polynomial((0,) * 64 + (1,)), q, 0, 1)

    def test_no_mutation_of_polynomial_operands(self):
        p = Polynomial((3, -1, 2))
        q = Polynomial((2, 0, 4))
        before = (p.coefficients, q.coefficients)
        definite_integral(p, -2, 3)
        substitution_integral(p, q, -1, 1)
        integration_by_parts(p, q, -1, 1)
        self.assertEqual((p.coefficients, q.coefficients), before)


if __name__ == "__main__":
    unittest.main()
