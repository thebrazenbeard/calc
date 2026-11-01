"""Independent algebraic witnesses for small exact Gröbner computations."""

from fractions import Fraction
import unittest

from calc.grobner import (
    MultiPolynomial, ResourceLimitError, grobner_basis, leading_term,
    normal_form, s_polynomial, variables,
)


class MultiPolynomialTests(unittest.TestCase):
    def setUp(self):
        self.x, self.y = variables(2)

    def test_exact_product_and_evaluation(self):
        p = (self.x + Fraction(1, 2) * self.y) * (self.x - self.y)
        self.assertEqual(p.evaluate((2, 4)), Fraction(-8))
        self.assertEqual(p, self.x**2 - Fraction(1, 2)*self.x*self.y - Fraction(1, 2)*self.y**2)

    def test_distributivity_and_nontrivial_fraction(self):
        f, g, h = self.x+Fraction(1, 2)*self.y, self.x-3, self.y**2+2
        self.assertEqual(f*(g+h), f*g+f*h)
        self.assertEqual((f*g).evaluate((Fraction(2, 3), -1)),
                         f.evaluate((Fraction(2, 3), -1))*g.evaluate((Fraction(2, 3), -1)))

    def test_canonical_zero_and_equality(self):
        self.assertEqual(self.x-self.x, MultiPolynomial.zero(2))
        self.assertEqual((self.x-self.x).terms, ())
        self.assertEqual(self.x+self.y, self.y+self.x)
        self.assertNotEqual(self.x, variables(3)[0])

    def test_leading_order_distinction(self):
        p = self.x*self.y + self.y**3
        self.assertEqual(leading_term(p, "lex")[0], (1, 1))
        self.assertEqual(leading_term(p, "grlex")[0], (0, 3))

    def test_reject_floating_coefficients_and_bad_exponents(self):
        with self.assertRaises(TypeError):
            MultiPolynomial({(1, 0): 0.1}, nvars=2)
        with self.assertRaises(ValueError):
            MultiPolynomial({(-1, 0): 1}, nvars=2)
        with self.assertRaises(ValueError):
            MultiPolynomial({(0,): 1}, nvars=2)
        with self.assertRaises(ValueError):
            self.x ** -1


class GrobnerBasisTests(unittest.TestCase):
    def setUp(self):
        self.x, self.y = variables(2)

    def test_buchberger_known_zero_dimensional_ideal(self):
        # x²=y, xy=1 implies x=y² and y³=1 over Q.
        basis = grobner_basis((self.x**2-self.y, self.x*self.y-1), order="lex")
        self.assertEqual(set(basis), {self.x-self.y**2, self.y**3-1})
        for f in (self.x**2-self.y, self.x*self.y-1):
            self.assertEqual(normal_form(f, basis), MultiPolynomial.zero(2))
        self.assertEqual(normal_form(self.x+self.y, basis), self.y**2+self.y)
        for i, f in enumerate(basis):
            for g in basis[i+1:]:
                self.assertFalse(normal_form(s_polynomial(f, g), basis))

    def test_s_polynomial_exposes_nonbasis(self):
        f, g = self.x**2-self.y, self.x*self.y-1
        self.assertTrue(normal_form(s_polynomial(f, g), (f, g)))

    def test_single_generator_and_unit_ideal(self):
        self.assertEqual(grobner_basis((2*self.x + 2,), order="lex"), (self.x+1,))
        self.assertEqual(grobner_basis((self.x+1, 5,), order="lex"), (MultiPolynomial.one(2),))
        self.assertEqual(normal_form(self.y**4+7, (MultiPolynomial.one(2),)), MultiPolynomial.zero(2))

    def test_zero_and_mixed_variables_rejected(self):
        with self.assertRaises(ValueError):
            grobner_basis((self.x-self.x,))
        with self.assertRaises(ValueError):
            grobner_basis((self.x, variables(3)[0]))
        with self.assertRaises(ValueError):
            normal_form(self.x, (self.y,), order="unknown")

    def test_reduction_respects_remainder_not_pseudoequality(self):
        basis = grobner_basis((self.x*self.y-1, self.y**2-1))
        self.assertEqual(normal_form(self.x*self.y-1, basis), MultiPolynomial.zero(2))
        self.assertNotEqual(normal_form(self.x, basis), MultiPolynomial.zero(2))

    def test_three_variable_criterion_and_generator_membership(self):
        x, y, z = variables(3)
        generators = (x*z-y, y*z-1, x+y-1)
        basis = grobner_basis(generators, max_pairs=200)
        for f in generators:
            self.assertEqual(normal_form(f, basis), MultiPolynomial.zero(3))
        for i, f in enumerate(basis):
            for g in basis[i+1:]:
                self.assertEqual(normal_form(s_polynomial(f, g), basis), MultiPolynomial.zero(3))

    def test_grlex_basis_and_reduced_leading_terms(self):
        generators = (self.x**2-self.y, self.x*self.y-1)
        basis = grobner_basis(generators, order="grlex")
        for f in generators:
            self.assertFalse(normal_form(f, basis, order="grlex"))
        for i, f in enumerate(basis):
            fm, fc = leading_term(f, order="grlex")
            self.assertEqual(fc, 1)
            for g in basis[i+1:]:
                self.assertFalse(normal_form(s_polynomial(f, g, order="grlex"), basis, order="grlex"))
            for j, g in enumerate(basis):
                if i == j:
                    continue
                gm = leading_term(g, order="grlex")[0]
                self.assertFalse(any(all(a <= b for a,b in zip(gm, monomial))
                                     for monomial, _ in f.terms if monomial != fm))

    def test_explicit_resource_budget(self):
        with self.assertRaises(ResourceLimitError):
            grobner_basis((self.x**2-self.y, self.x*self.y-1), max_pairs=1)
        with self.assertRaises(ValueError):
            grobner_basis((self.x,), max_pairs=0)
        with self.assertRaises(ResourceLimitError):
            grobner_basis((self.x**2-self.y, self.x*self.y-1), max_basis=2)

    def test_exact_division_step_budget_boundary(self):
        self.assertEqual(normal_form(self.x, (self.x,), max_steps=1), MultiPolynomial.zero(2))
        with self.assertRaises(ResourceLimitError):
            normal_form(self.x+self.y, (self.x,), max_steps=1)


if __name__ == "__main__":
    unittest.main()
