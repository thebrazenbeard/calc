"""Independent exact witnesses for oriented Euclidean exterior algebra."""

from fractions import Fraction
from itertools import combinations
import unittest

from calc.forms import Form, basis_form, hodge_star, inner_product, norm_squared, wedge


class ExteriorAlgebraTests(unittest.TestCase):
    def test_3d_hodge_basis_signs(self):
        dx = basis_form(3, (0,))
        dy = basis_form(3, (1,))
        dz = basis_form(3, (2,))
        self.assertEqual(hodge_star(dx), basis_form(3, (1, 2)))
        self.assertEqual(hodge_star(dy), basis_form(3, (0, 2), -1))
        self.assertEqual(hodge_star(dz), basis_form(3, (0, 1)))
        self.assertEqual(hodge_star(wedge(dx, dy)), dz)
        self.assertEqual(hodge_star(Form(3, 0, {(): 1})), basis_form(3, (0, 1, 2)))
        self.assertEqual(hodge_star(basis_form(3, (0, 1, 2))), Form(3, 0, {(): 1}))

    def test_star_square_sign_across_dimensions(self):
        for n in range(1, 6):
            for k in range(n + 1):
                indices = tuple(range(k))
                f = basis_form(n, indices, Fraction(3, 7))
                expected = f if (k * (n - k)) % 2 == 0 else basis_form(n, indices, Fraction(-3, 7))
                with self.subTest(n=n, k=k):
                    self.assertEqual(hodge_star(hodge_star(f)), expected)

    def test_exhaustive_basis_hodge_pairing_and_orientation(self):
        for n in range(1, 6):
            volume = basis_form(n, tuple(range(n)))
            for k in range(n + 1):
                for basis in combinations(range(n), k):
                    one = basis_form(n, basis)
                    with self.subTest(n=n, basis=basis):
                        self.assertEqual(wedge(one, hodge_star(one)), volume)
                        self.assertEqual(wedge(one, hodge_star(one, orientation=-1)), -volume)
                        self.assertEqual(norm_squared(hodge_star(one)), Fraction(1))
                        self.assertEqual(hodge_star(hodge_star(one, orientation=-1), orientation=-1),
                                         ((-1) ** (k * (n - k))) * one)

    def test_exhaustive_graded_commutativity(self):
        for n in range(1, 5):
            for p in range(n + 1):
                for q in range(n + 1):
                    sign = (-1) ** (p * q)
                    for left_indices in combinations(range(n), p):
                        for right_indices in combinations(range(n), q):
                            a = basis_form(n, left_indices, Fraction(-2, 3))
                            b = basis_form(n, right_indices, 5)
                            with self.subTest(n=n, p=p, q=q, left=left_indices, right=right_indices):
                                self.assertEqual(wedge(a, b), sign * wedge(b, a))

    def test_wedge_sign_associativity_and_nilpotence(self):
        e0, e1, e2 = (basis_form(3, (i,)) for i in range(3))
        self.assertEqual(wedge(e0, e1), basis_form(3, (0, 1)))
        self.assertEqual(wedge(e1, e0), basis_form(3, (0, 1), -1))
        self.assertEqual(wedge(e0, e0), Form(3, 2, {}))
        self.assertEqual(wedge(wedge(e0, e1), e2), wedge(e0, wedge(e1, e2)))
        self.assertEqual(wedge(wedge(e2, e1), e0), basis_form(3, (0, 1, 2), -1))

    def test_grading_and_wedge_zero_above_top_degree(self):
        top = basis_form(3, (0, 1, 2))
        e0 = basis_form(3, (0,))
        product = wedge(top, e0)
        self.assertEqual(product, Form(3, 4, {}))
        with self.assertRaises(ValueError):
            hodge_star(product)

    def test_rational_bilinearity_and_inner_product(self):
        a = Form(3, 1, {(0,): Fraction(1, 2), (1,): 2})
        b = Form(3, 1, {(0,): 3, (1,): Fraction(-1, 3)})
        self.assertEqual(inner_product(a, b), Fraction(5, 6))
        self.assertEqual(norm_squared(a), Fraction(17, 4))
        self.assertEqual(inner_product(a, basis_form(3, (1, 2))), 0)
        self.assertEqual(wedge(a, basis_form(3, (2,))).as_dict(),
                         {(0, 2): Fraction(1, 2), (1, 2): Fraction(2)})

    def test_exact_linear_arithmetic_and_wedge_distributivity(self):
        a = Form(4, 1, {(0,): Fraction(2, 3), (1,): Fraction(-1, 2)})
        b = Form(4, 1, {(0,): Fraction(-2, 3), (2,): 3})
        c = Form(4, 2, {(2, 3): Fraction(5, 7), (0, 2): 1})
        self.assertEqual(a + b, Form(4, 1, {(1,): Fraction(-1, 2), (2,): 3}))
        self.assertEqual(a - b, a + (-b))
        self.assertEqual((a * Fraction(3, 2)) / 3, a / 2)
        self.assertEqual(wedge(a + b, c), wedge(a, c) + wedge(b, c))
        self.assertEqual(wedge(c, a + b), wedge(c, a) + wedge(c, b))
        self.assertEqual(wedge(wedge(a, b), c), wedge(a, wedge(b, c)))
        with self.assertRaises(ValueError):
            _ = a + c
        with self.assertRaises(ZeroDivisionError):
            _ = a / 0

    def test_hodge_pairing_identity_and_orientation(self):
        a = Form(3, 2, {(0, 1): 2, (1, 2): -3})
        b = Form(3, 2, {(0, 1): 5, (0, 2): 7, (1, 2): 11})
        volume = wedge(a, hodge_star(b))
        self.assertEqual(volume, Form(3, 3, {(0, 1, 2): Fraction(-23)}))
        self.assertEqual(inner_product(a, b), -23)
        self.assertEqual(wedge(a, hodge_star(b, orientation=-1)),
                         Form(3, 3, {(0, 1, 2): 23}))
        self.assertEqual(hodge_star(hodge_star(a, orientation=-1), orientation=-1), a)

    def test_input_is_immutable_and_normalized(self):
        incoming = {(0,): 2, (1,): Fraction(0)}
        f = Form(3, 1, incoming)
        incoming[(0,)] = 9
        self.assertEqual(f.as_dict(), {(0,): Fraction(2)})
        copy = f.as_dict()
        copy[(0,)] = 99
        self.assertEqual(f.as_dict(), {(0,): Fraction(2)})

    def test_invalid_dimension_or_form_coefficients(self):
        for constructor in (
            lambda: Form(0, 0, {(): 1}),
            lambda: Form(3, -1, {}),
            lambda: Form(3, 1, {(0, 1): 1}),
            lambda: Form(3, 1, {(2, 2): 1}),
            lambda: Form(3, 1, {(3,): 1}),
            lambda: Form(3, 1, {(True,): 1}),
            lambda: Form(3, 1, {(0,): 0.1}),
            lambda: Form(3, 7, {(0,): 1}),
        ):
            with self.subTest(constructor=constructor), self.assertRaises((TypeError, ValueError)):
                constructor()
        with self.assertRaises(ValueError):
            wedge(Form(2, 1, {(0,): 1}), Form(3, 1, {(1,): 1}))
        with self.assertRaises(ValueError):
            inner_product(Form(2, 0, {(): 1}), Form(3, 0, {(): 1}))
        with self.assertRaises(ValueError):
            hodge_star(Form(3, 1, {(0,): 1}), orientation=0)


if __name__ == "__main__":
    unittest.main()
