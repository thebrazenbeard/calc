"""Independent exact-rational examples for bounded affine algebraic geometry."""

from fractions import Fraction
from itertools import repeat
from random import Random
import unittest

from calc.affine_geometry import AffineSystem, partial_derivative
from calc.grobner import MultiPolynomial, ResourceLimitError, variables


class AffineGeometryTests(unittest.TestCase):
    def test_circle_integer_grid_is_four_axis_points(self):
        x, y = variables(2)
        circle = AffineSystem((x*x + y*y - 1,))
        points = circle.points_on_grid(((-2, -1, 0, 1, 2), (-2, -1, 0, 1, 2)))
        expected = ((Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1)),
                    (Fraction(0), Fraction(1)), (Fraction(1), Fraction(0)))
        self.assertEqual(points, expected)
        self.assertFalse(circle.contains((Fraction(3, 5), Fraction(3, 5))))

    def test_parabola_has_exact_rational_points(self):
        x, y = variables(2)
        curve = AffineSystem((y-x*x,))
        grid = ((-2, -1, 0, 1, 2), (0, 1, 4))
        self.assertEqual(curve.points_on_grid(grid),
                         ((Fraction(-2), Fraction(4)),
                          (Fraction(-1), Fraction(1)),
                          (Fraction(0), Fraction(0)),
                          (Fraction(1), Fraction(1)),
                          (Fraction(2), Fraction(4))))

    def test_grid_search_is_not_an_emptiness_oracle(self):
        x, = variables(1)
        system = AffineSystem((2*x-1,))
        self.assertEqual(system.points_on_grid(((0, 1, 2),)), ())
        self.assertTrue(system.contains((Fraction(1, 2),)))
        self.assertEqual(system.points_on_grid(((Fraction(1, 2),),)),
                         ((Fraction(1, 2),),))

    def test_grid_preserves_input_order_and_deduplicates(self):
        x, y = variables(2)
        system = AffineSystem((x-y,))
        self.assertEqual(system.points_on_grid(((2, 1, 2, 0), (1, 2, 0, 1))),
                         ((Fraction(2), Fraction(2)),
                          (Fraction(1), Fraction(1)),
                          (Fraction(0), Fraction(0))))

    def test_zero_equation_system_is_full_rational_affine_space(self):
        zero = MultiPolynomial.zero(2)
        system = AffineSystem((zero,))
        self.assertTrue(system.contains((7, Fraction(-3, 4))))
        self.assertEqual(len(system.points_on_grid(((0, 1), (2, 3)))), 4)
        self.assertFalse(system.ideal_contains(MultiPolynomial.one(2)))
        self.assertTrue(system.ideal_contains(zero))

    def test_unit_ideal_is_empty_zero_set_and_contains_every_polynomial(self):
        x, y = variables(2)
        system = AffineSystem((x, MultiPolynomial.one(2)))
        self.assertFalse(system.contains((0, 0)))
        self.assertEqual(system.points_on_grid(((0, 1), (0, 1))), ())
        self.assertTrue(system.ideal_contains(y**2 + 7))

    def test_grobner_ideal_membership_not_just_zero_set_vanishing(self):
        x, = variables(1)
        double_origin = AffineSystem((x*x,))
        self.assertTrue(double_origin.contains((0,)))
        self.assertTrue(double_origin.ideal_contains(x**3))
        self.assertFalse(double_origin.ideal_contains(x))
        self.assertFalse(double_origin.ideal_contains(1 + x))

    def test_nonlinear_two_variable_elimination_and_order(self):
        x, y = variables(2)
        system = AffineSystem((x-y, y*y-1))
        self.assertTrue(system.ideal_contains(x*x-1, order="lex"))
        self.assertTrue(system.ideal_contains(x*x-1, order="grlex"))
        self.assertFalse(system.ideal_contains(x-1))
        self.assertEqual(system.points_on_grid(((-1, 0, 1), (-1, 0, 1))),
                         ((Fraction(-1), Fraction(-1)),
                          (Fraction(1), Fraction(1))))

    def test_multivariate_ideal_closure(self):
        x, y = variables(2)
        f = x*x-y
        system = AffineSystem((f, y*y-1))
        self.assertTrue(system.ideal_contains(f*(x+y+2)))
        self.assertTrue(system.ideal_contains(x**4-1))
        self.assertFalse(system.ideal_contains(x+1))

    def test_linear_derivatives_independent_oracles(self):
        x, y = variables(2)
        f = 7*x*x*y + 3*x*y*y + 5
        self.assertEqual(partial_derivative(f, 0), 14*x*y + 3*y*y)
        self.assertEqual(partial_derivative(f, 1), 7*x*x + 6*x*y)
        self.assertEqual(partial_derivative(f, 0).evaluate((2, 3)), Fraction(111))

    def test_jacobian_of_circle_and_cusp(self):
        x, y = variables(2)
        circle = AffineSystem((x*x+y*y-1,))
        self.assertEqual(circle.jacobian_at((Fraction(3, 5), Fraction(4, 5))),
                         ((Fraction(6, 5), Fraction(8, 5)),))
        cusp = AffineSystem((y*y-x**3,))
        self.assertEqual(cusp.jacobian_at((0, 0)), ((Fraction(0), Fraction(0)),))
        self.assertEqual(cusp.jacobian_at((1, 1)), ((Fraction(-3), Fraction(2)),))

    def test_seeded_exact_linear_solution_witnesses(self):
        rng = Random(21048)
        x, y = variables(2)
        for _ in range(70):
            a = rng.choice((-4, -3, -2, -1, 1, 2, 3, 4))
            b = rng.choice((-4, -3, -2, -1, 1, 2, 3, 4))
            if a+b == 0:
                continue
            c = rng.randint(-15, 15)
            answer = Fraction(c, a+b)
            system = AffineSystem((a*x+b*y-c, x-y))
            self.assertTrue(system.contains((answer, answer)))
            self.assertEqual(system.points_on_grid(((answer, answer+1), (answer+1, answer))),
                             ((answer, answer),))
            self.assertFalse(system.contains((answer+1, answer)))

    def test_reject_inexact_coordinates_and_wrong_shape(self):
        x, y = variables(2)
        system = AffineSystem((x-y,))
        with self.assertRaises(TypeError):
            system.contains((0.0, 0))
        with self.assertRaises(TypeError):
            system.contains((True, 0))
        with self.assertRaises(ValueError):
            system.contains((0,))
        with self.assertRaises(ValueError):
            system.jacobian_at((0, 0, 0))
        with self.assertRaises(ValueError):
            system.points_on_grid(((0,),))
        with self.assertRaises(TypeError):
            system.points_on_grid(((0,), range(3)))

    def test_reject_bad_polynomials_and_order(self):
        x, y = variables(2)
        z, = variables(1)
        with self.assertRaises(TypeError):
            AffineSystem(("x",))
        with self.assertRaises(ValueError):
            AffineSystem(())
        with self.assertRaises(ValueError):
            AffineSystem((x, z))
        with self.assertRaises(ValueError):
            AffineSystem((x,)).ideal_contains(z)
        with self.assertRaises(ValueError):
            AffineSystem((x,)).ideal_contains(x, order="invalid")
        with self.assertRaises(TypeError):
            partial_derivative("x", 0)
        with self.assertRaises(ValueError):
            partial_derivative(x, True)
        with self.assertRaises(ValueError):
            partial_derivative(y, 2)

    def test_resource_budgets_and_bounded_iterables(self):
        x, y = variables(2)
        with self.assertRaises(ResourceLimitError):
            AffineSystem(repeat(x))
        with self.assertRaises(ResourceLimitError):
            AffineSystem((x**33,))
        with self.assertRaises(ResourceLimitError):
            AffineSystem((variables(5)[0],))
        with self.assertRaises(ResourceLimitError):
            AffineSystem((x-y,)).contains((1 << 513, 0))
        with self.assertRaises(ResourceLimitError):
            AffineSystem((x-y,)).points_on_grid(((0,) * 65, (0,)))
        z = variables(3)[2]
        with self.assertRaises(ResourceLimitError):
            AffineSystem((z,)).points_on_grid((tuple(range(22)),) * 3)

    def test_empty_grid_axis_is_an_empty_search_not_a_theorem(self):
        x, y = variables(2)
        system = AffineSystem((x+y,))
        self.assertEqual(system.points_on_grid(((), (0, 1))), ())


if __name__ == "__main__":
    unittest.main()
