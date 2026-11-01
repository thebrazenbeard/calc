"""Independent convex-analysis witnesses for bounded nonsmooth operations."""

from fractions import Fraction
import unittest

from calc.nonsmooth import (
    ConvexPiecewiseLinear,
    absolute_subgradient,
    hinge_subgradient,
    soft_threshold,
    prox_hinge,
    moreau_envelope_absolute,
)


class ConvexPiecewiseLinearTests(unittest.TestCase):
    def setUp(self):
        self.absolute = ConvexPiecewiseLinear(((1, 0), (-1, 0)))
        self.hinge = ConvexPiecewiseLinear(((0, 0), (1, 0)))

    def test_absolute_values_and_subgradients(self):
        self.assertEqual(self.absolute.value(-3), 3)
        self.assertEqual(self.absolute.value(0), 0)
        self.assertEqual(self.absolute.value(4), 4)
        self.assertEqual(self.absolute.subdifferential(-3), (Fraction(-1), Fraction(-1)))
        self.assertEqual(self.absolute.subdifferential(0), (Fraction(-1), Fraction(1)))
        self.assertEqual(self.absolute.subdifferential(4), (Fraction(1), Fraction(1)))

    def test_hinge_kink_and_smooth_points(self):
        self.assertEqual(self.hinge.subdifferential(-2), (Fraction(0), Fraction(0)))
        self.assertEqual(self.hinge.subdifferential(0), (Fraction(0), Fraction(1)))
        self.assertEqual(self.hinge.subdifferential(3), (Fraction(1), Fraction(1)))

    def test_intersection_of_nontrivial_affine_pieces(self):
        f = ConvexPiecewiseLinear(((2, -1), (-1, 2), (0, 0)))
        self.assertEqual(f.value(1), 1)
        self.assertEqual(f.subdifferential(1), (Fraction(-1), Fraction(2)))
        self.assertEqual(f.subdifferential(0), (Fraction(-1), Fraction(-1)))
        self.assertEqual(f.subdifferential(3), (Fraction(2), Fraction(2)))

    def test_exact_rational_kink(self):
        f = ConvexPiecewiseLinear(((1, 0), (-2, 1)))
        self.assertEqual(f.subdifferential(Fraction(1, 3)), (Fraction(-2), Fraction(1)))
        self.assertEqual(f.value(Fraction(1, 3)), Fraction(1, 3))

    def test_global_supporting_line_inequality(self):
        # g is a convex subgradient iff f(y) >= f(x)+g*(y-x) for all y.
        x = Fraction(0)
        for g in (Fraction(-1), Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(1)):
            for y in (Fraction(-7, 3), Fraction(-1), Fraction(0), Fraction(2, 3), Fraction(4)):
                self.assertGreaterEqual(self.absolute.value(y), self.absolute.value(x) + g*(y-x))
        # g=2 is outside the subgradient interval and fails at y=1.
        self.assertLess(self.absolute.value(1), self.absolute.value(x) + 2*(1-x))

    def test_directional_derivative_at_kink(self):
        f = ConvexPiecewiseLinear(((2, -1), (-1, 2)))
        self.assertEqual(f.directional_derivative(1, 3), 6)
        self.assertEqual(f.directional_derivative(1, -3), 3)
        self.assertEqual(f.directional_derivative(1, 0), 0)

    def test_clarke_equals_convex_subdifferential_on_this_class(self):
        self.assertEqual(self.absolute.clarke_subdifferential(0), self.absolute.subdifferential(0))
        self.assertEqual(self.hinge.clarke_subdifferential(0), self.hinge.subdifferential(0))

    def test_convexity_midpoint_inequality(self):
        f = ConvexPiecewiseLinear(((4, -2), (-1, 2), (0, 0)))
        for x in (-3, 0, 1, 5):
            for y in (-4, -1, 2, 7):
                middle = Fraction(x+y, 2)
                self.assertLessEqual(f.value(middle), (f.value(x) + f.value(y))/2)

    def test_invalid_affine_inputs_and_budgets(self):
        with self.assertRaises(ValueError):
            ConvexPiecewiseLinear(())
        # Arbitrary, potentially nonconvex callbacks are not an accepted API.
        with self.assertRaises(TypeError):
            ConvexPiecewiseLinear((lambda x: -abs(x),))
        with self.assertRaises((TypeError, ValueError)):
            ConvexPiecewiseLinear(((1,),))
        with self.assertRaises(TypeError):
            ConvexPiecewiseLinear(((1.0, 0),))
        with self.assertRaises(TypeError):
            self.absolute.value(True)
        with self.assertRaises(ValueError):
            ConvexPiecewiseLinear(((1, 0),)*1025)
        with self.assertRaises(ValueError):
            ConvexPiecewiseLinear(((1 << 5000, 0),))


class ProximalTests(unittest.TestCase):
    def test_absolute_and_hinge_subgradient_helpers(self):
        self.assertEqual(absolute_subgradient(0), (Fraction(-1), Fraction(1)))
        self.assertEqual(absolute_subgradient(-2), (Fraction(-1), Fraction(-1)))
        self.assertEqual(hinge_subgradient(0), (Fraction(0), Fraction(1)))
        self.assertEqual(hinge_subgradient(-2), (Fraction(0), Fraction(0)))

    def test_soft_threshold_closed_form(self):
        self.assertEqual(soft_threshold(-3, 2), -1)
        self.assertEqual(soft_threshold(1, 2), 0)
        self.assertEqual(soft_threshold(3, 2), 1)
        self.assertEqual(soft_threshold(Fraction(5, 3), Fraction(2, 3)), 1)
        self.assertEqual(soft_threshold(2, 0), 2)

    def test_hinge_prox_closed_form(self):
        self.assertEqual(prox_hinge(-2, 3), -2)
        self.assertEqual(prox_hinge(1, 3), 0)
        self.assertEqual(prox_hinge(5, 3), 2)
        self.assertEqual(prox_hinge(Fraction(1, 2), Fraction(1, 2)), 0)

    def test_moreau_envelope_absolute(self):
        self.assertEqual(moreau_envelope_absolute(Fraction(1, 2), 1), Fraction(1, 8))
        self.assertEqual(moreau_envelope_absolute(2, 1), Fraction(3, 2))
        self.assertEqual(moreau_envelope_absolute(-2, 1), Fraction(3, 2))

    def test_proximal_stationarity_and_global_minimum(self):
        for v, t in ((-4, 1), (-1, 2), (0, 1), (1, 2), (4, 1)):
            x = soft_threshold(v, t)
            lo, hi = absolute_subgradient(x)
            witness = Fraction(v-x, t)
            self.assertLessEqual(lo, witness)
            self.assertLessEqual(witness, hi)
            f = lambda y: abs(y) + (y-v)**2/(2*t)
            for y in range(-6, 7):
                self.assertLessEqual(f(x), f(Fraction(y)))

    def test_hinge_proximal_stationarity_and_global_minimum(self):
        for v, t in ((-3, 2), (1, 3), (5, 3), (3, 1)):
            x = prox_hinge(v, t)
            low, high = hinge_subgradient(x)
            stationarity = Fraction(v-x, t)
            self.assertLessEqual(low, stationarity)
            self.assertLessEqual(stationarity, high)
            objective = lambda y: t * max(Fraction(0), y) + (y-v)**2 / 2
            for y in range(-6, 8):
                self.assertLessEqual(objective(x), objective(Fraction(y)))

    def test_proximal_validation(self):
        for f in (soft_threshold, prox_hinge):
            with self.assertRaises(ValueError):
                f(1, -1)
            with self.assertRaises(TypeError):
                f(1, 0.1)
            with self.assertRaises(TypeError):
                f(False, 1)
        with self.assertRaises(ValueError):
            moreau_envelope_absolute(1, 0)


if __name__ == "__main__":
    unittest.main()
