"""Independent exact witnesses for bounded distributional calculus."""

from fractions import Fraction as Q
import itertools
import unittest

from calc.distributions import DiracMeasure, PiecewiseConstant


class StepDistributionTests(unittest.TestCase):
    def test_indicator_derivative_has_opposite_boundary_atoms(self):
        step = PiecewiseConstant([(1, 3, 1)])
        self.assertEqual(step.weak_derivative().atoms, ((Q(1), Q(1)), (Q(3), Q(-1))))
        self.assertEqual(step.weak_derivative().total_mass(), 0)

    def test_negative_step_and_non_unit_rational_height(self):
        step = PiecewiseConstant([(Q(-3, 2), Q(1, 2), Q(-5, 3))])
        self.assertEqual(step.integral(), Q(-10, 3))
        self.assertEqual(step.weak_derivative().atoms,
                         ((Q(-3, 2), Q(-5, 3)), (Q(1, 2), Q(5, 3))))

    def test_touching_equal_segments_are_canonicalized(self):
        step = PiecewiseConstant([(0, 1, 2), (1, 2, 2), (2, 3, 2)])
        self.assertEqual(step.intervals, ((Q(0), Q(3), Q(2)),))
        self.assertEqual(step.weak_derivative().atoms, ((Q(0), Q(2)), (Q(3), Q(-2))))

    def test_touching_unequal_segments_have_differences_at_jumps(self):
        step = PiecewiseConstant([(0, 2, 3), (2, 3, -1), (3, 5, 4)])
        self.assertEqual(step.weak_derivative().atoms,
                         ((Q(0), Q(3)), (Q(2), Q(-4)), (Q(3), Q(5)), (Q(5), Q(-4))))

    def test_gap_creates_two_separate_boundaries(self):
        step = PiecewiseConstant([(-2, -1, 2), (1, 2, -3)])
        self.assertEqual(step.weak_derivative().atoms,
                         ((Q(-2), Q(2)), (Q(-1), Q(-2)),
                          (Q(1), Q(-3)), (Q(2), Q(3))))
        self.assertEqual(step.integral(), -1)

    def test_half_open_pointwise_values(self):
        step = PiecewiseConstant([(0, 1, 7), (2, 4, 9)])
        self.assertEqual([step.value_at(x) for x in (-1, 0, Q(1, 2), 1, 2, 4)],
                         [0, 7, 7, 0, 9, 0])

    def test_exact_polynomial_pairing_and_fundamental_identity(self):
        step = PiecewiseConstant([(0, 2, 3)])
        self.assertEqual(step.pair_polynomial([0, 0, 1]), 8)  # ∫_0² 3x² dx
        self.assertEqual(step.weak_derivative().pair_polynomial([0, 1]), -6)
        self.assertEqual(step.pair_polynomial([1]), step.integral())

    def test_integration_by_parts_many_coefficients(self):
        step = PiecewiseConstant([(-3, -1, Q(1, 3)), (0, 2, -5), (2, 8, Q(7, 9))])
        for coefficients in ([1], [1, 2], [3, Q(2, 7), -1, 4, Q(-7, 8)], []):
            self.assertEqual(step.integration_by_parts_residual(coefficients), Q(0))

    def test_integrate_derivative_constant_germ(self):
        for step in (PiecewiseConstant(), PiecewiseConstant([(1, 2, 1)]),
                     PiecewiseConstant([(-2, 0, -5), (0, 3, 3)])):
            self.assertEqual(step.weak_derivative().pair_polynomial([1]), 0)

    def test_translation_commutes_with_derivative(self):
        step = PiecewiseConstant([(-2, -1, -3), (0, 4, Q(1, 7))])
        self.assertEqual(step.shift(Q(5, 4)).weak_derivative().atoms,
                         step.weak_derivative().shift(Q(5, 4)).atoms)

    def test_scaling_commutes_with_derivative(self):
        step = PiecewiseConstant([(0, 1, 4), (1, 2, -5)])
        self.assertEqual(step.scale(Q(-3, 2)).weak_derivative().atoms,
                         step.weak_derivative().scale(Q(-3, 2)).atoms)
        self.assertEqual(step.scale(0).weak_derivative().atoms, ())

    def test_atom_coalescing_order_and_variation(self):
        mu = DiracMeasure([(2, 3), (1, -2), (2, -4), (1, 2), (3, 0), (5, Q(1, 2))])
        self.assertEqual(mu.atoms, ((Q(2), Q(-1)), (Q(5), Q(1, 2))))
        self.assertEqual(mu.total_mass(), Q(-1, 2))
        self.assertEqual(mu.total_variation(), Q(3, 2))
        self.assertEqual(mu.pair_polynomial([1, 2, 3]), Q(-17) + Q(43, 1))

    def test_empty_zero_function_and_measure(self):
        step = PiecewiseConstant()
        self.assertEqual(step.integral(), 0)
        self.assertEqual(step.value_at(0), 0)
        self.assertEqual(step.weak_derivative().atoms, ())
        self.assertEqual(step.pair_polynomial([2, 3]), 0)
        self.assertEqual(DiracMeasure().pair_polynomial([1, 2]), 0)

    def test_unsorted_overlap_zero_width_rejected(self):
        for bad in ([(0, 0, 1)], [(2, 1, 1)], [(0, 2, 1), (1, 3, 1)],
                    [(2, 3, 1), (0, 1, 1)]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                PiecewiseConstant(bad)

    def test_unsupported_float_bool_and_malformed_structure(self):
        for bad in ([(0.0, 1, 1)], [(0, 1, True)], [(0, 1, 1.5)], [(0, 1)],
                    [2], [(0, 1, 1, 2)]):
            with self.subTest(bad=bad), self.assertRaises((TypeError, ValueError)):
                PiecewiseConstant(bad)
        for bad in ([(0.0, 1)], [(2, True)], [(1,)], [(1, 2, 3)]):
            with self.subTest(bad=bad), self.assertRaises((TypeError, ValueError)):
                DiracMeasure(bad)

    def test_invalid_polynomial_and_shift_rejected(self):
        step = PiecewiseConstant([(0, 1, 1)])
        for bad in ([0.25], [True], [1] * 66):
            with self.subTest(bad=str(bad)[:20]), self.assertRaises((TypeError, ValueError)):
                step.pair_polynomial(bad)
        with self.assertRaises(TypeError):
            step.shift(1.0)
        with self.assertRaises(TypeError):
            step.scale(True)

    def test_finite_resource_bounds(self):
        with self.assertRaises(ValueError):
            PiecewiseConstant((i, i + Q(1, 2), 1) for i in range(2049))
        with self.assertRaises(ValueError):
            DiracMeasure((i, 1) for i in range(4097))
        with self.assertRaises(ValueError):
            PiecewiseConstant([(0, 1, 2**4097)])
        with self.assertRaises(ValueError):
            PiecewiseConstant([(0, 1, 1)]).pair_polynomial(itertools.repeat(1))


if __name__ == "__main__":
    unittest.main()
