"""Independent analytic and algebraic witnesses for wave-two modules.

These checks deliberately construct their oracles from external mathematical
identities rather than calling another implementation of the same algorithm.
"""

import math
import unittest

from calc import autodiff, geometry, stochastic, structures, symbolic


class WaveTwoContractTests(unittest.TestCase):
    def test_symbolic_product_chain_and_domain(self):
        expr = symbolic.parse_expression("x**2 * sin(x) + log(x)")
        derivative = expr.diff("x")
        for x in (0.25, 0.9, 2.0, 5.0):
            expected = 2*x*math.sin(x) + x*x*math.cos(x) + 1/x
            self.assertAlmostEqual(derivative.eval({"x": x}), expected, delta=2e-12)
        with self.assertRaises((ValueError, ZeroDivisionError)):
            derivative.eval({"x": 0.0})

    def test_symbolic_parser_must_not_execute_python(self):
        for source in (
            "__import__('os').system('true')", "x.__class__",
            "[x for x in (1,2)]", "(lambda x: x)(3)", "2**1000",
            "sin(x, y)", "x[0]", "f'{x}'",
        ):
            with self.subTest(source=source), self.assertRaises((ValueError, TypeError)):
                symbolic.parse_expression(source)

    def test_forward_dual_gradient_exact_chain(self):
        def f(values):
            x, y = values
            return autodiff.exp(x*y) + autodiff.sin(x)
        x, y = 0.4, -0.8
        expected = (y*math.exp(x*y)+math.cos(x), x*math.exp(x*y))
        actual = autodiff.gradient(f, (x, y))
        for a, b in zip(actual, expected):
            self.assertAlmostEqual(a, b, delta=1e-12)

    def test_pullback_of_affine_metric(self):
        def transform(point):
            x, y = point
            return (2*x+y, x-3*y)
        def euclidean(_point):
            return ((1, 0), (0, 1))
        # A = [[2,1],[1,-3]], so A^T A = [[5,-1],[-1,10]].
        result = geometry.pullback_metric(transform, euclidean, (0.2, -0.7))
        for row, expected_row in zip(result, ((5, -1), (-1, 10))):
            for actual, expected in zip(row, expected_row):
                self.assertAlmostEqual(actual, expected, delta=1e-9)

    def test_metric_inner_rejects_overflow_instead_of_returning_infinity(self):
        # All inputs are finite, but the mathematical result is outside the
        # representable floating-point range. An unmarked inf is misleading.
        with self.assertRaises((ValueError, OverflowError)):
            geometry.metric_inner(
                lambda _point: ((1.0,),), (0.0,), (1e308,), (1e308,)
            )

    def test_exterior_derivative_squared_vanishes_for_smooth_scalar(self):
        def f(point):
            x, y = point
            return x*x*y + math.sin(x*y)
        def one_form(point):
            return geometry.exterior_derivative_0(f, point, h=1e-4)
        two_form = geometry.exterior_derivative_1(one_form, (0.4, 0.7), h=1e-4)
        self.assertAlmostEqual(two_form[(0, 1)], 0.0, delta=1e-6)

    def test_stochastic_deterministic_drift_without_diffusion(self):
        path = stochastic.euler_maruyama(
            lambda t, y: 3, lambda t, y: 0,
            2, -1, 0.125, 8, increments=(0.7,)*8,
        )
        self.assertAlmostEqual(path[-1][0], 3.0)
        self.assertAlmostEqual(path[-1][1], 2.0)

    def test_stochastic_rejects_unrepresentable_time_progress(self):
        # Positive dt is not sufficient if t0+dt == t0 in binary64.
        # Returning the same timestamp with a changed state violates IVP time.
        with self.assertRaises(ValueError):
            stochastic.euler_maruyama(
                lambda t, y: 1.0, lambda t, y: 0.0,
                1e20, 0.0, 1e-9, 2, increments=(0.0, 0.0),
            )

    def test_prime_field_and_tropical_semiring_identities(self):
        field = structures.PrimeField(11)
        for a in range(11):
            for b in range(11):
                for c in range(11):
                    self.assertEqual(field.mul(a, field.add(b,c)),
                                     field.add(field.mul(a,b), field.mul(a,c)))
        tadd, tmul = structures.tropical_add, structures.tropical_multiply
        for a,b,c in ((1,2,3), (-7,4,2), (math.inf,8,-4)):
            self.assertEqual(tmul(a,tadd(b,c)),tadd(tmul(a,b),tmul(a,c)))


if __name__ == "__main__":
    unittest.main()
