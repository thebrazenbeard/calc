"""Independent exact witnesses for sequence-based discrete calculus."""

from fractions import Fraction
import random
import unittest

from calc.discrete import (
    backward_differences,
    central_differences,
    cumulative_sum,
    discrete_integral,
    forward_differences,
    power_sum,
    summation_by_parts,
)


class DifferenceTests(unittest.TestCase):
    def test_forward_polynomial(self):
        values = [0, 1, 4, 9, 16]
        self.assertEqual(forward_differences(values), (1, 3, 5, 7))
        self.assertEqual(forward_differences(values, order=2), (2, 2, 2))
        self.assertEqual(forward_differences(values, order=3), (0, 0))

    def test_step_normalization(self):
        # f(x) = x^2 sampled at x = 0,2,4,6.
        values = (0, 4, 16, 36)
        self.assertEqual(forward_differences(values, step=2), (2, 6, 10))
        self.assertEqual(forward_differences(values, order=2, step=2), (2, 2))

    def test_descending_grid_and_backward_alignment(self):
        # f(x)=x**2 sampled at x=4,3,2,1,0: h=-1.
        descending = (16, 9, 4, 1, 0)
        self.assertEqual(forward_differences(descending, step=-1), (7, 5, 3, 1))
        # Backward at x=3,2,1,0 and forward at x=4,3,2,1.
        self.assertEqual(backward_differences(descending, step=-1), (7, 5, 3, 1))
        self.assertEqual(forward_differences(descending, order=2, step=-1), (2, 2, 2))

    def test_backward_and_central_stencils(self):
        values = (0, 1, 4, 9, 16)
        self.assertEqual(backward_differences(values), (1, 3, 5, 7))
        self.assertEqual(central_differences(values), (2, 4, 6))

    def test_exact_rational(self):
        q = Fraction(1, 3)
        self.assertEqual(forward_differences((q, 2*q, 4*q), step=q), (1, 2))
        self.assertEqual(central_differences((q, 2*q, 4*q), step=q),
                         (Fraction(3, 2),))

    def test_zero_order_and_short_sequence(self):
        self.assertEqual(forward_differences((3, 5), order=0), (3, 5))
        self.assertEqual(forward_differences((3, 5), order=3), ())
        self.assertEqual(central_differences((3, 5)), ())


class IntegralTests(unittest.TestCase):
    def test_antiderivative_and_signed_bounds(self):
        values = (2, -1, 3, 0)
        self.assertEqual(cumulative_sum(values, initial=5), (5, 7, 6, 9, 9))
        self.assertEqual(discrete_integral(values, 1, 3), 2)
        self.assertEqual(discrete_integral(values, 3, 1), -2)
        self.assertEqual(discrete_integral(values, 2, 2), 0)

    def test_empty_and_exact_fractional_sums(self):
        self.assertEqual(cumulative_sum(()), (0,))
        self.assertEqual(discrete_integral((), 0, 0), 0)
        values = (Fraction(1, 2), Fraction(-2, 3), Fraction(1, 6))
        self.assertEqual(discrete_integral(values, 0, 3), 0)

    def test_discrete_fundamental_theorem_seeded(self):
        rng = random.Random(617)
        for _ in range(50):
            seq = tuple(rng.randrange(-10, 11) for _ in range(rng.randrange(0, 15)))
            prefix = cumulative_sum(seq, initial=-7)
            self.assertEqual(tuple(prefix[i+1]-prefix[i] for i in range(len(seq))), seq)
            for i in range(len(seq)+1):
                self.assertEqual(discrete_integral(seq, 0, i), prefix[i]+7)

    def test_summation_by_parts_known(self):
        # Sum u_k (v_{k+1}-v_k) = [uv] - Sum v_{k+1}(u_{k+1}-u_k)
        lhs, boundary, correction = summation_by_parts((1, 2, 4), (0, 3, 5))
        self.assertEqual((lhs, boundary, correction), (7, 20, 13))
        self.assertEqual(lhs, boundary-correction)

    def test_summation_by_parts_seeded(self):
        rng = random.Random(904)
        for n in range(10):
            u = tuple(Fraction(rng.randrange(-9, 10), 3) for _ in range(n+1))
            v = tuple(Fraction(rng.randrange(-9, 10), 5) for _ in range(n+1))
            lhs, boundary, correction = summation_by_parts(u, v)
            self.assertEqual(lhs, boundary-correction)


class FaulhaberTests(unittest.TestCase):
    def test_known_power_sums(self):
        self.assertEqual(power_sum(5, 0), 5)
        self.assertEqual(power_sum(5, 1), 10)
        self.assertEqual(power_sum(5, 2), 30)
        self.assertEqual(power_sum(5, 3), 100)
        self.assertEqual(power_sum(4, 4), 98)
        self.assertEqual(power_sum(6, 5), 4425)

    def test_brute_force_independent_oracle(self):
        for n in range(13):
            for exponent in range(11):
                with self.subTest(n=n, exponent=exponent):
                    self.assertEqual(power_sum(n, exponent), sum(k**exponent for k in range(n)))

    def test_zero_and_large_index(self):
        self.assertEqual(power_sum(0, 10), 0)
        n = 10**8
        self.assertEqual(power_sum(n, 1), n*(n-1)//2)

    def test_large_polynomial_and_telescoping_witness(self):
        for n in (0, 1, 2, 100, 10**9, 10**12-1):
            self.assertEqual(power_sum(n, 2), n*(n-1)*(2*n-1)//6)
            self.assertEqual(power_sum(n, 3), (n*(n-1)//2)**2)
            for p in range(11):
                with self.subTest(n=n, power=p):
                    self.assertEqual(power_sum(n+1, p)-power_sum(n, p), n**p)


class InvalidInputs(unittest.TestCase):
    def test_not_exact_or_nonfinite(self):
        with self.assertRaises(TypeError):
            forward_differences((1.5, 2.5))
        with self.assertRaises(TypeError):
            cumulative_sum((True, 2))
        with self.assertRaises(TypeError):
            power_sum(True, 2)

    def test_stencil_bounds_and_zero_step(self):
        for call in (lambda: forward_differences((1, 2), step=0),
                     lambda: backward_differences((1,), order=-1),
                     lambda: central_differences((1, 2, 3), step=0),
                     lambda: forward_differences((1, 2), order=100)):
            with self.assertRaises(ValueError):
                call()

    def test_integral_bad_bounds_and_incompatible_parts(self):
        with self.assertRaises(IndexError):
            discrete_integral((1, 2), 0, 3)
        with self.assertRaises(ValueError):
            summation_by_parts((1, 2), (1,))

    def test_resource_caps(self):
        with self.assertRaises(ValueError):
            cumulative_sum((1,) * 10001)
        with self.assertRaises(ValueError):
            power_sum(10**13, 1)
        with self.assertRaises(ValueError):
            power_sum(4, 11)


if __name__ == "__main__":
    unittest.main()
