"""Adversarial invariants for transforms, independently derived from definitions.

These tests intentionally avoid calling implementation helpers or using its
intermediate results as the expected oracle.
"""

import cmath
import math
import unittest

from calc.transforms import (
    dft, idft, circular_convolution, linear_convolution, z_transform_finite,
)


class HostileFourierContracts(unittest.TestCase):
    def assertComplexSequenceNear(self, observed, expected, tolerance=2e-11):
        self.assertEqual(len(observed), len(expected))
        for index, (actual, oracle) in enumerate(zip(observed, expected)):
            self.assertLessEqual(abs(actual - oracle), tolerance,
                                 f"coefficient {index}: {actual!r} != {oracle!r}")

    def test_length_one_and_non_power_of_two_round_trip(self):
        self.assertEqual(dft([3+2j]), (3+2j,))
        self.assertEqual(idft([3+2j]), (3+2j,))
        values = (1-2j, -3+5j, 0.25, 2j, -0.75+0.5j, 3, -4j)
        self.assertComplexSequenceNear(idft(dft(values)), values)

    def test_shift_theorem_has_negative_phase(self):
        # Circular right shift by m changes each DFT bin by exp(-2*pi*i*k*m/N).
        values = [1+2j, 3-1j, -2, 0.5j, -4+3j]
        n, shift = len(values), 2
        shifted = values[-shift:] + values[:-shift]
        original = dft(values)
        expected = [original[k] * cmath.exp(-2j * math.pi * k * shift / n)
                    for k in range(n)]
        self.assertComplexSequenceNear(dft(shifted), expected)

    def test_modulation_theorem_and_parseval_for_complex_samples(self):
        x = [2+1j, -1-3j, 0.5+2j, -4, 3-1j, -0.75j]
        n, shift = len(x), 2
        modulated = [value * cmath.exp(2j * math.pi * shift * j / n)
                     for j, value in enumerate(x)]
        spectrum = dft(x)
        expected = [spectrum[(k - shift) % n] for k in range(n)]
        self.assertComplexSequenceNear(dft(modulated), expected)
        self.assertAlmostEqual(sum(abs(value)**2 for value in x),
                               sum(abs(value)**2 for value in spectrum) / n,
                               places=10)

    def test_dc_cancellation_is_not_naive_summation(self):
        # The actual DC bin is 1.0, despite terms spanning 32 decades.
        values = [1e16, 1.0, -1e16, 0.0]
        self.assertEqual(dft(values)[0], 1+0j)


class HostileConvolutionContracts(unittest.TestCase):
    def assertNear(self, actual, expected, tolerance=1e-11):
        self.assertEqual(len(actual), len(expected))
        for index, (a, e) in enumerate(zip(actual, expected)):
            self.assertLessEqual(abs(a-e), tolerance, f"at index {index}")

    def test_impulse_identity_and_linear_commutativity(self):
        signal = (3+2j, -4, 0.5j, 7)
        self.assertEqual(circular_convolution(signal, (1,0,0,0)), signal)
        b = (2-1j, -3+2j)
        self.assertNear(linear_convolution(signal,b), linear_convolution(b,signal))
        self.assertNear(linear_convolution(signal,(1,)), signal)

    def test_linear_convolution_equals_zero_padded_circular_convolution(self):
        left, right = (1+2j, -3, 2j), (2, 4-1j, -1, 0.5j)
        length = len(left)+len(right)-1
        padded_left = left + (0,)*(length-len(left))
        padded_right = right + (0,)*(length-len(right))
        self.assertNear(linear_convolution(left,right),
                        circular_convolution(padded_left,padded_right))

    def test_convolution_associativity_for_small_complex_signals(self):
        a=(1+1j,2); b=(-2j,3); c=(4,-1j)
        self.assertNear(linear_convolution(linear_convolution(a,b),c),
                        linear_convolution(a,linear_convolution(b,c)))


class HostileFiniteZTransformContracts(unittest.TestCase):
    def test_complex_z_and_negative_real_z(self):
        # 3-4/z+5/z^2 at z=-2 = 3+2+5/4.
        self.assertEqual(z_transform_finite((3,-4,5),-2), 6.25+0j)
        # 1+2/z+3/z^2 evaluated at i equals -2-2i.
        self.assertAlmostEqual(z_transform_finite((1,2,3),1j).real, -2)
        self.assertAlmostEqual(z_transform_finite((1,2,3),1j).imag, -2)

    def test_trailing_zeros_must_not_force_spurious_overflow(self):
        # The finite polynomial in z^-1 is identically constant 1 here.
        # An enormous negative exponent with zero coefficient is not a term.
        self.assertEqual(z_transform_finite((1,0,0),1e-300), 1+0j)
        self.assertEqual(z_transform_finite((0,0,0,0),1e-300), 0+0j)

    def test_small_coefficients_can_cancel_exponent_range(self):
        # z^-2 alone is outside float range, but x[2]*z^-2 is about 1e100.
        # All nonzero summands and the exact finite sum are representable.
        actual = z_transform_finite((1,0,1e-300),1e-200)
        self.assertLess(abs(actual - (1e100+0j))/1e100, 1e-12)
        complex_actual = z_transform_finite((1,0,1e-300),1e-200j)
        self.assertLess(abs(complex_actual - (-1e100+0j))/1e100, 1e-12)

    def test_properly_rejects_nonzero_overflow_and_invalid_samples(self):
        with self.assertRaises(ValueError):
            z_transform_finite((1,2,3), 1e-300)
        for invalid in ([], [float('nan')], [True], ['1']):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                z_transform_finite(invalid, 2)


if __name__ == '__main__':
    unittest.main()
