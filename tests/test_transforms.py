"""Analytic and independently computed DFT/convolution witnesses."""

import cmath
import math
import unittest

from calc.transforms import circular_convolution, dft, idft, linear_convolution, z_transform_finite


class DiscreteFourierTests(unittest.TestCase):
    def assertSequenceComplexAlmostEqual(self, actual, expected, tol=1e-11):
        self.assertEqual(len(actual), len(expected))
        for a, b in zip(actual, expected):
            self.assertLessEqual(abs(a - b), tol, msg=f"{a!r} != {b!r}")

    def test_known_four_point_transform(self):
        # Independent evaluation of roots of unity 1, -i, -1, i.
        self.assertSequenceComplexAlmostEqual(
            dft([1, 2, 3, 4]),
            (10, -2 + 2j, -2, -2 - 2j),
        )

    def test_constant_and_impulse(self):
        self.assertSequenceComplexAlmostEqual(dft([3, 3, 3, 3]), (12, 0, 0, 0))
        self.assertSequenceComplexAlmostEqual(dft([1, 0, 0, 0]), (1, 1, 1, 1))

    def test_complex_frequency_basis(self):
        n, frequency = 9, 3
        signal = [cmath.exp(2j * math.pi * frequency * j / n) for j in range(n)]
        expected = [n if k == frequency else 0 for k in range(n)]
        self.assertSequenceComplexAlmostEqual(dft(signal), expected, tol=1e-10)

    def test_complex_roundtrip(self):
        values = (2 + 3j, -4 + 2j, 0, 0.5 - 1j, -2.5)
        self.assertSequenceComplexAlmostEqual(idft(dft(values)), values)

    def test_parseval_energy(self):
        samples = [1, 2, 3, 4]
        transformed = dft(samples)
        self.assertAlmostEqual(sum(abs(x)**2 for x in samples),
                               sum(abs(x)**2 for x in transformed) / len(samples),
                               places=11)

    def test_invalid_inputs(self):
        for invalid in ([], [1, float("nan")], [float("inf")], [complex(0, math.inf)], [True, 2], ["1", 2]):
            with self.subTest(invalid=invalid), self.assertRaises((TypeError, ValueError)):
                dft(invalid)
        with self.assertRaises((TypeError, ValueError)):
            idft([1, complex(float("nan"), 0)])
        with self.assertRaises(ValueError):
            dft([0] * 4097)
        with self.assertRaises(ValueError):
            dft([1e308] * 5)


class ConvolutionTests(unittest.TestCase):
    def assertSequenceComplexAlmostEqual(self, actual, expected, tol=1e-10):
        self.assertEqual(len(actual), len(expected))
        for a, b in zip(actual, expected):
            self.assertLessEqual(abs(a - b), tol, msg=f"{a!r} != {b!r}")

    def test_hand_computed_circular_convolution(self):
        self.assertSequenceComplexAlmostEqual(
            circular_convolution([1, 2, 3], [4, 5, 6]), [31, 31, 28],
        )

    def test_hand_computed_linear_convolution(self):
        self.assertSequenceComplexAlmostEqual(
            linear_convolution([1, 2, 3], [4, 5]), [4, 13, 22, 15],
        )

    def test_convolution_theorem(self):
        x, y = [1, -2, 3, 0.5], [2j, 3, -1j, 4]
        left = dft(circular_convolution(x, y))
        dx, dy = dft(x), dft(y)
        self.assertSequenceComplexAlmostEqual(left, [a * b for a, b in zip(dx, dy)])

    def test_invalid_shapes_and_values(self):
        with self.assertRaises(ValueError):
            circular_convolution([1, 2], [3])
        with self.assertRaises(ValueError):
            linear_convolution([], [2])
        with self.assertRaises(TypeError):
            linear_convolution([1, None], [2])

    def test_finite_z_transform(self):
        # 1 + 2 z^-1 + 3 z^-2 at z=2 gives 1 + 1 + 3/4.
        self.assertAlmostEqual(z_transform_finite([1, 2, 3], 2), 2.75)
        self.assertAlmostEqual(z_transform_finite([1j, 1], 1j), 0)
        with self.assertRaises(ValueError):
            z_transform_finite([1, 2], 0)
        with self.assertRaises(ValueError):
            z_transform_finite([1, 2, 3], 1e-300)


if __name__ == "__main__":
    unittest.main()
