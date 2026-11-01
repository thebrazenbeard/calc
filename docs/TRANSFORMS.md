# Discrete transforms: Fourier, convolution, and finite z-transform

**Status:** implemented as narrow, auditable examples in `src/calc/transforms.py`. This is not a general integral-transform package. The code is original, dependency-free Python with deterministic tests in `tests/test_transforms.py`.

## 1. Discrete Fourier transform

Given a nonempty sequence of length `N`, the forward transform in this repository is

\[
X_k=\sum_{n=0}^{N-1}x_n\,e^{-2\pi i kn/N},\quad k=0,\ldots,N-1.
\]

The inverse uses the opposite sign and normalization

\[
x_n=\frac{1}{N}\sum_{k=0}^{N-1}X_k\,e^{+2\pi i kn/N}.
\]

This is the same default normalization convention used by NumPy and SciPy. The output is always complex. For `x=(1,2,3,4)`, the exact forward result is `(10, -2+2i, -2, -2-2i)`; `idft(dft(x))` recovers `x` up to floating-point roundoff. Parseval's identity under this convention is

\[
\sum_n|x_n|^2=\frac{1}{N}\sum_k|X_k|^2.
\]

The implementation uses the direct `O(N²)` sum with compensated real and imaginary accumulation, not an FFT. It reduces integer phases modulo `N` before computing trigonometric functions. Inputs are limited to 4096 finite samples to avoid accidentally requesting extremely expensive quadratic work. Floating-point accuracy is not certified for poorly conditioned data.

## 2. Convolution

For equal-length sequences, `circular_convolution(a,b)` implements

\[
c_k=\sum_{j=0}^{N-1}a_jb_{(k-j)\bmod N}.
\]

The convolution theorem here is `dft(circular_convolution(a,b)) == dft(a)*dft(b)`, within numerical roundoff. For `a=(1,2,3)` and `b=(4,5,6)`, the circular result is `(31,31,28)`.

`linear_convolution` instead assumes zero outside each supplied sequence and produces `N+M-1` outputs. For `(1,2,3)` convolved with `(4,5)`, the result is `(4,13,22,15)`. The two operations must not be conflated: circular convolution aliases linear-convolution tails unless sufficient zero-padding is used.

## 3. Finite unilateral z-transform

`z_transform_finite(x,z)` evaluates only the finite sum

\[
X(z)=\sum_{n=0}^{N-1}x_n z^{-n},\qquad z\ne0.
\]

For `x=(1,2,3)`, `z=2`, the value is `1+1+3/4=11/4`. This does not infer the z-transform of an *infinite* sequence, a rational generating function, poles, or a region of convergence. Very small nonzero `z` can overflow the finite range and is rejected.

## Usage

```python
from calc.transforms import dft, idft, circular_convolution, z_transform_finite

assert abs(dft([1, 2, 3, 4])[1] - (-2 + 2j)) < 1e-12
assert abs(idft(dft([1, 2, 3]))[2] - 3) < 1e-12
assert circular_convolution([1, 2], [3, 4]) == (11+0j, 10+0j)
assert z_transform_finite([1, 2, 3], 2) == 2.75+0j
```

## Limits and future work

- No continuous Fourier, Laplace, Mellin, or general inverse integral transform; no analytic integration or distribution theory.
- No fast Fourier transform, frequency-axis generation, nonuniform sampling, windowing, or spectral-leakage correction.
- Finite-input complex arithmetic only; no arbitrary precision or validated error bound.
- Further development should add stable FFT algorithms or bridge to NumPy/SciPy with clearly versioned optional dependencies, not silently relabel this direct DFT as an FFT.

## Primary references (links, not copied content)

1. [NumPy: Discrete Fourier Transform, definitions and normalization](https://numpy.org/doc/stable/reference/routines.fft.html).
2. [SciPy: Discrete Fourier transforms and inverse algorithms](https://docs.scipy.org/doc/scipy/reference/fft.html).
3. [SciPy: Fourier transform tutorial](https://docs.scipy.org/doc/scipy/tutorial/fft.html).

The mathematics and API were written independently. NumPy/SciPy are references for conventions and production-grade implementations; this repository does not vendor their code.
