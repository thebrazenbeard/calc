"""Small dependency-free discrete transforms and convolution.

Conventions match NumPy/SciPy's default one-dimensional DFT: the forward
transform is unnormalized with a negative complex exponent; the inverse uses
the positive exponent and divides by N. This is a direct O(N**2) DFT for
education and small signals, *not* a fast Fourier transform or a continuous
Fourier/Laplace integral. Finite complex inputs are required.
"""

from __future__ import annotations

import math
from numbers import Number
from typing import Sequence

__all__ = [
    "dft", "idft", "circular_convolution", "linear_convolution",
    "z_transform_finite",
]

_MAX_SAMPLES = 4096


def _finite_complex(value: object, name: str) -> complex:
    if isinstance(value, bool) or not isinstance(value, Number):
        raise TypeError(f"{name} must be a finite real or complex number")
    try:
        result = complex(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{name} cannot be represented as a finite complex number") from exc
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError(f"{name} must be finite")
    return result


def _samples(values: Sequence[Number], name: str) -> tuple[complex, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError(f"{name} must be a finite sequence of numbers")
    n = len(values)
    if n == 0:
        raise ValueError(f"{name} must be nonempty")
    if n > _MAX_SAMPLES:
        raise ValueError(f"{name} is limited to {_MAX_SAMPLES} samples for O(N²) algorithms")
    return tuple(_finite_complex(value, f"{name}[{i}]") for i, value in enumerate(values))


def _sum_complex(terms: Sequence[complex]) -> complex:
    """Use compensated real/imaginary summation and reject overflow."""
    try:
        value = complex(
            math.fsum(term.real for term in terms),
            math.fsum(term.imag for term in terms),
        )
    except OverflowError as exc:
        raise ValueError("result overflows finite floating-point arithmetic") from exc
    return _finite_complex(value, "result")


def _dft(values: Sequence[Number], *, inverse: bool) -> tuple[complex, ...]:
    x = _samples(values, "values")
    n = len(x)
    sign = 1 if inverse else -1
    scale = 1 / n if inverse else 1
    output: list[complex] = []
    for k in range(n):
        terms: list[complex] = []
        for j, sample in enumerate(x):
            # Reduce the integer phase modulo N before forming its float;
            # this avoids avoidable argument-reduction error for large j*k.
            theta = math.tau * ((j * k) % n) / n
            unit = complex(math.cos(theta), sign * math.sin(theta))
            terms.append(sample * unit)
        output.append(_finite_complex(scale * _sum_complex(terms), "transform output"))
    return tuple(output)


def dft(values: Sequence[Number]) -> tuple[complex, ...]:
    """Return ``X[k] = Σ_j x[j] exp(-2πi jk/N)`` for N nonempty samples.

    The result is complex, even for real input. Cost is O(N²) time and O(N)
    working memory; use a dedicated FFT package for production-scale signals.
    """
    return _dft(values, inverse=False)


def idft(values: Sequence[Number]) -> tuple[complex, ...]:
    """Return ``x[j] = (1/N) Σ_k X[k] exp(+2πi jk/N)``."""
    return _dft(values, inverse=True)


def circular_convolution(
    left: Sequence[Number], right: Sequence[Number],
) -> tuple[complex, ...]:
    """Length-N cyclic convolution ``c[k]=Σ_j a[j] b[(k-j) mod N]``.

    Both arrays must have the same positive length. No zero-padding is done;
    this is *not* ordinary linear convolution.
    """
    a = _samples(left, "left")
    b = _samples(right, "right")
    if len(a) != len(b):
        raise ValueError("circular convolution requires equal input lengths")
    n = len(a)
    return tuple(_sum_complex([a[j] * b[(k - j) % n] for j in range(n)]) for k in range(n))


def linear_convolution(
    left: Sequence[Number], right: Sequence[Number],
) -> tuple[complex, ...]:
    """Finite linear convolution with zero extension, length N+M-1."""
    a = _samples(left, "left")
    b = _samples(right, "right")
    n, m = len(a), len(b)
    result = []
    for k in range(n + m - 1):
        start, stop = max(0, k - m + 1), min(n, k + 1)
        result.append(_sum_complex([a[j] * b[k-j] for j in range(start, stop)]))
    return tuple(result)


def z_transform_finite(values: Sequence[Number], z: Number) -> complex:
    """Evaluate the finite unilateral z-transform ``Σ_{n=0}^{N-1} x[n] z^-n``.

    Only the supplied finite sequence is represented. This does not infer a
    generating function or region of convergence for an infinite sequence.
    """
    x = _samples(values, "values")
    z = _finite_complex(z, "z")
    if z == 0:
        raise ValueError("z must be nonzero for this finite z-transform convention")
    terms: list[complex] = []
    for i, value in enumerate(x):
        if value == 0:
            # A mathematically zero term must not evaluate a potentially
            # overflowing inverse power just to multiply it by zero.
            continue
        try:
            term = value * z ** (-i)
        except (ZeroDivisionError, OverflowError):
            term = complex(math.nan, 0)
        if not (math.isfinite(term.real) and math.isfinite(term.imag)):
            # z**(-i) may overflow before a tiny coefficient cancels its
            # magnitude. Dividing the coefficient progressively avoids that
            # false failure when the final term is representable.
            term = value
            for _ in range(i):
                try:
                    term /= z
                except (ZeroDivisionError, OverflowError) as exc:
                    raise ValueError("z-transform term exceeds finite floating-point range") from exc
                if not (math.isfinite(term.real) and math.isfinite(term.imag)):
                    raise ValueError("z-transform term exceeds finite floating-point range")
        terms.append(_finite_complex(term, f"z-transform term {i}"))
    return _sum_complex(terms)
