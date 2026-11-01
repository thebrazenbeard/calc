"""Bounded Chebyshev--Lobatto polynomial interpolation and differentiation.

This is a numerical interpolation building block, not a PDE solver or a
general spectral-method framework. Degree N uses N+1 nodes on [-1, 1],
with N bounded to make the dense O(N^2) differentiation matrix tractable.
The returned binary64 values are approximations, not certified enclosures.
"""

from __future__ import annotations

from collections.abc import Iterable
from itertools import islice
import math
from numbers import Real

_MAX_DEGREE = 128


def _degree(value: object) -> int:
    if type(value) is not int:
        raise TypeError("Chebyshev degree must be an integer (not bool)")
    if not 1 <= value <= _MAX_DEGREE:
        raise ValueError(f"Chebyshev degree must be between 1 and {_MAX_DEGREE}")
    return value


def _real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} cannot be represented as a finite float") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _samples(values: Iterable[Real]) -> tuple[float, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError("node samples must be a finite numeric iterable")
    try:
        iterator = iter(values)
    except TypeError as exc:
        raise TypeError("node samples must be an iterable") from exc
    # Never fully consume an unbounded generator.
    samples = tuple(islice(iterator, _MAX_DEGREE + 2))
    if not 2 <= len(samples) <= _MAX_DEGREE + 1:
        raise ValueError("node samples must have between 2 and 129 entries")
    return tuple(_real(value, "node sample") for value in samples)


def chebyshev_lobatto_nodes(degree: int) -> tuple[float, ...]:
    """Return descending ``x_j=cos(pi*j/N)``, ``j=0..N``, on ``[-1,1]``.

    Endpoint values are set exactly to +/-1; all interior nodes use libm.
    """
    n = _degree(degree)
    return (1.0, *(0.0 if 2 * j == n else math.cos(math.pi * j / n)
                   for j in range(1, n)), -1.0)


def _barycentric_weights(degree: int) -> tuple[float, ...]:
    # Common scale cancels in the second barycentric formula and matrix rows.
    return tuple((0.5 if j in (0, degree) else 1.0) * (-1 if j % 2 else 1)
                 for j in range(degree + 1))


def barycentric_interpolate(values: Iterable[Real], x: Real) -> float:
    """Evaluate the degree-N interpolating polynomial at x in [-1,1].

    ``values[j]`` belongs to ``chebyshev_lobatto_nodes(N)[j]``. Interpolation
    is exact at the stored *floating-point nodes*, returning the supplied
    value unchanged. Elsewhere it uses scaled second-form barycentric sums.
    Near a node the result is approximate; no rigorous error bound is given.
    """
    samples = _samples(values)
    query = _real(x, "evaluation point")
    if not -1.0 <= query <= 1.0:
        raise ValueError("evaluation point must belong to [-1,1]")
    n = len(samples) - 1
    nodes = chebyshev_lobatto_nodes(n)
    for xi, yi in zip(nodes, samples):
        if query == xi:
            return yi
    weights = _barycentric_weights(n)
    terms = tuple(w / (query - xi) for w, xi in zip(weights, nodes))
    scale = max(abs(term) for term in terms)
    amplitude = max(abs(y) for y in samples)
    if amplitude == 0.0:
        return 0.0
    denominator = math.fsum(term / scale for term in terms)
    if denominator == 0.0:
        raise ValueError("barycentric denominator underflowed")
    numerator = math.fsum((term / scale) * (y / amplitude)
                          for term, y in zip(terms, samples))
    answer = (numerator / denominator) * amplitude
    if not math.isfinite(answer):
        raise ValueError("interpolated value is outside finite float range")
    return answer


def chebyshev_differentiation_matrix(degree: int) -> tuple[tuple[float, ...], ...]:
    """Return dense ``D`` with ``(D @ f)_i = p'(x_i)`` for nodal polynomial p.

    For i != j, ``D_ij = (w_j/w_i)/(x_i-x_j)`` with barycentric weights w;
    the diagonal is the negative sum of the other entries in its row. The
    order matches the descending Chebyshev-Lobatto nodes.
    """
    n = _degree(degree)
    nodes = chebyshev_lobatto_nodes(n)
    weights = _barycentric_weights(n)
    rows: list[tuple[float, ...]] = []
    for i in range(n + 1):
        row = [0.0] * (n + 1)
        for j in range(n + 1):
            if i != j:
                row[j] = (weights[j] / weights[i]) / (nodes[i] - nodes[j])
        row[i] = -math.fsum(row)
        rows.append(tuple(row))
    return tuple(rows)


def differentiate_at_nodes(values: Iterable[Real]) -> tuple[float, ...]:
    """Differentiate the interpolating polynomial at its Chebyshev nodes.

    Uses off-diagonal differences ``D_ij*(f_j-f_i)``. This form annihilates
    constant data exactly and avoids row-sum cancellation on near-constant
    data. Output remains ordinary binary64 and may overflow for large input.
    """
    samples = _samples(values)
    matrix = chebyshev_differentiation_matrix(len(samples) - 1)
    output = []
    for i, row in enumerate(matrix):
        try:
            derivative = math.fsum(row[j] * (samples[j] - samples[i])
                                   for j in range(len(samples)) if j != i)
        except (OverflowError, ValueError) as exc:
            raise ValueError("interpolant derivative overflowed") from exc
        if not math.isfinite(derivative):
            raise ValueError("interpolant derivative is outside finite float range")
        output.append(derivative)
    return tuple(output)


__all__ = [
    "chebyshev_lobatto_nodes", "barycentric_interpolate",
    "chebyshev_differentiation_matrix", "differentiate_at_nodes",
]
