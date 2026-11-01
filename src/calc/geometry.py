"""Small numerical differential-geometry and differential-forms toolkit.

Every derivative is a centered *finite difference*, not symbolic or automatic
differentiation. The metric need only be smooth, symmetric, and nondegenerate
at the evaluation point; positive definiteness is not required.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
import math

Vector = tuple[float, ...]
Matrix = tuple[Vector, ...]


def _vector(values: Sequence[float], *, name: str = "vector") -> Vector:
    try:
        result = tuple(float(value) for value in values)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite numeric vector") from exc
    if not result or any(not math.isfinite(x) for x in result):
        raise ValueError(f"{name} must be a nonempty finite numeric vector")
    return result


def _matrix(values: Sequence[Sequence[float]], n: int, *, symmetric: bool = False) -> Matrix:
    try:
        matrix = tuple(_vector(row, name="matrix row") for row in values)
    except (TypeError, ValueError) as exc:
        raise ValueError("matrix must contain finite numeric rows") from exc
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError(f"expected {n} by {n} matrix")
    if symmetric:
        for i in range(n):
            for j in range(i):
                if not math.isclose(matrix[i][j], matrix[j][i], rel_tol=1e-9, abs_tol=1e-12):
                    raise ValueError("metric must be symmetric")
    return matrix


def _step(h: float, coordinate: float) -> float:
    try:
        value = float(h)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("step must be a positive finite number") from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError("step must be a positive finite number")
    delta = value * max(1.0, abs(coordinate))
    if not math.isfinite(delta) or coordinate + delta == coordinate:
        raise ValueError("step is not representable at this coordinate")
    return delta


def _shift(point: Vector, index: int, delta: float) -> Vector:
    result = list(point)
    result[index] += delta
    return tuple(result)


def _inverse(values: Matrix) -> Matrix:
    """Gauss-Jordan matrix inversion, with a relative pivot threshold."""
    n = len(values)
    aug = [list(row) + [float(i == j) for j in range(n)] for i, row in enumerate(values)]
    scale = max(abs(value) for row in values for value in row)
    if not scale or not math.isfinite(scale):
        raise ValueError("metric is degenerate")
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        pivot = aug[pivot_row][col]
        if abs(pivot) <= 1e-13 * scale:
            raise ValueError("metric is degenerate or ill-conditioned")
        aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
        factor = aug[col][col]
        aug[col] = [value / factor for value in aug[col]]
        for row in range(n):
            if row == col:
                continue
            weight = aug[row][col]
            aug[row] = [a - weight * b for a, b in zip(aug[row], aug[col])]
    result = tuple(tuple(row[n:]) for row in aug)
    if any(not math.isfinite(value) for row in result for value in row):
        raise ValueError("metric inverse is nonfinite")
    return result


def metric_inner(
    metric: Callable[[Vector], Sequence[Sequence[float]]],
    point: Sequence[float],
    left: Sequence[float],
    right: Sequence[float],
) -> float:
    """Calculate ``g_ij(x) left^i right^j`` in a coordinate basis."""
    x, u, v = _vector(point, name="point"), _vector(left, name="left"), _vector(right, name="right")
    n = len(x)
    if len(u) != n or len(v) != n:
        raise ValueError("tangent vectors must have the same dimension as point")
    g = _matrix(metric(x), n, symmetric=True)
    return sum(u[i] * g[i][j] * v[j] for i in range(n) for j in range(n))


def coordinate_jacobian(
    transform: Callable[[Vector], Sequence[float]],
    point: Sequence[float],
    *,
    h: float = 1e-5,
) -> Matrix:
    """Numerically calculate ``J[a][i] = ∂ transform[a] / ∂ x[i]``."""
    x = _vector(point, name="point")
    y = _vector(transform(x), name="transformed point")
    rows = [[0.0] * len(x) for _ in y]
    for i, coordinate in enumerate(x):
        delta = _step(h, coordinate)
        plus = _vector(transform(_shift(x, i, delta)), name="transformed point")
        minus = _vector(transform(_shift(x, i, -delta)), name="transformed point")
        if len(plus) != len(y) or len(minus) != len(y):
            raise ValueError("transform output dimension must be constant")
        for a in range(len(y)):
            rows[a][i] = (plus[a] - minus[a]) / (2.0 * delta)
    return tuple(tuple(row) for row in rows)


def pullback_metric(
    transform: Callable[[Vector], Sequence[float]],
    target_metric: Callable[[Vector], Sequence[Sequence[float]]],
    point: Sequence[float],
    *,
    h: float = 1e-5,
) -> Matrix:
    """Calculate ``J^T g(transform(x)) J`` for a coordinate map."""
    x = _vector(point, name="point")
    y = _vector(transform(x), name="transformed point")
    g = _matrix(target_metric(y), len(y), symmetric=True)
    jac = coordinate_jacobian(transform, x, h=h)
    n = len(x)
    return tuple(
        tuple(sum(jac[a][i] * g[a][b] * jac[b][j]
                  for a in range(len(y)) for b in range(len(y)))
              for j in range(n))
        for i in range(n)
    )


def christoffel_symbols(
    metric: Callable[[Vector], Sequence[Sequence[float]]],
    point: Sequence[float],
    *,
    h: float = 1e-5,
) -> tuple[Matrix, ...]:
    """Compute Levi-Civita ``Gamma[upper][lower1][lower2]``.

    Only valid where the metric is smooth, symmetric, nondegenerate and the
    chosen coordinates are regular. ``h`` is a relative finite-difference
    scale; there is no certified error bound.
    """
    x = _vector(point, name="point")
    n = len(x)
    g = _matrix(metric(x), n, symmetric=True)
    inverse = _inverse(g)
    dg = [[[0.0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    # dg[coordinate][row][col] = derivative of metric[row][col]
    for coordinate, position in enumerate(x):
        delta = _step(h, position)
        plus = _matrix(metric(_shift(x, coordinate, delta)), n, symmetric=True)
        minus = _matrix(metric(_shift(x, coordinate, -delta)), n, symmetric=True)
        for i in range(n):
            for j in range(n):
                dg[coordinate][i][j] = (plus[i][j] - minus[i][j]) / (2.0 * delta)
    result = tuple(
        tuple(
            tuple(0.5 * sum(inverse[upper][ell] *
                            (dg[j][ell][k] + dg[k][ell][j] - dg[ell][j][k])
                            for ell in range(n)) for k in range(n))
            for j in range(n)
        )
        for upper in range(n)
    )
    if any(not math.isfinite(value) for plane in result for row in plane for value in row):
        raise ValueError("connection coefficients are nonfinite")
    return result


def exterior_derivative_0(
    field: Callable[[Vector], float],
    point: Sequence[float],
    *,
    h: float = 1e-5,
) -> Vector:
    """Coordinate components of the 1-form ``df`` for a scalar field."""
    x = _vector(point, name="point")
    result = []
    for i, position in enumerate(x):
        delta = _step(h, position)
        result.append((float(field(_shift(x, i, delta))) -
                       float(field(_shift(x, i, -delta)))) / (2.0 * delta))
    return _vector(result, name="d field")


def exterior_derivative_1(
    form: Callable[[Vector], Sequence[float]],
    point: Sequence[float],
    *,
    h: float = 1e-5,
) -> dict[tuple[int, int], float]:
    """Two-form components ``(dω)_{ij}=∂_i ω_j-∂_j ω_i`` for ``i<j``.

    Coordinates use zero-based axis indexes. The result has one entry for
    each strictly increasing pair and is the empty dict in dimension one.
    """
    x = _vector(point, name="point")
    n = len(x)
    if len(_vector(form(x), name="one-form")) != n:
        raise ValueError("one-form needs one component per coordinate")
    derivatives: list[Vector] = []
    for i, position in enumerate(x):
        delta = _step(h, position)
        plus = _vector(form(_shift(x, i, delta)), name="one-form")
        minus = _vector(form(_shift(x, i, -delta)), name="one-form")
        if len(plus) != n or len(minus) != n:
            raise ValueError("one-form output dimension must be constant")
        derivatives.append(tuple((plus[j] - minus[j]) / (2.0 * delta) for j in range(n)))
    return {(i, j): derivatives[i][j] - derivatives[j][i]
            for i in range(n) for j in range(i + 1, n)}
