"""Bounded real matrix-calculus operations with explicit precision contracts.

Exact arithmetic uses ``Fraction`` for integer/rational entries. Supplying any
float selects finite binary floating-point arithmetic; mixed input therefore
does not silently promise exactness. Numerical derivatives are *estimates*,
not certified gradients or general automatic differentiation.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from fractions import Fraction
import math

Exact = int | Fraction
RealScalar = Exact | float
Matrix = tuple[tuple[Fraction | float, ...], ...]
FloatMatrix = tuple[tuple[float, ...], ...]


def _matrix(values: Sequence[Sequence[RealScalar]], name: str) -> Matrix:
    try:
        rows = tuple(tuple(row) for row in values)
    except TypeError as exc:
        raise TypeError(f"{name} must be a rectangular numeric matrix") from exc
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{name} must be a nonempty rectangular matrix")
    has_float = False
    for row in rows:
        for value in row:
            if type(value) is float:
                if not math.isfinite(value):
                    raise ValueError(f"{name} entries must be finite")
                has_float = True
            elif type(value) is int or isinstance(value, Fraction):
                continue
            else:
                raise TypeError(f"{name} entries must be int, Fraction or finite float")
    if has_float:
        try:
            result = tuple(tuple(float(v) for v in row) for row in rows)
        except (ValueError, OverflowError) as exc:
            raise ValueError(f"{name} entries are outside finite float range") from exc
        if any(not math.isfinite(v) for row in result for v in row):
            raise ValueError(f"{name} entries are outside finite float range")
        return result
    return tuple(tuple(Fraction(v) for v in row) for row in rows)


def _float_matrix(values: Sequence[Sequence[RealScalar]], name: str) -> FloatMatrix:
    exact_or_float = _matrix(values, name)
    try:
        result = tuple(tuple(float(v) for v in row) for row in exact_or_float)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} cannot be represented as finite floats") from exc
    if any(not math.isfinite(v) for row in result for v in row):
        raise ValueError(f"{name} cannot be represented as finite floats")
    return result


def _finite_float(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a real number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be representable as a real float") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _sum_terms(values: Sequence[Fraction | float], exact: bool) -> Fraction | float:
    if exact:
        return sum(values, Fraction(0))
    try:
        result = math.fsum(float(v) for v in values)
    except (OverflowError, ValueError) as exc:
        raise ValueError("matrix arithmetic overflowed") from exc
    if not math.isfinite(result):
        raise ValueError("matrix arithmetic overflowed")
    return result


def _sum_products(pairs: Sequence[tuple[Fraction | float, Fraction | float]],
                  exact: bool) -> Fraction | float:
    """Sum products, preserving a finite answer when float terms overflow.

    Normal finite arithmetic uses ``math.fsum``. Only an exceptional float
    overflow triggers exact reconstruction of the supplied binary64 inputs.
    This does not turn ordinary inexact input values into symbolic values.
    """
    if exact:
        return sum((left * right for left, right in pairs), Fraction(0))
    terms = [float(left) * float(right) for left, right in pairs]
    if all(math.isfinite(term) for term in terms):
        try:
            result = math.fsum(terms)
        except (OverflowError, ValueError):
            pass
        else:
            if math.isfinite(result):
                return result
    total = sum((Fraction.from_float(float(left)) * Fraction.from_float(float(right))
                 for left, right in pairs), Fraction(0))
    try:
        value = float(total)
    except OverflowError as exc:
        raise ValueError("matrix arithmetic overflowed") from exc
    if not math.isfinite(value):
        raise ValueError("matrix arithmetic overflowed")
    return value


def _central_quotient(plus: float, minus: float, step: float, name: str) -> float:
    """Central-difference quotient with a binary64-exact overflow fallback."""
    numerator = plus - minus
    denominator = 2 * step
    if math.isfinite(numerator) and math.isfinite(denominator):
        value = numerator / denominator
        if math.isfinite(value):
            return value
    exact = ((Fraction.from_float(plus) - Fraction.from_float(minus)) /
             (2 * Fraction.from_float(step)))
    try:
        value = float(exact)
    except OverflowError as exc:
        raise ValueError(f"{name} is nonfinite") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} is nonfinite")
    return value


def frobenius_inner(
    left: Sequence[Sequence[RealScalar]], right: Sequence[Sequence[RealScalar]],
) -> Fraction | float:
    """Compute ``<A,B>_F = tr(A^T B) = Σ_ij A_ij B_ij`` over real matrices."""
    a, b = _matrix(left, "left"), _matrix(right, "right")
    if (len(a), len(a[0])) != (len(b), len(b[0])):
        raise ValueError("Frobenius inner product requires matching dimensions")
    exact = all(isinstance(v, Fraction) for row in a+b for v in row)
    return _sum_products([(a[i][j], b[i][j])
                          for i in range(len(a)) for j in range(len(a[0]))], exact)


def squared_frobenius_gradient(matrix: Sequence[Sequence[RealScalar]]) -> Matrix:
    """Exact gradient ``∇_X ||X||_F² = 2X`` for real matrices."""
    a = _matrix(matrix, "matrix")
    result = tuple(tuple(2 * value for value in row) for row in a)
    if any(isinstance(v, float) and not math.isfinite(v) for row in result for v in row):
        raise ValueError("squared Frobenius gradient overflowed")
    return result


def trace_axb_gradient(
    left: Sequence[Sequence[RealScalar]], right: Sequence[Sequence[RealScalar]],
) -> Matrix:
    """Gradient of ``F(X)=tr(A X B)`` with ``A:p×m, X:m×n, B:n×p``.

    Returns ``(B A)^T`` as an ``m×n`` matrix. The formula is valid for real
    matrices and the standard Frobenius differential pairing.
    """
    a, b = _matrix(left, "left"), _matrix(right, "right")
    p, m = len(a), len(a[0])
    n = len(b)
    if len(b[0]) != p:
        raise ValueError("A:p×m and B:n×p need matching trace dimension p")
    exact = all(isinstance(v, Fraction) for row in a+b for v in row)
    return tuple(tuple(_sum_products([(a[k][i], b[j][k]) for k in range(p)], exact)
                       for j in range(n)) for i in range(m))


def _step(step: float, scale: float) -> float:
    s = _finite_float(step, "step")
    if s <= 0:
        raise ValueError("step must be positive")
    h = s * max(1.0, scale)
    if not math.isfinite(h) or h == 0:
        raise ValueError("numerical step is not finite and representable")
    return h


def _shift_one(matrix: FloatMatrix, row: int, col: int, h: float) -> FloatMatrix:
    return tuple(tuple(v+h if (i,j)==(row,col) else v for j,v in enumerate(values))
                 for i,values in enumerate(matrix))


def matrix_gradient(
    function: Callable[[FloatMatrix], RealScalar],
    matrix: Sequence[Sequence[RealScalar]], *, step: float = 1e-5,
) -> FloatMatrix:
    """Central finite-difference ``∂f/∂X_ij`` for a real scalar matrix map.

    Each coordinate uses ``h = step*max(1, |X_ij|)``. The function receives an
    immutable tuple-of-tuples of finite floats. Smoothness is assumed; this
    routine offers no rigorous truncation/roundoff error bound.
    """
    if not callable(function):
        raise TypeError("function must be callable")
    x = _float_matrix(matrix, "matrix")
    s = _finite_float(step, "step")
    if s <= 0:
        raise ValueError("step must be positive")
    results: list[tuple[float, ...]] = []
    for i,row in enumerate(x):
        output = []
        for j,value in enumerate(row):
            h = _step(s, abs(value))
            if value+h == value or value-h == value:
                raise ValueError("matrix perturbation is not representable")
            if not math.isfinite(value+h) or not math.isfinite(value-h):
                raise ValueError("matrix perturbation leaves finite float range")
            plus = _shift_one(x, i, j, h)
            minus = _shift_one(x, i, j, -h)
            f_plus = _finite_float(function(plus), "function output")
            f_minus = _finite_float(function(minus), "function output")
            derivative = _central_quotient(f_plus, f_minus, h, "matrix gradient")
            output.append(derivative)
        results.append(tuple(output))
    return tuple(results)


def matrix_jvp(
    transform: Callable[[FloatMatrix], Sequence[Sequence[RealScalar]]],
    matrix: Sequence[Sequence[RealScalar]],
    tangent: Sequence[Sequence[RealScalar]], *, step: float = 1e-5,
) -> FloatMatrix:
    """Central finite-difference Jacobian-vector product ``D F(X)[H]``.

    Supports real matrix-to-matrix maps with a stable output shape. The step
    is scaled relative to the largest input/tangent coordinate, and any
    unrepresentable nonzero directional perturbation is rejected.
    """
    if not callable(transform):
        raise TypeError("transform must be callable")
    x = _float_matrix(matrix, "matrix")
    direction = _float_matrix(tangent, "tangent")
    if (len(x),len(x[0])) != (len(direction),len(direction[0])):
        raise ValueError("tangent dimensions must match matrix")
    x_scale = max(abs(v) for row in x for v in row)
    d_scale = max(abs(v) for row in direction for v in row)
    h = _step(step, x_scale) / max(1.0, d_scale)
    if not math.isfinite(h) or h == 0:
        raise ValueError("directional step is not representable")
    plus_rows = []
    minus_rows = []
    for xrow,drow in zip(x,direction):
        plus_row = []
        minus_row = []
        for value,delta in zip(xrow,drow):
            offset = h*delta
            if not math.isfinite(offset):
                raise ValueError("directional perturbation is nonfinite")
            plus,minus = value+offset,value-offset
            if not math.isfinite(plus) or not math.isfinite(minus):
                raise ValueError("directional perturbation is nonfinite")
            if delta and (plus==value or minus==value):
                raise ValueError("nonzero tangent perturbation is unrepresentable")
            plus_row.append(plus)
            minus_row.append(minus)
        plus_rows.append(tuple(plus_row))
        minus_rows.append(tuple(minus_row))
    plus_result = _float_matrix(transform(tuple(plus_rows)), "transform output")
    minus_result = _float_matrix(transform(tuple(minus_rows)), "transform output")
    if (len(plus_result),len(plus_result[0])) != (len(minus_result),len(minus_result[0])):
        raise ValueError("transform output shape changes under perturbation")
    values = tuple(tuple(_central_quotient(plus_result[i][j], minus_result[i][j],
                                          h, "matrix Jacobian-vector product")
                         for j in range(len(plus_result[0]))) for i in range(len(plus_result)))
    if any(not math.isfinite(v) for row in values for v in row):
        raise ValueError("matrix Jacobian-vector product is nonfinite")
    return values
