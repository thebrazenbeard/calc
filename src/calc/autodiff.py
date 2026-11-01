"""Small, explicit forward-mode automatic differentiation over real scalars.

``Dual(value, tangent)`` represents ``value + tangent * ε`` where ε² = 0.
Propagating the tangent through each operator applies the chain rule exactly
at the level of its floating-point operations. It is *not* symbolic calculus,
and does not differentiate opaque ``math.*`` calls or arbitrary Python code.

The public multivariate interface is ``f(tuple_of_duals)``. ``gradient`` and
``jacobian`` evaluate the function once per input dimension, so callables
should be deterministic and side-effect-free.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import cos as _cos
from math import exp as _exp
from math import isfinite, log as _log
from math import sin as _sin
from math import sqrt as _sqrt
from numbers import Real
from types import NotImplementedType


def _real(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(f"{name} must be finite") from error
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True, slots=True)
class Dual:
    """Real-valued dual number carrying one directional derivative."""

    value: float
    tangent: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _real(self.value, "value"))
        object.__setattr__(self, "tangent", _real(self.tangent, "tangent"))

    @staticmethod
    def _coerce(other: object) -> Dual | NotImplementedType:
        if isinstance(other, Dual):
            return other
        if isinstance(other, Real) and not isinstance(other, bool):
            return Dual(other, 0)
        return NotImplemented

    def __add__(self, other: object) -> Dual:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        return Dual(self.value + rhs.value, self.tangent + rhs.tangent)

    __radd__ = __add__

    def __neg__(self) -> Dual:
        return Dual(-self.value, -self.tangent)

    def __pos__(self) -> Dual:
        return self

    def __sub__(self, other: object) -> Dual:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self + (-rhs)

    def __rsub__(self, other: object) -> Dual:
        lhs = self._coerce(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs - self

    def __mul__(self, other: object) -> Dual:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        return Dual(self.value * rhs.value,
                    self.tangent * rhs.value + self.value * rhs.tangent)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> Dual:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        if rhs.value == 0:
            raise ZeroDivisionError("dual division by zero")
        return Dual(self.value / rhs.value,
                    (self.tangent * rhs.value - self.value * rhs.tangent)
                    / (rhs.value * rhs.value))

    def __rtruediv__(self, other: object) -> Dual:
        lhs = self._coerce(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs / self

    def __pow__(self, exponent: object) -> Dual:
        exp_dual = self._coerce(exponent)
        if exp_dual is NotImplemented:
            return NotImplemented
        n = exp_dual.value
        x = self.value
        # Real-valued x**y with a varying exponent requires positive x.
        if exp_dual.tangent != 0:
            if x <= 0:
                raise ValueError("variable exponent requires positive base")
            result = x**n
            return Dual(result, result * (exp_dual.tangent * _log(x)
                                          + n * self.tangent / x))
        if x == 0:
            if n < 0:
                raise ZeroDivisionError("zero raised to a negative exponent")
            if n == 0:
                raise ValueError("0**0 is undefined")
            if not n.is_integer() and n != 1:
                raise ValueError("fractional powers are not two-sided real differentiable at zero")
        if x < 0 and not n.is_integer():
            raise ValueError("noninteger powers of negative numbers are not real")
        result = x**n
        # Avoid 0**(-1) in the derivative of x**0 at nonzero x.
        tangent = 0 if n == 0 else n * x**(n - 1) * self.tangent
        return Dual(result, tangent)

    def __rpow__(self, other: object) -> Dual:
        lhs = self._coerce(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs**self


def sin(x: Dual | Real) -> Dual:
    """Sine with propagated tangent (input in radians)."""
    d = Dual._coerce(x)
    if d is NotImplemented:
        raise TypeError("sin expects a dual or real number")
    return Dual(_sin(d.value), _cos(d.value) * d.tangent)


def cos(x: Dual | Real) -> Dual:
    """Cosine with propagated tangent (input in radians)."""
    d = Dual._coerce(x)
    if d is NotImplemented:
        raise TypeError("cos expects a dual or real number")
    return Dual(_cos(d.value), -_sin(d.value) * d.tangent)


def exp(x: Dual | Real) -> Dual:
    """Base-e exponential with propagated tangent."""
    d = Dual._coerce(x)
    if d is NotImplemented:
        raise TypeError("exp expects a dual or real number")
    y = _exp(d.value)
    return Dual(y, y * d.tangent)


def log(x: Dual | Real) -> Dual:
    """Natural logarithm; domain x > 0."""
    d = Dual._coerce(x)
    if d is NotImplemented:
        raise TypeError("log expects a dual or real number")
    if d.value <= 0:
        raise ValueError("log requires a positive real argument")
    return Dual(_log(d.value), d.tangent / d.value)


def sqrt(x: Dual | Real) -> Dual:
    """Principal square root with finite derivative; domain x > 0."""
    d = Dual._coerce(x)
    if d is NotImplemented:
        raise TypeError("sqrt expects a dual or real number")
    if d.value <= 0:
        raise ValueError("sqrt derivative is not finite at nonpositive argument")
    y = _sqrt(d.value)
    return Dual(y, d.tangent / (2 * y))


def _output(output: Dual | Real) -> Dual:
    d = Dual._coerce(output)
    if d is NotImplemented:
        raise TypeError("differentiated function must return a real number or Dual")
    return d


def value_and_derivative(f: Callable[[Dual], Dual | Real], x: Real) -> tuple[float, float]:
    """Return ``(f(x), f'(x))`` using one forward-mode evaluation."""
    result = _output(f(Dual(x, 1)))
    return (result.value, result.tangent)


def derivative(f: Callable[[Dual], Dual | Real], x: Real) -> float:
    """Compute derivative of a scalar-to-scalar dual-compatible function."""
    return value_and_derivative(f, x)[1]


def _point(values: Sequence[Real], name: str) -> tuple[float, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError(f"{name} must be a sequence of finite reals")
    result = tuple(_real(v, name) for v in values)
    if not result:
        raise ValueError(f"{name} must be nonempty")
    return result


def _seed(point: tuple[float, ...], tangent: tuple[float, ...]) -> tuple[Dual, ...]:
    return tuple(Dual(value, d) for value, d in zip(point, tangent))


def gradient(f: Callable[[tuple[Dual, ...]], Dual | Real],
             point: Sequence[Real]) -> tuple[float, ...]:
    """Gradient of a scalar-valued function, evaluated once per dimension."""
    values = _point(point, "point")
    return tuple(
        _output(f(_seed(values, tuple(float(i == axis) for i in range(len(values)))))).tangent
        for axis in range(len(values))
    )


def directional_derivative(f: Callable[[tuple[Dual, ...]], Dual | Real],
                           point: Sequence[Real], direction: Sequence[Real]) -> float:
    """One-pass directional derivative, not normalized by direction magnitude."""
    values = _point(point, "point")
    tangent = _point(direction, "direction")
    if len(values) != len(tangent):
        raise ValueError("point and direction dimensions must match")
    return _output(f(_seed(values, tangent))).tangent


def jacobian(f: Callable[[tuple[Dual, ...]], Sequence[Dual | Real]],
             point: Sequence[Real]) -> tuple[tuple[float, ...], ...]:
    """Rows of a vector-valued function's Jacobian; repeated evaluations."""
    values = _point(point, "point")
    columns = []
    output_count: int | None = None
    for axis in range(len(values)):
        result = f(_seed(values, tuple(float(i == axis) for i in range(len(values)))))
        if isinstance(result, (str, bytes)) or not isinstance(result, Sequence):
            raise TypeError("jacobian function must return a sequence of real or Dual values")
        out = tuple(_output(item) for item in result)
        if output_count is None:
            if not out:
                raise ValueError("jacobian function output must be nonempty")
            output_count = len(out)
        elif len(out) != output_count:
            raise ValueError("jacobian function output shape changed between evaluations")
        columns.append(tuple(item.tangent for item in out))
    return tuple(tuple(column[row] for column in columns) for row in range(output_count))


__all__ = [
    "Dual", "sin", "cos", "exp", "log", "sqrt", "value_and_derivative",
    "derivative", "gradient", "directional_derivative", "jacobian",
]
