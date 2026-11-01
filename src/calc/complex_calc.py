"""Small numerical primitives for complex-variable calculus.

These functions **do not prove holomorphicity** or certify a contour integral.
The Cauchy–Riemann check is a local finite-difference diagnostic; numerical
quadrature cannot reliably detect arbitrary unsampled singularities.

All functions are dependency-free and accept user-supplied callables, which
execute normally with the caller's privileges. Do not pass untrusted code.
"""

from __future__ import annotations

import math
from numbers import Complex, Real
from typing import Callable


ComplexFunction = Callable[[complex], complex]
Path = Callable[[float], complex]


def _finite_real(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be real")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _finite_complex(value: Complex, name: str) -> complex:
    if isinstance(value, (bool, str, bytes)):
        raise ValueError(f"{name} must be a finite complex or real number")
    try:
        result = complex(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite complex or real number") from exc
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        raise ValueError(f"{name} must be finite")
    return result


def _step(point: complex, step: Real | None) -> float:
    if step is None:
        return 1e-5 * max(1.0, abs(point))
    result = _finite_real(step, "step")
    if result <= 0:
        raise ValueError("step must be positive")
    return result


def complex_derivative(
    function: ComplexFunction,
    point: Complex,
    *,
    direction: Complex = 1 + 0j,
    step: Real | None = None,
) -> complex:
    """Approximate the directional complex difference quotient.

    Calculates ``(f(z+h*d)-f(z-h*d))/(2*h*d)`` with unit direction ``d``.
    For a truly holomorphic function, the limit is direction-independent.
    A single numerical direction **is not** evidence of holomorphicity.
    """
    z = _finite_complex(point, "point")
    d = _finite_complex(direction, "direction")
    if d == 0:
        raise ValueError("direction must be nonzero")
    d /= abs(d)
    h = _step(z, step)
    numerator = _finite_complex(function(z + h*d), "function output")
    numerator -= _finite_complex(function(z - h*d), "function output")
    return _finite_complex(numerator / (2*h*d), "derivative")


def cauchy_riemann_residual(
    function: ComplexFunction,
    point: Complex,
    *,
    step: Real | None = None,
) -> tuple[float, float]:
    """Return approximate ``(u_x-v_y, u_y+v_x)`` for ``f=u+iv``.

    Near-zero residuals are necessary at differentiable holomorphic points,
    but not a sufficient holomorphicity certificate over a neighborhood.
    """
    z = _finite_complex(point, "point")
    h = _step(z, step)
    dx = (
        _finite_complex(function(z + h), "function output")
        - _finite_complex(function(z - h), "function output")
    ) / (2*h)
    dy = (
        _finite_complex(function(z + 1j*h), "function output")
        - _finite_complex(function(z - 1j*h), "function output")
    ) / (2*h)
    return (
        _finite_real(dx.real - dy.imag, "Cauchy–Riemann residual"),
        _finite_real(dy.real + dx.imag, "Cauchy–Riemann residual"),
    )


def _panels(value: int) -> int:
    if type(value) is not int or value < 2 or value % 2 or value > 1_000_000:
        raise ValueError("panels must be an even integer between 2 and 1,000,000")
    return value


def contour_integral(
    function: ComplexFunction,
    path: Path,
    path_derivative: Path,
    start: Real,
    stop: Real,
    *,
    panels: int = 512,
) -> complex:
    """Numerically evaluate ``∫ f(γ(t)) γ'(t) dt`` by composite Simpson.

    The caller supplies both the parametrized path and its derivative. Panel
    count is even. Discontinuities, cusps, poles and rapid oscillations require
    independent convergence checks, singularity analysis and/or split paths.
    """
    a, b = _finite_real(start, "start"), _finite_real(stop, "stop")
    n = _panels(panels)
    if a == b:
        return 0j
    h = (b-a) / n

    def integrand(t: float) -> complex:
        position = _finite_complex(path(t), "path output")
        derivative = _finite_complex(path_derivative(t), "path derivative")
        try:
            value = function(position)
        except ZeroDivisionError as exc:
            raise ValueError("integrand singular on sampled path") from exc
        result = _finite_complex(value, "integrand") * derivative
        return _finite_complex(result, "integrand times path derivative")

    total = integrand(a) + integrand(b)
    for i in range(1, n):
        total += (4 if i % 2 else 2) * integrand(a + i*h)
    return _finite_complex(total * h / 3, "contour integral")


def winding_number(
    path: Path,
    path_derivative: Path,
    center: Complex = 0j,
    *,
    start: Real = 0.0,
    stop: Real = 2*math.pi,
    panels: int = 512,
    tolerance: Real = 1e-3,
) -> int:
    """Estimate integer winding via ``(1/2πi) ∮ dz/(z-center)``.

    Require the sampled quadrature estimate to be near an integer, and the
    endpoints to coincide numerically. This does not certify a closed path is
    pole-free between sample locations. Use analytic/topological proof for
    rigorous winding computations, especially near the center.
    """
    a, b = _finite_real(start, "start"), _finite_real(stop, "stop")
    if a == b:
        raise ValueError("winding path must have a nonempty parameter interval")
    n = _panels(panels)
    tol = _finite_real(tolerance, "tolerance")
    if tol <= 0 or tol >= 0.25:
        raise ValueError("tolerance must be between 0 and 0.25")
    c = _finite_complex(center, "center")
    first = _finite_complex(path(a), "path start")
    last = _finite_complex(path(b), "path end")
    scale = max(1.0, abs(first), abs(last), abs(c))
    if abs(first-last) > min(tol, 1e-6) * scale:
        raise ValueError("winding number requires a closed path")
    for i in range(n+1):
        if abs(_finite_complex(path(a + (b-a)*i/n), "path output")-c) <= 1e-12 * scale:
            raise ValueError("path meets the winding center at a sampled point")
    integral = contour_integral(
        lambda z: 1/(z-c), path, path_derivative, a, b, panels=n
    )
    estimate = integral / (2j*math.pi)
    nearest = round(estimate.real)
    if abs(estimate-nearest) > tol:
        raise ValueError("estimated winding is not close to an integer; refine or inspect path")
    return nearest
