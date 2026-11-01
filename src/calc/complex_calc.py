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
        result = 1e-5 * max(1.0, abs(point))
        if not math.isfinite(result):
            raise ValueError("default step is not finite at this point")
        return result
    result = _finite_real(step, "step")
    if result <= 0:
        raise ValueError("step must be positive")
    return result


def _perturbed(point: complex, delta: complex) -> tuple[complex, complex]:
    plus, minus = point + delta, point - delta
    if plus == point or minus == point or plus == minus:
        raise ValueError("step is not representable at this point and direction")
    return plus, minus


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
    norm = abs(d)
    if not math.isfinite(norm):
        raise ValueError("direction magnitude is not representable")
    d /= norm
    h = _step(z, step)
    plus, minus = _perturbed(z, h*d)
    numerator = _finite_complex(function(plus), "function output")
    numerator -= _finite_complex(function(minus), "function output")
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
    xp, xm = _perturbed(z, h)
    yp, ym = _perturbed(z, 1j*h)
    dx = (
        _finite_complex(function(xp), "function output")
        - _finite_complex(function(xm), "function output")
    ) / (2*h)
    dy = (
        _finite_complex(function(yp), "function output")
        - _finite_complex(function(ym), "function output")
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
    interval = b - a
    if not math.isfinite(interval):
        raise ValueError("contour parameter interval exceeds finite range")
    h = interval / n
    if h == 0 or a + h == a or b - h == b:
        raise ValueError("contour parameter step is not representable")

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
    local_radius = max(1.0, abs(first-c), abs(last-c))
    # Absolute closure with a small allowance for representable coordinate
    # rounding, not a relative tolerance inflated by a distant origin.
    closure_tolerance = max(
        min(tol, 1e-6),
        4 * max(math.ulp(v) for v in
                (first.real, first.imag, last.real, last.imag)),
    )
    if abs(first-last) > closure_tolerance:
        raise ValueError("winding number requires a closed path")
    for i in range(n+1):
        sample = _finite_complex(path(a + (b-a)*i/n), "path output")
        roundoff = max(math.ulp(v) for v in
                       (sample.real, sample.imag, c.real, c.imag))
        if abs(sample-c) <= max(1e-12 * local_radius, 8 * roundoff):
            raise ValueError("path meets the winding center at a sampled point")
    integral = contour_integral(
        lambda z: 1/(z-c), path, path_derivative, a, b, panels=n
    )
    estimate = integral / (2j*math.pi)
    nearest = round(estimate.real)
    if abs(estimate-nearest) > tol:
        raise ValueError("estimated winding is not close to an integer; refine or inspect path")
    return nearest
