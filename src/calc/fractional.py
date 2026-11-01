"""Small, explicit fractional-derivative examples (not a general solver).

The power-function formulas are closed-form left-terminal-zero Caputo and
Riemann–Liouville derivatives of order 0 < alpha <= 1. The sampled finite
Grünwald–Letnikov convolution is a *Riemann–Liouville-type* approximation.
These definitions are not interchangeable for functions with f(0) != 0.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from numbers import Real


def _real(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number, not bool")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is outside finite floating-point range") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _power_parameters(alpha: Real, beta: Real, x: Real) -> tuple[float, float, float]:
    order = _real(alpha, "alpha")
    power = _real(beta, "beta")
    point = _real(x, "x")
    if not 0 < order <= 1 or power < 0 or point < 0:
        raise ValueError("requires 0 < alpha <= 1, beta >= 0, x >= 0")
    return order, power, point


def _power_result(order: float, power: float, point: float) -> float:
    """Evaluate Gamma(power+1)/Gamma(power-order+1) * x**(power-order)."""
    if order == 1 and point == 1:
        return power
    exponent = power - order
    try:
        if order == 1:
            # Γ(beta+1)/Γ(beta) = beta exactly. Differencing two enormous
            # log-gamma values needlessly loses many significant digits.
            log_coefficient = math.log(power)
        elif power >= 1_000_000:
            # Γ(z+a)/Γ(z) ≈ z**a exp(a(a-1)/(2z)), with z=beta+1-alpha.
            # The neglected relative term is O(z^-2) for fixed 0<a<1.
            z = power + 1 - order
            log_coefficient = order * math.log(z) + order * (order - 1) / (2 * z)
        else:
            log_coefficient = math.lgamma(power + 1) - math.lgamma(exponent + 1)
    except (ValueError, OverflowError) as exc:
        raise ValueError("fractional power coefficient cannot be evaluated") from exc
    if not math.isfinite(log_coefficient):
        raise ValueError("fractional power coefficient is nonfinite")
    if point == 0:
        if exponent > 0:
            return 0.0
        if exponent < 0:
            raise ValueError("fractional derivative diverges at x=0")
        log_result = log_coefficient
    else:
        log_result = log_coefficient + exponent * math.log(point)
    if not math.isfinite(log_result):
        raise ValueError("fractional derivative magnitude is nonfinite")
    try:
        result = math.exp(log_result)
    except OverflowError as exc:
        raise ValueError("fractional derivative exceeds finite float range") from exc
    if not math.isfinite(result):
        raise ValueError("fractional derivative is nonfinite")
    return result


def caputo_power_derivative(alpha: Real, beta: Real, x: Real) -> float:
    """Caputo ``D_C^alpha(t^beta)`` at x, left terminal 0.

    For 0 < alpha < 1 and beta > 0 this is
    ``Gamma(beta+1)/Gamma(beta-alpha+1) * x**(beta-alpha)``; the constant
    beta=0 has derivative 0. ``alpha=1`` is the ordinary first derivative.
    Finite x=0 limits are returned, and singular endpoints are rejected.
    """
    order, power, point = _power_parameters(alpha, beta, x)
    if power == 0:
        return 0.0
    return _power_result(order, power, point)


def riemann_liouville_power_derivative(alpha: Real, beta: Real, x: Real) -> float:
    """Left Riemann–Liouville ``D_RL^alpha(t^beta)`` at x, terminal 0.

    In contrast to Caputo, the fractional derivative of 1 is
    ``x**(-alpha)/Gamma(1-alpha)`` when 0 < alpha < 1 and x > 0.
    ``alpha=1`` reduces to the ordinary first derivative.
    """
    order, power, point = _power_parameters(alpha, beta, x)
    if order == 1 and power == 0:
        return 0.0
    return _power_result(order, power, point)


def grunwald_letnikov_derivative(
    function: Callable[[float], Real], x: Real, alpha: Real, *, steps: int = 1024
) -> float:
    """Left-zero Grünwald–Letnikov difference approximation at x > 0.

    Uses a uniform grid ``h=x/steps`` and weights
    ``w_k=(-1)**k * binomial(alpha,k)``, with recurrence
    ``w_(k+1) = w_k * (k-alpha)/(k+1)``. Returns
    ``h**(-alpha) * sum(w_k * f(x-k*h))`` for 0 < alpha < 1.

    This samples the full history [0,x], and approximates the
    Riemann–Liouville derivative under regularity/convergence conditions.
    It is **not** a general Caputo algorithm, certified error estimator,
    or stable method for singular/discontinuous input functions.
    """
    if not callable(function):
        raise TypeError("function must be callable")
    point = _real(x, "x")
    order = _real(alpha, "alpha")
    if point <= 0 or not 0 < order < 1:
        raise ValueError("requires x > 0 and 0 < alpha < 1")
    if type(steps) is not int or not 1 <= steps <= 200_000:
        raise ValueError("steps must be an integer from 1 to 200000")
    h = point / steps
    if h == 0 or not math.isfinite(h) or point - h == point:
        raise ValueError("fractional grid step is not representable")

    weight = 1.0
    terms: list[float] = []
    for k in range(steps + 1):
        # This form avoids making the left endpoint slightly negative from
        # accumulated floating-point subtraction.
        sample = point * ((steps - k) / steps)
        value = _real(function(sample), "function output")
        term = weight * value
        if not math.isfinite(term):
            raise ValueError("fractional convolution term is nonfinite")
        terms.append(term)
        weight *= (k - order) / (k + 1)
    try:
        weighted_sum = math.fsum(terms)
    except OverflowError as exc:
        raise ValueError("fractional convolution overflows") from exc
    if weighted_sum == 0:
        return 0.0
    magnitude = math.log(abs(weighted_sum)) - order * math.log(h)
    try:
        result = math.copysign(math.exp(magnitude), weighted_sum)
    except OverflowError as exc:
        raise ValueError("fractional derivative exceeds finite float range") from exc
    if not math.isfinite(result):
        raise ValueError("fractional derivative is nonfinite")
    return result
