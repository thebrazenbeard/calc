"""Exact-rational piecewise-linear finite elements for a bounded 1D Poisson problem.

This module solves ``-u'' = f`` on a finite interval with homogeneous Dirichlet
conditions, using a nonuniform rational mesh and a *piecewise constant* load.
Each element is integrated exactly. It is intentionally not a general PDE/FEM
framework, arbitrary forcing-function integrator, or error-certification tool.
"""

from __future__ import annotations

from bisect import bisect_right
from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice

Rational = int | Fraction

_MAX_NODES = 513
_MAX_BITS = 4096


def _rational(value: object, label: str) -> Fraction:
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError(f"{label} must be an exact int or Fraction")
    result = Fraction(value)
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > _MAX_BITS:
        raise ValueError(f"{label} exceeds the rational bit limit")
    return result


def _checked(value: Fraction, label: str) -> Fraction:
    if max(value.numerator.bit_length(), value.denominator.bit_length()) > _MAX_BITS:
        raise ValueError(f"{label} exceeds the rational bit limit")
    return value


def _nodes(nodes: Iterable[Rational]) -> tuple[Fraction, ...]:
    if isinstance(nodes, (str, bytes)):
        raise TypeError("nodes must be a finite iterable of rational coordinates")
    try:
        raw = tuple(islice(iter(nodes), _MAX_NODES + 1))
    except TypeError as exc:
        raise TypeError("nodes must be an iterable of rational coordinates") from exc
    if not 2 <= len(raw) <= _MAX_NODES:
        raise ValueError("mesh must contain between 2 and 513 nodes")
    values = tuple(_rational(v, "node") for v in raw)
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("mesh coordinates must be strictly increasing")
    return values


def _forcing(load: Rational | Iterable[Rational], elements: int) -> tuple[Fraction, ...]:
    if type(load) is int or isinstance(load, Fraction):
        return (_rational(load, "forcing"),) * elements
    if isinstance(load, (str, bytes)):
        raise TypeError("forcing must be rational or an iterable of rationals")
    try:
        raw = tuple(islice(iter(load), elements + 1))
    except TypeError as exc:
        raise TypeError("forcing must be rational or an iterable of rationals") from exc
    if len(raw) != elements:
        raise ValueError("piecewise-constant forcing needs one value per element")
    return tuple(_rational(v, "forcing") for v in raw)


@dataclass(frozen=True, slots=True)
class Poisson1DSystem:
    """Interior-node Galerkin system ``K u = b`` with tridiagonal stiffness."""

    nodes: tuple[Fraction, ...]
    forcing: tuple[Fraction, ...]
    diagonal: tuple[Fraction, ...]
    upper: tuple[Fraction, ...]
    load: tuple[Fraction, ...]

    @property
    def degrees_of_freedom(self) -> int:
        return len(self.diagonal)

    def residual(self, interior_values: Iterable[Rational]) -> tuple[Fraction, ...]:
        """Return the exact residual ``K v - b`` at interior vertices."""
        try:
            raw = tuple(islice(iter(interior_values), self.degrees_of_freedom + 1))
        except TypeError as exc:
            raise TypeError("interior values must be an iterable") from exc
        if len(raw) != self.degrees_of_freedom:
            raise ValueError("one interior value is required per degree of freedom")
        v = tuple(_rational(value, "interior value") for value in raw)
        residual = []
        for i, d in enumerate(self.diagonal):
            value = d * v[i] - self.load[i]
            if i:
                value += self.upper[i - 1] * v[i - 1]
            if i + 1 < len(v):
                value += self.upper[i] * v[i + 1]
            residual.append(_checked(value, "residual"))
        return tuple(residual)

    def solve(self) -> Poisson1DSolution:
        """Solve the positive-definite tridiagonal system in O(n) rational ops."""
        n = self.degrees_of_freedom
        if not n:
            return Poisson1DSolution(self, (Fraction(0), Fraction(0)))
        d = list(self.diagonal)
        rhs = list(self.load)
        for i in range(1, n):
            if d[i - 1] <= 0:
                raise ArithmeticError("nonpositive FEM pivot")
            factor = _checked(self.upper[i - 1] / d[i - 1], "elimination factor")
            d[i] = _checked(d[i] - factor * self.upper[i - 1], "FEM pivot")
            rhs[i] = _checked(rhs[i] - factor * rhs[i - 1], "FEM right-hand side")
        if d[-1] <= 0:
            raise ArithmeticError("nonpositive FEM pivot")
        v = [Fraction(0)] * n
        v[-1] = _checked(rhs[-1] / d[-1], "FEM solution")
        for i in range(n - 2, -1, -1):
            v[i] = _checked((rhs[i] - self.upper[i] * v[i + 1]) / d[i], "FEM solution")
        return Poisson1DSolution(self, (Fraction(0), *v, Fraction(0)))


@dataclass(frozen=True, slots=True)
class Poisson1DSolution:
    """Continuous, piecewise-affine interpolant of the FEM nodal solution."""

    system: Poisson1DSystem
    nodal_values: tuple[Fraction, ...]

    @property
    def nodes(self) -> tuple[Fraction, ...]:
        return self.system.nodes

    def evaluate(self, point: Rational) -> Fraction:
        """Evaluate the exact affine interpolant inside the mesh interval."""
        x = _rational(point, "point")
        if not self.nodes[0] <= x <= self.nodes[-1]:
            raise ValueError("point must lie inside the mesh interval")
        i = min(bisect_right(self.nodes, x) - 1, len(self.nodes) - 2)
        h = self.nodes[i + 1] - self.nodes[i]
        t = (x - self.nodes[i]) / h
        return _checked((1 - t) * self.nodal_values[i] + t * self.nodal_values[i + 1],
                        "FEM evaluation")

    def energy_seminorm_squared(self) -> Fraction:
        """Return the exact integral of the squared piecewise constant slope."""
        total = Fraction(0)
        for x0, x1, u0, u1 in zip(self.nodes, self.nodes[1:],
                                  self.nodal_values, self.nodal_values[1:]):
            total = _checked(total + (u1 - u0) ** 2 / (x1 - x0), "FEM energy")
        return total


def assemble_poisson_1d(
    nodes: Iterable[Rational], forcing: Rational | Iterable[Rational] = 1,
) -> Poisson1DSystem:
    """Assemble P1 weak form of ``-u''=f``, ``u(a)=u(b)=0``.

    For element length ``h`` and constant forcing ``f_e`` the local stiffness
    is ``[[1,-1],[-1,1]]/h`` and local load is ``f_e*h*(1,1)/2``.
    Only interior vertices are unknown, giving an SPD tridiagonal system.
    """
    x = _nodes(nodes)
    f = _forcing(forcing, len(x) - 1)
    n = len(x) - 2
    h = tuple(_checked(b - a, "element length") for a, b in zip(x, x[1:]))
    diag = tuple(_checked(1 / h[i] + 1 / h[i + 1], "stiffness") for i in range(n))
    upper = tuple(_checked(-1 / h[i + 1], "stiffness") for i in range(n - 1))
    load = tuple(_checked((f[i] * h[i] + f[i + 1] * h[i + 1]) / 2,
                          "element load") for i in range(n))
    return Poisson1DSystem(x, f, diag, upper, load)


def solve_poisson_1d(
    nodes: Iterable[Rational], forcing: Rational | Iterable[Rational] = 1,
) -> Poisson1DSolution:
    """Assemble and solve the finite-element 1D Poisson problem."""
    return assemble_poisson_1d(nodes, forcing).solve()
