"""Exact elementary algebra over the rational numbers.

No optional dependencies, expression evaluation, or implicit floating-point
rounding. Polynomial coefficients run in ascending powers: ``(3, 2, 1)``
represents ``3 + 2*x + x**2``. Matrix inputs are rectangular nested sequences.
All public results are immutable and exact; integer entries become ``Fraction``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence


Exact = int | Fraction
Matrix = tuple[tuple[Fraction, ...], ...]


def _fraction(value: Exact) -> Fraction:
    """Convert only integers and rational values (never implicit floats)."""
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    raise TypeError("exact algebra accepts integers and Fraction, not floats")


@dataclass(frozen=True, init=False)
class Polynomial:
    """Immutable polynomial with low-to-high rational coefficients.

    The zero polynomial is represented canonically by ``(Fraction(0),)``
    and has degree ``-1``. Derivation and integration are exact.
    """

    coefficients: tuple[Fraction, ...]

    def __init__(self, coefficients: Iterable[Exact] = ()) -> None:
        normalized = [_fraction(value) for value in coefficients]
        while len(normalized) > 1 and not normalized[-1]:
            normalized.pop()
        if not normalized:
            normalized = [Fraction(0)]
        object.__setattr__(self, "coefficients", tuple(normalized))

    @property
    def degree(self) -> int:
        """Return the highest nonzero exponent, or -1 for zero."""
        return -1 if self.coefficients == (Fraction(0),) else len(self.coefficients) - 1

    def __neg__(self) -> Polynomial:
        return Polynomial(-c for c in self.coefficients)

    def __add__(self, other: Polynomial) -> Polynomial:
        if not isinstance(other, Polynomial):
            return NotImplemented
        return Polynomial(
            (self.coefficients[i] if i < len(self.coefficients) else 0)
            + (other.coefficients[i] if i < len(other.coefficients) else 0)
            for i in range(max(len(self.coefficients), len(other.coefficients)))
        )

    def __sub__(self, other: Polynomial) -> Polynomial:
        if not isinstance(other, Polynomial):
            return NotImplemented
        return self + (-other)

    def __mul__(self, other: Polynomial) -> Polynomial:
        if not isinstance(other, Polynomial):
            return NotImplemented
        result = [Fraction(0) for _ in range(len(self.coefficients) + len(other.coefficients) - 1)]
        for i, a in enumerate(self.coefficients):
            for j, b in enumerate(other.coefficients):
                result[i + j] += a * b
        return Polynomial(result)

    def __pow__(self, exponent: int) -> Polynomial:
        """Return a polynomial raised to a nonnegative integer power."""
        if type(exponent) is not int or exponent < 0:
            raise ValueError("polynomial exponent must be a nonnegative integer")
        result = Polynomial((1,))
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent //= 2
        return result

    def __divmod__(self, other: Polynomial) -> tuple[Polynomial, Polynomial]:
        """Euclidean division over rational polynomial coefficients.

        Returns quotient and remainder with remainder degree smaller than the
        divisor's. The zero divisor raises ``ZeroDivisionError``.
        """
        if not isinstance(other, Polynomial):
            return NotImplemented
        if other.degree < 0:
            raise ZeroDivisionError("cannot divide by zero polynomial")
        quotient = [Fraction(0)] * max(1, self.degree - other.degree + 1)
        remainder = self
        while remainder.degree >= other.degree:
            offset = remainder.degree - other.degree
            scale = remainder.coefficients[-1] / other.coefficients[-1]
            quotient[offset] += scale
            remainder -= Polynomial([0] * offset + [scale * c for c in other.coefficients])
        return Polynomial(quotient), remainder

    def gcd(self, other: Polynomial) -> Polynomial:
        """Return the monic polynomial greatest common divisor over Q[x].

        The greatest common divisor of two zero polynomials is zero.
        """
        if not isinstance(other, Polynomial):
            raise TypeError("polynomial gcd requires another Polynomial")
        a, b = self, other
        while b.degree >= 0:
            _, remainder = divmod(a, b)
            a, b = b, remainder
        if a.degree < 0:
            return a
        leading = a.coefficients[-1]
        return Polynomial(c / leading for c in a.coefficients)

    def __call__(self, x: Exact) -> Fraction:
        """Evaluate by Horner's method without losing rational precision."""
        argument = _fraction(x)
        result = Fraction(0)
        for coefficient in reversed(self.coefficients):
            result = result * argument + coefficient
        return result

    def differentiate(self, order: int = 1) -> Polynomial:
        """Return the ``order``-th formal derivative (order >= 0)."""
        if type(order) is not int or order < 0:
            raise ValueError("derivative order must be a nonnegative integer")
        result = self
        for _ in range(order):
            result = Polynomial(i * c for i, c in enumerate(result.coefficients) if i)
        return result

    def integrate(self, constant: Exact = 0) -> Polynomial:
        """Return an antiderivative with the specified constant term."""
        return Polynomial([_fraction(constant)] + [
            c / (i + 1) for i, c in enumerate(self.coefficients)
        ])


def _matrix(values: Sequence[Sequence[Exact]]) -> Matrix:
    """Validate and normalize a nonempty rectangular rational matrix."""
    rows = tuple(tuple(_fraction(value) for value in row) for row in values)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("matrix must be nonempty and rectangular")
    return rows


def transpose(values: Sequence[Sequence[Exact]]) -> Matrix:
    """Return a matrix transposed without changing its input."""
    matrix = _matrix(values)
    return tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                 for j in range(len(matrix[0])))


def matmul(left: Sequence[Sequence[Exact]], right: Sequence[Sequence[Exact]]) -> Matrix:
    """Multiply conformable matrices exactly; reject dimension mismatches."""
    a, b = _matrix(left), _matrix(right)
    if len(a[0]) != len(b):
        raise ValueError("incompatible dimensions for matrix multiplication")
    return tuple(
        tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), Fraction(0))
              for j in range(len(b[0])))
        for i in range(len(a))
    )


def trace(values: Sequence[Sequence[Exact]]) -> Fraction:
    """Return the sum of diagonal entries of a square matrix."""
    matrix = _matrix(values)
    if len(matrix) != len(matrix[0]):
        raise ValueError("trace requires a square matrix")
    return sum((row[i] for i, row in enumerate(matrix)), Fraction(0))


def matrix_rank(values: Sequence[Sequence[Exact]]) -> int:
    """Compute exact row rank of a rectangular matrix by elimination."""
    matrix = _matrix(values)
    a = [list(row) for row in matrix]
    rows, cols = len(a), len(a[0])
    pivot_row = 0
    for col in range(cols):
        index = next((r for r in range(pivot_row, rows) if a[r][col]), None)
        if index is None:
            continue
        a[pivot_row], a[index] = a[index], a[pivot_row]
        pivot = a[pivot_row][col]
        for r in range(pivot_row + 1, rows):
            factor = a[r][col] / pivot
            for j in range(col + 1, cols):
                a[r][j] -= factor * a[pivot_row][j]
            a[r][col] = Fraction(0)
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(values: Sequence[Sequence[Exact]]) -> Fraction:
    """Compute a square matrix determinant by rational elimination.

    A singular matrix returns zero. No rounding tolerance is involved.
    """
    source = _matrix(values)
    n = len(source)
    if any(len(row) != n for row in source):
        raise ValueError("determinant requires a square matrix")
    a = [list(row) for row in source]
    sign = 1
    result = Fraction(1)
    for col in range(n):
        pivot_index = next((r for r in range(col, n) if a[r][col]), None)
        if pivot_index is None:
            return Fraction(0)
        if pivot_index != col:
            a[col], a[pivot_index] = a[pivot_index], a[col]
            sign *= -1
        pivot = a[col][col]
        result *= pivot
        for row in range(col + 1, n):
            factor = a[row][col] / pivot
            for j in range(col + 1, n):
                a[row][j] -= factor * a[col][j]
            a[row][col] = Fraction(0)
    return result * sign


def solve_linear_system(
    coefficients: Sequence[Sequence[Exact]],
    constants: Sequence[Exact],
) -> tuple[Fraction, ...]:
    """Solve ``A*x=b`` exactly when square A has a unique solution.

    Raises ``ValueError`` for dimension mismatch or a singular system.
    """
    matrix = _matrix(coefficients)
    n = len(matrix)
    if len(matrix[0]) != n or len(constants) != n:
        raise ValueError("linear system requires square A and matching b")
    a = [list(row) + [_fraction(b)] for row, b in zip(matrix, constants)]
    for col in range(n):
        pivot_index = next((row for row in range(col, n) if a[row][col]), None)
        if pivot_index is None:
            raise ValueError("linear system is singular")
        if pivot_index != col:
            a[col], a[pivot_index] = a[pivot_index], a[col]
        pivot = a[col][col]
        a[col] = [x / pivot for x in a[col]]
        for row in range(n):
            if row != col:
                factor = a[row][col]
                a[row] = [a[row][j] - factor * a[col][j] for j in range(n + 1)]
    return tuple(a[row][-1] for row in range(n))


def inverse_matrix(values: Sequence[Sequence[Exact]]) -> Matrix:
    """Return a square rational matrix inverse; reject singular matrices."""
    matrix = _matrix(values)
    n = len(matrix)
    if len(matrix[0]) != n:
        raise ValueError("matrix inverse requires a square matrix")
    columns = [solve_linear_system(matrix, [int(row == col) for row in range(n)])
               for col in range(n)]
    return tuple(tuple(columns[col][row] for col in range(n)) for row in range(n))


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return ``(g,x,y)`` with nonnegative gcd ``g=a*x+b*y``."""
    if type(a) is not int or type(b) is not int:
        raise TypeError("extended_gcd requires integers")
    old_r, r = abs(a), abs(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, old_s if a >= 0 else -old_s, old_t if b >= 0 else -old_t


def mod_inverse(a: int, modulus: int) -> int:
    """Return inverse of ``a`` modulo ``modulus > 1``, if one exists."""
    if type(modulus) is not int or modulus <= 1:
        raise ValueError("modulus must be an integer greater than 1")
    g, x, _ = extended_gcd(a, modulus)
    if g != 1:
        raise ValueError("modular inverse does not exist")
    return x % modulus


def chinese_remainder(congruences: Iterable[tuple[int, int]]) -> tuple[int, int]:
    """Generalized Chinese remainder theorem, even for noncoprime moduli.

    Returns ``(residue, least_positive_modulus)`` with normalized residue.
    Inconsistent congruences raise ``ValueError``. An empty system is ``(0,1)``.
    """
    residue, modulus = 0, 1
    for requested_residue, next_modulus in congruences:
        if type(requested_residue) is not int or type(next_modulus) is not int:
            raise TypeError("congruences must have integer residue and modulus")
        if next_modulus <= 0:
            raise ValueError("moduli must be positive integers")
        a = requested_residue % next_modulus
        g, _, _ = extended_gcd(modulus, next_modulus)
        difference = a - residue
        if difference % g:
            raise ValueError("congruences are inconsistent")
        reduced_modulus = next_modulus // g
        if reduced_modulus == 1:
            factor = 0
        else:
            factor = ((difference // g) * mod_inverse(modulus // g, reduced_modulus)) % reduced_modulus
        residue = (residue + modulus * factor) % (modulus * reduced_modulus)
        modulus *= reduced_modulus
    return residue, modulus
