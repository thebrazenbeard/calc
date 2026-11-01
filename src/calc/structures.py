"""Finite algebraic structures with explicit, verifiable contracts.

The operations in this module are deliberately small and exact where possible.
They are examples of group and field axioms, Boolean algebra, and the min-plus
tropical semiring, not generic implementations of all abstract algebra.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


def _strict_integer(value: object, name: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer (not bool)")
    return value


class FiniteGroup:
    """Validated finite group with elements labelled by integers ``0..n-1``.

    ``table[a][b]`` is the product ``a*b``. Construction verifies closure,
    associativity, identity and all inverses directly in O(n³) time. This
    exhaustive validation is practical for small educational groups only.
    """

    def __init__(self, table: Sequence[Sequence[int]]) -> None:
        if not isinstance(table, (tuple, list)):
            raise TypeError("table must be a sequence of rows")
        n = len(table)
        if n == 0 or n > 128:
            raise ValueError("group order must be between 1 and 128")
        rows: list[tuple[int, ...]] = []
        for row in table:
            if not isinstance(row, (tuple, list)) or len(row) != n:
                raise ValueError("group table must be square")
            converted = tuple(_strict_integer(x, "table entry") for x in row)
            if any(x < 0 or x >= n for x in converted):
                raise ValueError("group multiplication must be closed")
            rows.append(converted)
        products = tuple(rows)
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    if products[products[a][b]][c] != products[a][products[b][c]]:
                        raise ValueError(f"multiplication not associative at ({a},{b},{c})")
        identities = [e for e in range(n)
                      if all(products[e][x] == x and products[x][e] == x for x in range(n))]
        if len(identities) != 1:
            raise ValueError("group must have a unique two-sided identity")
        identity = identities[0]
        inverses: list[int] = []
        for a in range(n):
            valid = [b for b in range(n)
                     if products[a][b] == identity and products[b][a] == identity]
            if len(valid) != 1:
                raise ValueError(f"element {a} lacks a unique inverse")
            inverses.append(valid[0])
        self._table = products
        self._identity = identity
        self._inverse = tuple(inverses)

    @property
    def order(self) -> int:
        return len(self._table)

    @property
    def identity(self) -> int:
        return self._identity

    @property
    def is_abelian(self) -> bool:
        return all(self._table[a][b] == self._table[b][a]
                   for a in range(self.order) for b in range(self.order))

    def _element(self, value: int) -> int:
        _strict_integer(value, "group element")
        if not 0 <= value < self.order:
            raise ValueError("element outside group")
        return value

    def mul(self, a: int, b: int) -> int:
        """Multiply two element labels."""
        return self._table[self._element(a)][self._element(b)]

    def inverse(self, a: int) -> int:
        """Return the two-sided inverse label."""
        return self._inverse[self._element(a)]

    def element_order(self, a: int) -> int:
        """Smallest positive power of ``a`` equal to the group identity."""
        a = self._element(a)
        product = self.identity
        for exponent in range(1, self.order + 1):
            product = self.mul(product, a)
            if product == self.identity:
                return exponent
        raise RuntimeError("validated group violated Lagrange's theorem")

    def cyclic_generators(self) -> tuple[int, ...]:
        """All elements that generate the complete finite group."""
        return tuple(a for a in range(self.order) if self.element_order(a) == self.order)


def _prime_below_2_pow_64(number: int) -> bool:
    """Deterministic Miller–Rabin for all 0 <= n < 2**64."""
    if number < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if number % p == 0:
            return number == p
    d, s = number - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Sufficient for deterministic primality testing on unsigned 64-bit n:
    # Jim Sinclair / deterministic bases in known computational number theory.
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        a = base % number
        if a == 0:
            continue
        x = pow(a, d, number)
        if x in (1, number - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % number
            if x == number - 1:
                break
        else:
            return False
    return True


@dataclass(frozen=True)
class PrimeField:
    """Arithmetic on the finite field F_p for prime p < 2**64.

    Arguments are arbitrary integers, normalized modulo p; floats and bools
    are not silently coerced. Field instances operate on integer residues.
    """

    p: int

    def __post_init__(self) -> None:
        _strict_integer(self.p, "modulus")
        if not 2 <= self.p < 2**64 or not _prime_below_2_pow_64(self.p):
            raise ValueError("modulus must be a prime between 2 and 2**64-1")

    def _residue(self, value: int) -> int:
        return _strict_integer(value, "field element") % self.p

    def add(self, a: int, b: int) -> int:
        return (self._residue(a) + self._residue(b)) % self.p

    def sub(self, a: int, b: int) -> int:
        return (self._residue(a) - self._residue(b)) % self.p

    def mul(self, a: int, b: int) -> int:
        return (self._residue(a) * self._residue(b)) % self.p

    def inv(self, a: int) -> int:
        residue = self._residue(a)
        if residue == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        return pow(residue, -1, self.p)

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))

    def pow(self, a: int, exponent: int) -> int:
        _strict_integer(exponent, "exponent")
        residue = self._residue(a)
        if exponent < 0:
            residue = self.inv(residue)
        return pow(residue, abs(exponent), self.p)


def _boolean(value: bool) -> bool:
    if type(value) is not bool:
        raise TypeError("Boolean algebra operands must be bool, not int")
    return value


def boolean_join(a: bool, b: bool) -> bool:
    """Boolean join (logical OR)."""
    left, right = _boolean(a), _boolean(b)
    return left or right


def boolean_meet(a: bool, b: bool) -> bool:
    """Boolean meet (logical AND)."""
    left, right = _boolean(a), _boolean(b)
    return left and right


def boolean_complement(a: bool) -> bool:
    return not _boolean(a)


def boolean_xor(a: bool, b: bool) -> bool:
    return _boolean(a) != _boolean(b)


def _tropical_number(value: int | float) -> int | float:
    if type(value) not in (int, float):
        raise TypeError("tropical operands must be int or float")
    if isinstance(value, float) and (math.isnan(value) or value == -math.inf):
        raise ValueError("NaN and negative infinity are outside R union {+infinity}")
    return value


def tropical_add(a: int | float, b: int | float) -> int | float:
    """Min-plus addition: a ⊕ b = min(a,b); zero is +infinity."""
    return min(_tropical_number(a), _tropical_number(b))


def tropical_multiply(a: int | float, b: int | float) -> int | float:
    """Min-plus multiplication: a ⊗ b = a+b; one is 0."""
    a, b = _tropical_number(a), _tropical_number(b)
    if a == math.inf or b == math.inf:
        return math.inf
    result = a + b
    # A finite mathematical sum might overflow IEEE float; that is not an
    # exact result and must not silently become the absorbing element.
    if isinstance(result, float) and not math.isfinite(result):
        raise OverflowError("floating tropical product overflowed")
    return result
