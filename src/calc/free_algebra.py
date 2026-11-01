"""Small exact free associative (noncommutative) polynomial algebras over Q.

A monomial is a finite ordered *word* such as ``('x', 'y', 'x')``.
The empty word is the unit. Multiplication concatenates words: in general
``x*y != y*x``. This module supports a bounded algebraic demonstration,
not completed formal power series or arbitrary noncommutative Groebner bases.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import re


Exact = int | Fraction
Word = tuple[str, ...]
_MAX_GENERATORS = 8
_MAX_SYMBOL_LENGTH = 32
_MAX_WORD_LENGTH = 64
_MAX_TERMS = 4096
_MAX_PRODUCT_PAIRS = 100_000
_MAX_COEFFICIENT_BITS = 4096
_MAX_POWER = 32
_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")


class ResourceLimitError(RuntimeError):
    """A declared bounded-computation budget was exceeded."""


def _rational(value: Exact) -> Fraction:
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError("coefficients must be exact int or Fraction, never bool/float")
    rational = Fraction(value)
    if (rational.numerator.bit_length() > _MAX_COEFFICIENT_BITS or
            rational.denominator.bit_length() > _MAX_COEFFICIENT_BITS):
        raise ResourceLimitError("rational coefficient bit budget exceeded")
    return rational


@dataclass(frozen=True, slots=True)
class FreeAlgebra:
    """A finitely generated free associative algebra over rational numbers."""

    generators: tuple[str, ...]

    def __post_init__(self) -> None:
        names = self.generators
        if (not isinstance(names, tuple) or not 1 <= len(names) <= _MAX_GENERATORS):
            raise ValueError(f"generators must be a tuple of 1 to {_MAX_GENERATORS} names")
        if any(not isinstance(name, str) or len(name) > _MAX_SYMBOL_LENGTH
               or _NAME.fullmatch(name) is None for name in names):
            raise ValueError("generator names must be short ASCII identifiers")
        if len(set(names)) != len(names):
            raise ValueError("generator names must be distinct")

    def zero(self) -> FreePolynomial:
        return FreePolynomial(self)

    def one(self) -> FreePolynomial:
        return self.constant(1)

    def constant(self, coefficient: Exact) -> FreePolynomial:
        value = _rational(coefficient)
        return FreePolynomial(self, {(): value} if value else {})

    def monomial(self, word: Sequence[str], coefficient: Exact = 1) -> FreePolynomial:
        if not isinstance(word, (tuple, list)):
            raise TypeError("word must be a tuple or list of generator names")
        return FreePolynomial(self, {tuple(word): coefficient})

    def variable(self, name: str) -> FreePolynomial:
        if name not in self.generators:
            raise ValueError("unknown generator name")
        return self.monomial((name,))


@dataclass(frozen=True, init=False, slots=True)
class FreePolynomial:
    """Immutable sparse exact polynomial in ordered noncommuting words."""

    algebra: FreeAlgebra
    terms: tuple[tuple[Word, Fraction], ...]

    def __init__(self, algebra: FreeAlgebra,
                 terms: Mapping[Word, Exact] | None = None) -> None:
        if not isinstance(algebra, FreeAlgebra):
            raise TypeError("algebra must be a FreeAlgebra")
        if terms is None:
            terms = {}
        if not isinstance(terms, Mapping):
            raise TypeError("terms must be a mapping of words to coefficients")
        if len(terms) > _MAX_TERMS:
            raise ResourceLimitError("polynomial term budget exceeded")
        normalized: dict[Word, Fraction] = {}
        for word, raw in terms.items():
            if not isinstance(word, tuple):
                raise TypeError("monomial keys must be tuples of generator names")
            if len(word) > _MAX_WORD_LENGTH:
                raise ResourceLimitError("monomial word length budget exceeded")
            if any(not isinstance(symbol, str) or symbol not in algebra.generators
                   for symbol in word):
                raise ValueError("word contains an unknown generator")
            value = _rational(raw)
            if value:
                normalized[word] = value
        object.__setattr__(self, "algebra", algebra)
        object.__setattr__(self, "terms", tuple(sorted(normalized.items())))

    def as_dict(self) -> dict[Word, Fraction]:
        return dict(self.terms)

    def __bool__(self) -> bool:
        return bool(self.terms)

    def _coerce(self, other: object) -> FreePolynomial | object:
        if isinstance(other, FreePolynomial):
            if self.algebra != other.algebra:
                raise ValueError("polynomials belong to different free algebras")
            return other
        if type(other) is int or isinstance(other, Fraction):
            return self.algebra.constant(other)
        return NotImplemented

    def __neg__(self) -> FreePolynomial:
        return FreePolynomial(self.algebra, {word: -coefficient
                                             for word, coefficient in self.terms})

    def __add__(self, other: object) -> FreePolynomial:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        values = self.as_dict()
        for word, coefficient in rhs.terms:
            values[word] = _rational(values.get(word, Fraction(0)) + coefficient)
        return FreePolynomial(self.algebra, values)

    __radd__ = __add__

    def __sub__(self, other: object) -> FreePolynomial:
        rhs = self._coerce(other)
        return NotImplemented if rhs is NotImplemented else self + (-rhs)

    def __rsub__(self, other: object) -> FreePolynomial:
        lhs = self._coerce(other)
        return NotImplemented if lhs is NotImplemented else lhs - self

    def __mul__(self, other: object) -> FreePolynomial:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        if len(self.terms) * len(rhs.terms) > _MAX_PRODUCT_PAIRS:
            raise ResourceLimitError("noncommutative product pair budget exceeded")
        values: dict[Word, Fraction] = {}
        for lhs_word, lhs_coefficient in self.terms:
            for rhs_word, rhs_coefficient in rhs.terms:
                word = lhs_word + rhs_word
                if len(word) > _MAX_WORD_LENGTH:
                    raise ResourceLimitError("product word length budget exceeded")
                values[word] = _rational(values.get(word, Fraction(0)) +
                                         lhs_coefficient * rhs_coefficient)
                if len(values) > _MAX_TERMS:
                    raise ResourceLimitError("product term budget exceeded")
        return FreePolynomial(self.algebra, values)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> FreePolynomial:
        if type(exponent) is not int or not 0 <= exponent <= _MAX_POWER:
            raise ValueError(f"exponent must be an integer from 0 to {_MAX_POWER}")
        result = self.algebra.one()
        base = self
        n = exponent
        while n:
            if n & 1:
                result = result * base
            n >>= 1
            if n:
                base = base * base
        return result

    def commutator(self, other: FreePolynomial) -> FreePolynomial:
        """The Lie bracket ``[P, Q] = P Q - Q P`` in the free algebra."""
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            raise TypeError("commutator requires a polynomial or exact scalar")
        return self * rhs - rhs * self

    def derivative(self, generator: str) -> FreePolynomial:
        """Free derivation D with D(generator)=1 and D(other generator)=0.

        The Leibniz extension deletes *each occurrence* of the selected
        generator while preserving the order of all remaining symbols.
        """
        if generator not in self.algebra.generators:
            raise ValueError("unknown differentiation generator")
        values: dict[Word, Fraction] = {}
        for word, coefficient in self.terms:
            for index, symbol in enumerate(word):
                if symbol == generator:
                    derivative_word = word[:index] + word[index + 1:]
                    values[derivative_word] = _rational(
                        values.get(derivative_word, Fraction(0)) + coefficient)
                    if len(values) > _MAX_TERMS:
                        raise ResourceLimitError("derivation term budget exceeded")
        return FreePolynomial(self.algebra, values)

    def derivation(self, images: Mapping[str, FreePolynomial]) -> FreePolynomial:
        """Unique derivation D extending the supplied images of generators.

        Omitted generators map to zero. For each word ``w_1...w_k``, computes
        ``sum_i w_1...w_(i-1) D(w_i) w_(i+1)...w_k`` in the original order.
        """
        if not isinstance(images, Mapping):
            raise TypeError("images must map generator names to polynomials")
        normalized: dict[str, FreePolynomial] = {}
        for symbol, image in images.items():
            if symbol not in self.algebra.generators:
                raise ValueError("unknown derivation generator")
            if not isinstance(image, FreePolynomial):
                raise TypeError("each derivation image must be a FreePolynomial")
            if image.algebra != self.algebra:
                raise ValueError("derivation image belongs to a different free algebra")
            normalized[symbol] = image
        values: dict[Word, Fraction] = {}
        operations = 0
        for word, coefficient in self.terms:
            for index, symbol in enumerate(word):
                image = normalized.get(symbol)
                if not image:
                    continue
                operations += len(image.terms)
                if operations > _MAX_PRODUCT_PAIRS:
                    raise ResourceLimitError("derivation pair budget exceeded")
                for inserted, inserted_coefficient in image.terms:
                    derived_word = word[:index] + inserted + word[index + 1:]
                    if len(derived_word) > _MAX_WORD_LENGTH:
                        raise ResourceLimitError("derived word length budget exceeded")
                    values[derived_word] = _rational(
                        values.get(derived_word, Fraction(0)) + coefficient * inserted_coefficient)
                    if len(values) > _MAX_TERMS:
                        raise ResourceLimitError("derived term budget exceeded")
        return FreePolynomial(self.algebra, values)


__all__ = ["FreeAlgebra", "FreePolynomial", "ResourceLimitError"]
