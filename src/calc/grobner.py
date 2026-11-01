"""Exact, bounded Gröbner bases for small polynomial ideals over Q.

The ordering uses x0 > x1 > ...; ``lex`` and graded lexicographic ``grlex``
are supported. No floating point or arbitrary expression evaluation occurs.
This is an educational Buchberger implementation, not a large-scale CAS.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction


Monomial = tuple[int, ...]
Exact = int | Fraction
_MAX_VARIABLES = 8
_MAX_TERMS = 4096
_MAX_EXPONENT = 512
_MAX_COEFFICIENT_BITS = 4096
_MAX_PRODUCT_PAIRS = 100_000


class ResourceLimitError(RuntimeError):
    """A declared small-problem complexity or arithmetic budget was exceeded."""


def _fraction(value: Exact) -> Fraction:
    if type(value) is int or isinstance(value, Fraction):
        result = Fraction(value)
        if (result.numerator.bit_length() > _MAX_COEFFICIENT_BITS or
                result.denominator.bit_length() > _MAX_COEFFICIENT_BITS):
            raise ResourceLimitError("coefficient exceeds exact arithmetic bit budget")
        return result
    raise TypeError("polynomial coefficients must be int or Fraction, not float/bool")


def _dimension(nvars: int) -> int:
    if type(nvars) is not int or not 1 <= nvars <= _MAX_VARIABLES:
        raise ValueError(f"nvars must be an integer from 1 to {_MAX_VARIABLES}")
    return nvars


def _monomial(exponents: Sequence[int], nvars: int) -> Monomial:
    if not isinstance(exponents, (tuple, list)) or len(exponents) != nvars:
        raise ValueError("monomial must have one exponent per variable")
    if any(type(power) is not int or power < 0 for power in exponents):
        raise ValueError("monomial exponents must be nonnegative integers")
    if any(power > _MAX_EXPONENT for power in exponents):
        raise ResourceLimitError("monomial exponent budget exceeded")
    return tuple(exponents)


@dataclass(frozen=True, init=False, slots=True)
class MultiPolynomial:
    """Immutable sparse polynomial in ``nvars`` indeterminates over Q.

    ``terms`` are sorted tuples of ``(exponent_vector, Fraction)`` with zero
    coefficients deleted. Examples: ``{(2,0): 1, (0,1): -1}`` means x²-y.
    """

    nvars: int
    terms: tuple[tuple[Monomial, Fraction], ...]

    def __init__(self, terms: Mapping[Monomial, Exact] | None = None, *, nvars: int):
        n = _dimension(nvars)
        entries: dict[Monomial, Fraction] = {}
        if terms is not None:
            if not isinstance(terms, Mapping):
                raise TypeError("terms must be a monomial-to-coefficient mapping")
            if len(terms) > _MAX_TERMS:
                raise ResourceLimitError("polynomial term budget exceeded")
            for raw_monomial, raw_coefficient in terms.items():
                monomial = _monomial(raw_monomial, n)
                coefficient = _fraction(raw_coefficient)
                if coefficient:
                    entries[monomial] = coefficient
        object.__setattr__(self, "nvars", n)
        object.__setattr__(self, "terms", tuple(sorted(entries.items())))

    @classmethod
    def zero(cls, nvars: int) -> MultiPolynomial:
        return cls(nvars=nvars)

    @classmethod
    def one(cls, nvars: int) -> MultiPolynomial:
        return cls({(0,) * _dimension(nvars): 1}, nvars=nvars)

    @classmethod
    def constant(cls, value: Exact, nvars: int) -> MultiPolynomial:
        n = _dimension(nvars)
        coefficient = _fraction(value)
        return cls({(0,) * n: coefficient} if coefficient else {}, nvars=n)

    def __bool__(self) -> bool:
        return bool(self.terms)

    def _coerce(self, value: object) -> MultiPolynomial | object:
        if isinstance(value, MultiPolynomial):
            if self.nvars != value.nvars:
                raise ValueError("polynomials have different numbers of variables")
            return value
        if type(value) is int or isinstance(value, Fraction):
            return MultiPolynomial.constant(value, self.nvars)
        return NotImplemented

    def __neg__(self) -> MultiPolynomial:
        return MultiPolynomial({m: -c for m, c in self.terms}, nvars=self.nvars)

    def __add__(self, other: object) -> MultiPolynomial:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        result = dict(self.terms)
        for monomial, coefficient in rhs.terms:
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
        return MultiPolynomial(result, nvars=self.nvars)

    def __radd__(self, other: object) -> MultiPolynomial:
        return self.__add__(other)

    def __sub__(self, other: object) -> MultiPolynomial:
        rhs = self._coerce(other)
        return NotImplemented if rhs is NotImplemented else self + (-rhs)

    def __rsub__(self, other: object) -> MultiPolynomial:
        lhs = self._coerce(other)
        return NotImplemented if lhs is NotImplemented else lhs - self

    def __mul__(self, other: object) -> MultiPolynomial:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        if len(self.terms) * len(rhs.terms) > _MAX_PRODUCT_PAIRS:
            raise ResourceLimitError("polynomial multiplication pair budget exceeded")
        result: dict[Monomial, Fraction] = {}
        for a, ca in self.terms:
            for b, cb in rhs.terms:
                monomial = _monomial(tuple(x+y for x, y in zip(a, b)), self.nvars)
                result[monomial] = result.get(monomial, Fraction(0)) + ca*cb
        return MultiPolynomial(result, nvars=self.nvars)

    def __rmul__(self, other: object) -> MultiPolynomial:
        return self.__mul__(other)

    def __pow__(self, exponent: int) -> MultiPolynomial:
        if type(exponent) is not int or exponent < 0:
            raise ValueError("polynomial exponent must be a nonnegative integer")
        if exponent > _MAX_EXPONENT:
            raise ResourceLimitError("power budget exceeded")
        result, base = MultiPolynomial.one(self.nvars), self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def evaluate(self, values: Sequence[Exact]) -> Fraction:
        if not isinstance(values, (tuple, list)) or len(values) != self.nvars:
            raise ValueError("evaluation needs exactly nvars coordinate values")
        coordinates = tuple(_fraction(value) for value in values)
        return sum((coefficient * _product(x**degree for x, degree in zip(coordinates, monomial))
                    for monomial, coefficient in self.terms), Fraction(0))


def _product(values: Iterable[Fraction]) -> Fraction:
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def variables(nvars: int) -> tuple[MultiPolynomial, ...]:
    """Return ``x0, x1, ...`` in the chosen number of variables."""
    n = _dimension(nvars)
    return tuple(MultiPolynomial({tuple(int(i == j) for i in range(n)): 1}, nvars=n)
                 for j in range(n))


def _order(order: str) -> str:
    if order not in ("lex", "grlex"):
        raise ValueError("monomial order must be 'lex' or 'grlex'")
    return order


def _key(monomial: Monomial, order: str) -> tuple[int, ...]:
    return monomial if order == "lex" else (sum(monomial), *monomial)


def leading_term(poly: MultiPolynomial, order: str = "lex") -> tuple[Monomial, Fraction]:
    """Return leading exponent vector and coefficient; zero has none."""
    order = _order(order)
    if not poly:
        raise ValueError("the zero polynomial has no leading term")
    return max(poly.terms, key=lambda entry: _key(entry[0], order))


def _divides(divisor: Monomial, dividend: Monomial) -> bool:
    return all(a <= b for a, b in zip(divisor, dividend))


def _term_multiple(poly: MultiPolynomial, monomial: Monomial, coefficient: Fraction) -> MultiPolynomial:
    return MultiPolynomial({tuple(a+b for a, b in zip(m, monomial)): c*coefficient
                            for m, c in poly.terms}, nvars=poly.nvars)


def _budget(value: int, name: str) -> int:
    if type(value) is not int or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def normal_form(
    poly: MultiPolynomial, divisors: Sequence[MultiPolynomial],
    order: str = "lex", *, max_steps: int = 10_000, max_terms: int = _MAX_TERMS,
) -> MultiPolynomial:
    """Multivariate long-division remainder in the given divisor order.

    With a Gröbner basis the remainder is unique for the chosen monomial
    order. For arbitrary divisors it generally depends on their order.
    """
    order = _order(order)
    _budget(max_steps, "max_steps")
    _budget(max_terms, "max_terms")
    if not isinstance(poly, MultiPolynomial):
        raise TypeError("poly must be a MultiPolynomial")
    prepared = tuple(divisors)
    for divisor in prepared:
        if not isinstance(divisor, MultiPolynomial) or divisor.nvars != poly.nvars:
            raise ValueError("divisors must be polynomials in the same variables")
        if not divisor:
            raise ValueError("zero cannot be a polynomial divisor")
    if len(poly.terms) > max_terms or any(len(g.terms) > max_terms for g in prepared):
        raise ResourceLimitError("division input term budget exceeded")
    remainder = MultiPolynomial.zero(poly.nvars)
    current = poly
    for _ in range(max_steps):
        if not current:
            return remainder
        monomial, coefficient = leading_term(current, order)
        for divisor in prepared:
            dm, dc = leading_term(divisor, order)
            if _divides(dm, monomial):
                shift = tuple(a-b for a, b in zip(monomial, dm))
                current -= _term_multiple(divisor, shift, coefficient/dc)
                break
        else:
            head = MultiPolynomial({monomial: coefficient}, nvars=poly.nvars)
            remainder += head
            current -= head
        if len(current.terms) > max_terms or len(remainder.terms) > max_terms:
            raise ResourceLimitError("division term budget exceeded")
    if not current:
        return remainder
    raise ResourceLimitError("polynomial division step budget exceeded")


def s_polynomial(first: MultiPolynomial, second: MultiPolynomial, order: str = "lex") -> MultiPolynomial:
    """Buchberger S-polynomial cancelling the leading monomials."""
    order = _order(order)
    if first.nvars != second.nvars or not first or not second:
        raise ValueError("S-polynomial needs nonzero polynomials in the same variables")
    fm, fc = leading_term(first, order)
    gm, gc = leading_term(second, order)
    lcm = tuple(max(a, b) for a, b in zip(fm, gm))
    return (_term_multiple(first, tuple(a-b for a, b in zip(lcm, fm)), 1/fc)
            - _term_multiple(second, tuple(a-b for a, b in zip(lcm, gm)), 1/gc))


def _monic(poly: MultiPolynomial, order: str) -> MultiPolynomial:
    return _term_multiple(poly, (0,) * poly.nvars, 1/leading_term(poly, order)[1])


def grobner_basis(
    generators: Sequence[MultiPolynomial], order: str = "lex", *,
    max_pairs: int = 1000, max_steps: int = 10_000, max_terms: int = _MAX_TERMS,
    max_basis: int = 64,
) -> tuple[MultiPolynomial, ...]:
    """Return a monic, reduced Gröbner basis for a small ideal in Q[x].

    Buchberger's algorithm terminates for polynomial ideals over Q, but its
    worst-case resources can be enormous. Explicit budgets cause a visible
    ``ResourceLimitError``, never an apparently successful partial basis.
    """
    order = _order(order)
    for name, value in (("max_pairs", max_pairs), ("max_steps", max_steps),
                        ("max_terms", max_terms), ("max_basis", max_basis)):
        _budget(value, name)
    raw_items = tuple(generators)
    if not raw_items:
        raise ValueError("at least one generator is required")
    examples = tuple(g for g in raw_items if isinstance(g, MultiPolynomial))
    if not examples:
        raise TypeError("at least one generator must establish the variable dimension")
    n = examples[0].nvars
    items = tuple(MultiPolynomial.constant(g, n) if type(g) is int or isinstance(g, Fraction)
                  else g for g in raw_items)
    if not all(isinstance(g, MultiPolynomial) for g in items):
        raise TypeError("generators must be exact polynomials or scalar constants")
    if any(g.nvars != n for g in items):
        raise ValueError("all generators must use the same number of variables")
    basis: list[MultiPolynomial] = []
    for generator in items:
        if generator:
            monic = _monic(generator, order)
            if monic not in basis:
                basis.append(monic)
    if not basis:
        raise ValueError("zero generators do not define a nontrivial Gröbner task")
    one = MultiPolynomial.one(n)
    if any(g == one for g in basis):
        return (one,)
    if len(basis) > max_basis:
        raise ResourceLimitError("basis size budget exceeded")
    pending = deque((i, j) for i in range(len(basis)) for j in range(i+1, len(basis)))
    pairs_processed = 0
    while pending:
        if pairs_processed >= max_pairs:
            raise ResourceLimitError("Buchberger pair budget exceeded")
        i, j = pending.popleft()
        pairs_processed += 1
        remainder = normal_form(s_polynomial(basis[i], basis[j], order), basis,
                                order, max_steps=max_steps, max_terms=max_terms)
        if remainder:
            new = _monic(remainder, order)
            if new == one:
                return (one,)
            if len(basis) >= max_basis:
                raise ResourceLimitError("basis size budget exceeded")
            index = len(basis)
            pending.extend((k, index) for k in range(index))
            basis.append(new)

    # Minimal leading-monomial basis: discard terms with divisible leaders.
    minimal: list[MultiPolynomial] = []
    for i, f in enumerate(basis):
        fm = leading_term(f, order)[0]
        if any(i != j and _divides(leading_term(g, order)[0], fm)
               and (leading_term(g, order)[0] != fm or j < i)
               for j, g in enumerate(basis)):
            continue
        minimal.append(f)

    # Interreduce tails against the other (minimal) basis elements.
    reduced: list[MultiPolynomial] = []
    for i, f in enumerate(minimal):
        others = minimal[:i] + minimal[i+1:]
        remainder = normal_form(f, others, order,
                                max_steps=max_steps, max_terms=max_terms)
        if remainder:
            reduced.append(_monic(remainder, order))
    return tuple(sorted(set(reduced),
                        key=lambda g: _key(leading_term(g, order)[0], order),
                        reverse=True))
