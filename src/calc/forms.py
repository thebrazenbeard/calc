"""Exact exterior algebra and Euclidean Hodge star in a fixed oriented frame.

The coefficient field is the rational numbers.  A basis ``(i, j, ...)``
represents ``dx_i ∧ dx_j ∧ ...``, with zero-based, strictly increasing indices.
There is no coordinate dependence, manifold, exterior derivative, or connection
in this module: these are constant-coefficient forms on an oriented Euclidean
vector space with an orthonormal basis.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction


def _exact_scalar(value: int | Fraction) -> Fraction:
    """Convert an integer/rational, rejecting inexact floating-point inputs."""
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("coefficients must be int or Fraction (never float/bool)")
    return Fraction(value)


def _index_tuple(indices: tuple[int, ...], *, dimension: int, degree: int) -> None:
    if not isinstance(indices, tuple) or len(indices) != degree:
        raise ValueError("basis indices must be a tuple of the form's degree")
    previous = -1
    for index in indices:
        if isinstance(index, bool) or not isinstance(index, int):
            raise TypeError("basis indices must be integers (not booleans)")
        if index <= previous or index >= dimension:
            raise ValueError("basis indices must be strictly increasing and in range")
        previous = index


@dataclass(frozen=True, slots=True, init=False)
class Form:
    """An immutable sparse rational ``degree``-form in ``dimension`` dimensions.

    ``Form(3, 2, {(0, 1): Fraction(1, 2)})`` represents
    ``(1/2) dx_0 ∧ dx_1``.  Zero terms are removed and coefficients are stored
    in a canonical, sorted immutable tuple.  The empty basis ``()`` denotes
    the unit of the exterior algebra (a 0-form).

    Degrees higher than the ambient dimension are permitted *only* for a zero
    form, which is useful because the wedge product has a graded result even
    when its value vanishes identically.
    """

    dimension: int
    degree: int
    coefficients: tuple[tuple[tuple[int, ...], Fraction], ...] = field(repr=True)

    def __init__(
        self,
        dimension: int,
        degree: int,
        coefficients: Mapping[tuple[int, ...], int | Fraction] | None = None,
    ) -> None:
        if isinstance(dimension, bool) or not isinstance(dimension, int):
            raise TypeError("dimension must be an integer")
        if dimension < 1:
            raise ValueError("dimension must be positive")
        if isinstance(degree, bool) or not isinstance(degree, int):
            raise TypeError("degree must be an integer")
        if degree < 0:
            raise ValueError("degree must be nonnegative")
        if coefficients is None:
            coefficients = {}
        if not isinstance(coefficients, Mapping):
            raise TypeError("coefficients must be a mapping of index tuples to rationals")
        terms: list[tuple[tuple[int, ...], Fraction]] = []
        for indices, coefficient in coefficients.items():
            _index_tuple(indices, dimension=dimension, degree=degree)
            exact = _exact_scalar(coefficient)
            if exact:
                terms.append((indices, exact))
        if degree > dimension and terms:
            raise ValueError("nonzero exterior forms cannot exceed the ambient degree")
        object.__setattr__(self, "dimension", dimension)
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "coefficients", tuple(sorted(terms)))

    def as_dict(self) -> dict[tuple[int, ...], Fraction]:
        """Return a fresh mutable dictionary of the canonical coefficients."""
        return dict(self.coefficients)

    def __add__(self, other: Form) -> Form:
        if not isinstance(other, Form):
            return NotImplemented
        _same_grade(self, other)
        result = self.as_dict()
        for indices, coefficient in other.coefficients:
            result[indices] = result.get(indices, Fraction(0)) + coefficient
        return Form(self.dimension, self.degree, result)

    def __neg__(self) -> Form:
        return Form(self.dimension, self.degree, {i: -c for i, c in self.coefficients})

    def __sub__(self, other: Form) -> Form:
        if not isinstance(other, Form):
            return NotImplemented
        return self + (-other)

    def __mul__(self, scalar: int | Fraction) -> Form:
        exact = _exact_scalar(scalar)
        return Form(self.dimension, self.degree, {i: exact * c for i, c in self.coefficients})

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | Fraction) -> Form:
        exact = _exact_scalar(scalar)
        if not exact:
            raise ZeroDivisionError("division of a form by zero")
        return self * (1 / exact)


def _same_dimension(left: Form, right: Form) -> None:
    if left.dimension != right.dimension:
        raise ValueError("forms have different ambient dimensions")


def _same_grade(left: Form, right: Form) -> None:
    _same_dimension(left, right)
    if left.degree != right.degree:
        raise ValueError("forms have different degrees")


def _inversion_sign(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    """Sign of placing ``left`` followed by ``right`` in sorted order."""
    inversions = sum(i > j for i in left for j in right)
    return -1 if inversions % 2 else 1


def basis_form(
    dimension: int, indices: tuple[int, ...], coefficient: int | Fraction = 1
) -> Form:
    """Construct one sorted exterior basis monomial with an exact coefficient."""
    if not isinstance(indices, tuple):
        raise TypeError("indices must be a tuple")
    return Form(dimension, len(indices), {indices: coefficient})


def wedge(left: Form, right: Form) -> Form:
    """The associative, graded-anticommutative exterior product ``left ∧ right``."""
    _same_dimension(left, right)
    output_degree = left.degree + right.degree
    if output_degree > left.dimension:
        return Form(left.dimension, output_degree)
    result: dict[tuple[int, ...], Fraction] = {}
    for left_basis, left_coefficient in left.coefficients:
        for right_basis, right_coefficient in right.coefficients:
            # Repeated basis directions force the exterior product to zero.
            if not set(left_basis).isdisjoint(right_basis):
                continue
            basis = tuple(sorted((*left_basis, *right_basis)))
            value = left_coefficient * right_coefficient * _inversion_sign(left_basis, right_basis)
            result[basis] = result.get(basis, Fraction(0)) + value
    return Form(left.dimension, output_degree, result)


def hodge_star(form: Form, *, orientation: int = 1) -> Form:
    """Euclidean Hodge dual for an oriented *orthonormal* basis.

    For the standard orientation, ``basis ∧ star(basis)`` equals the positive
    volume form.  ``orientation=-1`` reverses the volume form, hence reverses
    every Hodge dual.  For a k-form, ``star(star(a))=(-1)^(k(n-k))*a`` when
    the *same* orientation is used on both calls.
    """
    if type(orientation) is not int or orientation not in (-1, 1):
        raise ValueError("orientation must be +1 or -1")
    if form.degree > form.dimension:
        raise ValueError("Hodge star is only defined for degrees 0 through n")
    result: dict[tuple[int, ...], Fraction] = {}
    for basis, coefficient in form.coefficients:
        basis_set = set(basis)
        complement = tuple(i for i in range(form.dimension) if i not in basis_set)
        result[complement] = orientation * _inversion_sign(basis, complement) * coefficient
    return Form(form.dimension, form.dimension - form.degree, result)


def inner_product(left: Form, right: Form) -> Fraction:
    """Coefficient-wise Euclidean inner product in the orthonormal wedge basis.

    Forms of different degrees are orthogonal by convention; ambient spaces
    of different dimensions cannot be compared.
    """
    _same_dimension(left, right)
    if left.degree != right.degree:
        return Fraction(0)
    right_terms = right.as_dict()
    return sum((coefficient * right_terms.get(basis, 0)
                for basis, coefficient in left.coefficients), Fraction(0))


def norm_squared(form: Form) -> Fraction:
    """Nonnegative rational squared Euclidean norm, without taking a root."""
    return inner_product(form, form)


__all__ = ["Form", "basis_form", "wedge", "hodge_star", "inner_product", "norm_squared"]
