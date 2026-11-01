"""Exact rational integration and conditional expectation on finite atom sets.

This module models *only* finite discrete measure spaces, not arbitrary
sigma-algebras, Lebesgue integration, or constructive Radon--Nikodym
derivatives on continuous domains. Its zero-measure conventions are explicit.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice

Rational = int | Fraction
_MAX_ATOMS = 128
_MAX_BITS = 2048
_MAX_RESULT_BITS = 65536
_MAX_LABEL = 128


def _rational(value: object, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an exact integer or Fraction")
    result = Fraction(value)
    if (result.numerator.bit_length() > _MAX_BITS or
            result.denominator.bit_length() > _MAX_BITS):
        raise ValueError(f"{name} exceeds the {_MAX_BITS}-bit input limit")
    return result


def _result(value: Fraction, name: str) -> Fraction:
    if (value.numerator.bit_length() > _MAX_RESULT_BITS or
            value.denominator.bit_length() > _MAX_RESULT_BITS):
        raise ValueError(f"{name} exceeds the exact result size limit")
    return value


def _atom(value: object) -> str:
    if not isinstance(value, str) or not value or len(value) > _MAX_LABEL:
        raise TypeError(f"atom labels must be nonempty strings of at most {_MAX_LABEL} characters")
    return value


def _observable(values: Mapping[str, Rational], atoms: tuple[str, ...],
                name: str) -> tuple[Fraction, ...]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{name} must map every atom to a rational value")
    if set(values) != set(atoms):
        raise ValueError(f"{name} must contain exactly the measure's atom labels")
    return tuple(_rational(values[atom], f"{name}[{atom!r}]") for atom in atoms)


@dataclass(frozen=True, slots=True, init=False)
class FiniteMeasure:
    """A finite labelled atom space with nonnegative rational weights.

    ``FiniteMeasure({'a': 1, 'b': Fraction(1, 2)})`` represents
    ``mu({a})=1, mu({b})=1/2``. Labels are sorted canonically; the original
    input mapping is never retained, so subsequent mutation cannot affect the
    measure. Zero-mass atoms are permitted and remain in the domain.
    """

    atoms: tuple[str, ...]
    weights: tuple[Fraction, ...]

    def __init__(self, weights: Mapping[str, Rational]):
        if not isinstance(weights, Mapping):
            raise TypeError("weights must be a finite mapping of string atoms to rational masses")
        if not 1 <= len(weights) <= _MAX_ATOMS:
            raise ValueError(f"weights must have 1 through {_MAX_ATOMS} atoms")
        atoms = tuple(sorted(_atom(atom) for atom in weights))
        values = tuple(_rational(weights[atom], f"mass[{atom!r}]") for atom in atoms)
        if any(value < 0 for value in values):
            raise ValueError("measure masses must be nonnegative")
        object.__setattr__(self, "atoms", atoms)
        object.__setattr__(self, "weights", values)

    @property
    def total_mass(self) -> Fraction:
        """Compute ``mu(Omega)`` exactly."""
        return _result(sum(self.weights, Fraction(0)), "total mass")

    @property
    def support(self) -> tuple[str, ...]:
        """Names of positive-mass atoms, excluding null atoms."""
        return tuple(atom for atom, mass in zip(self.atoms, self.weights) if mass)

    def mass(self, event: Iterable[str]) -> Fraction:
        """Mass of a finite event, with duplicate labels treated setwise."""
        if isinstance(event, (str, bytes)):
            raise TypeError("event must be an iterable of atom labels, not a string")
        try:
            sequence = tuple(islice(iter(event), _MAX_ATOMS + 1))
        except TypeError as exc:
            raise TypeError("event must be an iterable of atom labels") from exc
        if len(sequence) > _MAX_ATOMS:
            raise ValueError("event exceeds atom limit")
        chosen = {_atom(label) for label in sequence}
        if not chosen.issubset(self.atoms):
            raise ValueError("event contains an atom outside this measure's domain")
        return _result(sum((mass for atom, mass in zip(self.atoms, self.weights)
                            if atom in chosen), Fraction(0)), "event mass")

    def integrate(self, observable: Mapping[str, Rational]) -> Fraction:
        """Finite Lebesgue integral ``sum_x f(x)*mu({x})`` over the atoms."""
        values = _observable(observable, self.atoms, "observable")
        return _result(sum((value * mass for value, mass in zip(values, self.weights)),
                           Fraction(0)), "integral")

    def normalize(self) -> FiniteMeasure:
        """Return the induced probability distribution; require positive mass."""
        total = self.total_mass
        if not total:
            raise ValueError("zero-total-mass measure cannot be normalized")
        normalized = tuple(_result(weight / total, "normalized mass")
                           for weight in self.weights)
        # This is a validated algebraic result, not new untrusted user input.
        # It may legitimately need more bits than each original mass because
        # normalization combines denominators from multiple atoms.
        result = object.__new__(FiniteMeasure)
        object.__setattr__(result, "atoms", self.atoms)
        object.__setattr__(result, "weights", normalized)
        return result

    def expectation(self, observable: Mapping[str, Rational]) -> Fraction:
        """Expectation with respect to the normalized distribution."""
        total = self.total_mass
        if not total:
            raise ValueError("expectation requires positive total mass")
        return _result(self.integrate(observable) / total, "expectation")

    def variance(self, observable: Mapping[str, Rational]) -> Fraction:
        """Normalized finite variance ``E[(X-E[X])^2]`` without rounding."""
        values = _observable(observable, self.atoms, "observable")
        total = self.total_mass
        if not total:
            raise ValueError("variance requires positive total mass")
        mean = self.expectation(dict(zip(self.atoms, values)))
        numerator = sum((weight * (value - mean)**2
                         for weight, value in zip(self.weights, values)), Fraction(0))
        return _result(numerator / total, "variance")

    def density_of(self, target: FiniteMeasure) -> dict[str, Fraction]:
        """Return ``d(target)/d(self)`` for finite absolutely continuous pairs.

        If a reference atom has zero mass, the target must also give it zero
        mass. The density is conventionally set to zero there, since its
        value on a null atom is not uniquely determined.
        """
        if not isinstance(target, FiniteMeasure):
            raise TypeError("target must be a FiniteMeasure")
        if target.atoms != self.atoms:
            raise ValueError("measures must have the same underlying atom set")
        density: dict[str, Fraction] = {}
        for atom, ref_mass, target_mass in zip(self.atoms, self.weights, target.weights):
            if not ref_mass:
                if target_mass:
                    raise ValueError("target is not absolutely continuous with respect to reference")
                density[atom] = Fraction(0)
            else:
                density[atom] = _result(target_mass / ref_mass, "density")
        return density

    def conditional_expectation(
        self, observable: Mapping[str, Rational], partition: Iterable[Iterable[str]],
    ) -> dict[str, Fraction]:
        """Conditional expectation onto the sigma-algebra of a finite partition.

        Each nonempty cell must appear exactly once in a disjoint cover of
        all atoms. The function is constant on each cell, the weighted mean
        if the cell has positive measure. On a null cell, its value is chosen
        to be zero; conditional expectations are defined only almost surely.
        """
        values = _observable(observable, self.atoms, "observable")
        try:
            cells = tuple(islice(iter(partition), _MAX_ATOMS + 1))
        except TypeError as exc:
            raise TypeError("partition must be an iterable of cells") from exc
        if not 1 <= len(cells) <= _MAX_ATOMS:
            raise ValueError("partition has an invalid number of cells")
        known = set(self.atoms)
        seen: set[str] = set()
        result: dict[str, Fraction] = {}
        weights = dict(zip(self.atoms, self.weights))
        observations = dict(zip(self.atoms, values))
        for cell in cells:
            if isinstance(cell, (str, bytes)):
                raise TypeError("each partition cell must be an iterable of labels")
            try:
                members = tuple(islice(iter(cell), _MAX_ATOMS + 1))
            except TypeError as exc:
                raise TypeError("each partition cell must be iterable") from exc
            if not 1 <= len(members) <= _MAX_ATOMS:
                raise ValueError("each partition cell must be nonempty and bounded")
            chosen = {_atom(item) for item in members}
            if len(chosen) != len(members) or not chosen.issubset(known):
                raise ValueError("partition cells contain duplicates or unknown atoms")
            if seen.intersection(chosen):
                raise ValueError("partition cells must be disjoint")
            seen.update(chosen)
            mass = sum((weights[atom] for atom in chosen), Fraction(0))
            numerator = sum((weights[atom] * observations[atom] for atom in chosen), Fraction(0))
            average = _result(numerator / mass, "conditional mean") if mass else Fraction(0)
            result.update({atom: average for atom in chosen})
        if seen != known:
            raise ValueError("partition must cover every atom")
        return {atom: result[atom] for atom in self.atoms}


__all__ = ["FiniteMeasure"]
