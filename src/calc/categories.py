"""Finite, explicitly enumerated categories, functors and natural transformations.

Composition is written ``(after, before)``: ``compose(g, f)`` means ``g ∘ f``.
All data are finite combinatorial labels. This is *not* a general category
theory prover, a category of all sets, or an implementation of large categories.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import islice
from types import MappingProxyType

_MAX_OBJECTS = 32
_MAX_ARROWS = 64


def _bounded_tuple(values: Iterable, name: str, cap: int) -> tuple:
    try:
        items = tuple(islice(iter(values), cap + 1))
    except TypeError as exc:
        raise TypeError(f"{name} must be iterable") from exc
    if len(items) > cap:
        raise ValueError(f"{name} exceeds {cap} entries")
    return items


def _label(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or len(value) > 80 or value.strip() != value:
        raise ValueError(f"{name} must be a nonempty, trimmed label of at most 80 characters")
    return value


def _unique_labels(values: Iterable[str], name: str, cap: int) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be an iterable of labels, not a string")
    try:
        items = _bounded_tuple(values, name, cap)
    except ValueError as exc:
        raise ValueError(f"{name} count must be in 1..{cap}") from exc
    if not 1 <= len(items) <= cap:
        raise ValueError(f"{name} count must be in 1..{cap}")
    for item in items:
        _label(item, name)
    if len(items) != len(set(items)):
        raise ValueError(f"{name} labels must be unique")
    return items


def _copy_map(mapping: Mapping, name: str) -> dict:
    if not isinstance(mapping, Mapping):
        raise TypeError(f"{name} must be a mapping")
    caps = {"identities": _MAX_OBJECTS, "composition": _MAX_ARROWS ** 2,
            "on_objects": _MAX_OBJECTS, "on_arrows": _MAX_ARROWS,
            "components": _MAX_OBJECTS}
    if name in caps and len(mapping) > caps[name]:
        raise ValueError(f"{name} has too many entries")
    return dict(mapping)


@dataclass(frozen=True, slots=True)
class Arrow:
    """Named arrow ``source → target`` in a finite category."""

    name: str
    source: str
    target: str

    def __post_init__(self) -> None:
        for field in ("name", "source", "target"):
            _label(getattr(self, field), field)


@dataclass(frozen=True, slots=True, init=False)
class FiniteCategory:
    """Small category with a complete, validated composition table.

    ``composition[(g, f)]`` is required for every composable pair of arrows
    ``f:A→B`` and ``g:B→C`` and forbidden for noncomposable pairs. Constructor
    verifies identity laws and all composable triples' associativity.
    """

    objects: tuple[str, ...]
    arrows: tuple[Arrow, ...]
    identities: Mapping[str, str]
    composition: Mapping[tuple[str, str], str]

    def __init__(
        self,
        objects: Iterable[str],
        arrows: Iterable[Arrow],
        identities: Mapping[str, str],
        composition: Mapping[tuple[str, str], str],
    ) -> None:
        objs = _unique_labels(objects, "objects", _MAX_OBJECTS)
        if isinstance(arrows, (str, bytes)):
            raise TypeError("arrows must be an iterable of Arrow")
        morphs = _bounded_tuple(arrows, "arrows", _MAX_ARROWS)
        if not 1 <= len(morphs) <= _MAX_ARROWS or not all(isinstance(m, Arrow) for m in morphs):
            raise ValueError(f"arrows must contain 1..{_MAX_ARROWS} Arrow values")
        names = [m.name for m in morphs]
        if len(set(names)) != len(names):
            raise ValueError("arrow names must be unique")
        by_name = {m.name: m for m in morphs}
        for arrow in morphs:
            if arrow.source not in objs or arrow.target not in objs:
                raise ValueError(f"arrow {arrow.name} has an unknown endpoint")

        ids = _copy_map(identities, "identities")
        if set(ids) != set(objs):
            raise ValueError("identities must specify exactly one arrow per object")
        for obj, name in ids.items():
            if name not in by_name or (by_name[name].source, by_name[name].target) != (obj, obj):
                raise ValueError(f"identity of {obj} must be an endomorphism at that object")
        if len(set(ids.values())) != len(objs):
            raise ValueError("distinct objects must have distinct identity arrows")

        comp = _copy_map(composition, "composition")
        expected = {(g.name, f.name)
                    for f in morphs for g in morphs if f.target == g.source}
        if set(comp) != expected:
            missing = len(expected - set(comp))
            extra = len(set(comp) - expected)
            raise ValueError(f"composition table incomplete or ill-typed ({missing} missing, {extra} extra)")
        for (after, before), value in comp.items():
            g, f = by_name[after], by_name[before]
            h = by_name.get(value)
            if h is None or h.source != f.source or h.target != g.target:
                raise ValueError(f"composition {after} ∘ {before} has the wrong result type")
        for f in morphs:
            if comp[(f.name, ids[f.source])] != f.name:
                raise ValueError("right identity law fails")
            if comp[(ids[f.target], f.name)] != f.name:
                raise ValueError("left identity law fails")
        # At most 64^3 candidate triples, bounded independently of input data.
        for f in morphs:
            for g in morphs:
                if f.target != g.source:
                    continue
                gf = comp[(g.name, f.name)]
                for h in morphs:
                    if g.target == h.source:
                        left = comp[(h.name, gf)]
                        right = comp[(comp[(h.name, g.name)], f.name)]
                        if left != right:
                            raise ValueError(f"associativity fails for {h.name}, {g.name}, {f.name}")

        object.__setattr__(self, "objects", objs)
        object.__setattr__(self, "arrows", morphs)
        object.__setattr__(self, "identities", MappingProxyType(ids))
        object.__setattr__(self, "composition", MappingProxyType(comp))

    def arrow(self, name: str) -> Arrow:
        for arrow in self.arrows:
            if arrow.name == name:
                return arrow
        raise KeyError(name)

    def identity(self, obj: str) -> str:
        return self.identities[obj]

    def hom(self, source: str, target: str) -> tuple[str, ...]:
        if source not in self.objects or target not in self.objects:
            raise KeyError("source and target must belong to the category")
        return tuple(m.name for m in self.arrows
                     if m.source == source and m.target == target)

    def compose(self, after: str, before: str) -> str:
        """Return ``after ∘ before``; reject ill-typed pairs."""
        try:
            return self.composition[(after, before)]
        except KeyError as exc:
            raise ValueError("arrows are unknown or not composable in this category") from exc


@dataclass(frozen=True, slots=True, init=False)
class FiniteFunctor:
    """A covariant functor between two validated finite categories."""

    source: FiniteCategory
    target: FiniteCategory
    on_objects: Mapping[str, str]
    on_arrows: Mapping[str, str]

    def __init__(self, source: FiniteCategory, target: FiniteCategory,
                 on_objects: Mapping[str, str], on_arrows: Mapping[str, str]) -> None:
        if not isinstance(source, FiniteCategory) or not isinstance(target, FiniteCategory):
            raise TypeError("source and target must be finite categories")
        objects = _copy_map(on_objects, "on_objects")
        arrows = _copy_map(on_arrows, "on_arrows")
        if set(objects) != set(source.objects) or any(o not in target.objects for o in objects.values()):
            raise ValueError("functor object mapping must be total and target valid objects")
        if set(arrows) != {a.name for a in source.arrows}:
            raise ValueError("functor arrow mapping must be total")
        for a in source.arrows:
            try:
                mapped = target.arrow(arrows[a.name])
            except KeyError as exc:
                raise ValueError("functor maps an arrow to an unknown target arrow") from exc
            if mapped.source != objects[a.source] or mapped.target != objects[a.target]:
                raise ValueError("functor does not preserve arrow endpoints")
        for obj in source.objects:
            if arrows[source.identity(obj)] != target.identity(objects[obj]):
                raise ValueError("functor does not preserve identity")
        for (g, f), gf in source.composition.items():
            if arrows[gf] != target.compose(arrows[g], arrows[f]):
                raise ValueError("functor does not preserve composition")
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "on_objects", MappingProxyType(objects))
        object.__setattr__(self, "on_arrows", MappingProxyType(arrows))


@dataclass(frozen=True, slots=True, init=False)
class NaturalTransformation:
    """A family ``η_X:F(X)→G(X)`` satisfying ``G(f)∘η_X=η_Y∘F(f)``."""

    source_functor: FiniteFunctor
    target_functor: FiniteFunctor
    components: Mapping[str, str]

    def __init__(self, source_functor: FiniteFunctor,
                 target_functor: FiniteFunctor,
                 components: Mapping[str, str]) -> None:
        if not isinstance(source_functor, FiniteFunctor) or not isinstance(target_functor, FiniteFunctor):
            raise TypeError("both endpoints must be finite functors")
        f, g = source_functor, target_functor
        if f.source != g.source or f.target != g.target:
            raise ValueError("functors must have the same source and target categories")
        terms = _copy_map(components, "components")
        if set(terms) != set(f.source.objects):
            raise ValueError("transformation must give exactly one component per source object")
        for obj, name in terms.items():
            try:
                arrow = f.target.arrow(name)
            except KeyError as exc:
                raise ValueError("natural-transformation component is an unknown arrow") from exc
            if arrow.source != f.on_objects[obj] or arrow.target != g.on_objects[obj]:
                raise ValueError("component has incorrect endpoints")
        for a in f.source.arrows:
            left = f.target.compose(g.on_arrows[a.name], terms[a.source])
            right = f.target.compose(terms[a.target], f.on_arrows[a.name])
            if left != right:
                raise ValueError(f"naturality fails for arrow {a.name}")
        object.__setattr__(self, "source_functor", f)
        object.__setattr__(self, "target_functor", g)
        object.__setattr__(self, "components", MappingProxyType(terms))


def poset_category(objects: Iterable[str], less_equal: Iterable[tuple[str, str]]) -> FiniteCategory:
    """Build the thin category of a finite poset with full ``≤`` relation.

    ``less_equal`` must already contain reflexive and transitive closure.
    Posets are antisymmetric; arbitrary preorders are not accepted here.
    """
    objs = _unique_labels(objects, "objects", _MAX_OBJECTS)
    if isinstance(less_equal, (str, bytes)):
        raise TypeError("less_equal must be an iterable of pairs")
    try:
        relation = set(_bounded_tuple(less_equal, "poset relation", _MAX_ARROWS))
    except (TypeError, ValueError) as exc:
        if isinstance(exc, ValueError):
            raise ValueError("poset has too many morphisms for finite category cap") from exc
        raise TypeError("less_equal must contain hashable pairs") from exc
    if not all(isinstance(pair, tuple) and len(pair) == 2 and
               pair[0] in objs and pair[1] in objs for pair in relation):
        raise ValueError("relation must contain ordered pairs of the provided objects")
    if not all((obj, obj) in relation for obj in objs):
        raise ValueError("poset relation must be reflexive")
    if any(a != b and (b, a) in relation for a, b in relation):
        raise ValueError("poset relation must be antisymmetric")
    if any((a, c) not in relation for a, b in relation
           for d, c in relation if b == d):
        raise ValueError("poset relation must be transitive")
    if len(relation) > _MAX_ARROWS:
        raise ValueError("poset has too many morphisms for finite category cap")
    indices = {name: i for i, name in enumerate(objs)}
    names = {(a, b): f"m_{indices[a]}_{indices[b]}" for a, b in relation}
    arrows = [Arrow(names[(a, b)], a, b) for a, b in sorted(relation)]
    identities = {obj: names[(obj, obj)] for obj in objs}
    comp = {(names[(b, c)], names[(a, b)]): names[(a, c)]
            for a, b in relation for d, c in relation if b == d}
    return FiniteCategory(objs, arrows, identities, comp)


def monoid_category(elements: Iterable[str], identity: str,
                    multiplication: Mapping[tuple[str, str], str]) -> FiniteCategory:
    """Build a one-object category from a fully tabulated finite monoid.

    ``multiplication[(g, f)]`` specifies the categorical composite ``g∘f``.
    FiniteCategory checks that the operation is associative and has ``identity``.
    """
    names = _unique_labels(elements, "monoid elements", _MAX_ARROWS)
    if identity not in names:
        raise ValueError("monoid identity must be an element")
    return FiniteCategory(("•",), (Arrow(name, "•", "•") for name in names),
                          {"•": identity}, multiplication)


__all__ = ["Arrow", "FiniteCategory", "FiniteFunctor", "NaturalTransformation",
           "poset_category", "monoid_category"]
