"""A navigable taxonomy, not a claim to exhaustive mathematical coverage.

Status refers to this repository: ``implemented`` means executable examples
are provided; ``documented`` means the subject is indexed for study, not
that a complete algorithm or theorem prover exists here.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Topic:
    slug: str
    title: str
    family: str
    status: str = "documented"


_DATA: dict[str, tuple[str, ...]] = {
    "calculus-foundations": (
        "limits", "continuity", "differential", "integral", "real-analysis",
        "complex-analysis", "asymptotic-analysis", "measure-theory",
    ),
    "calculus-multivariable": (
        "multivariable", "vector", "tensor", "differential-geometry",
        "differential-forms", "exterior", "manifold", "riemannian",
    ),
    "calculus-dynamics": (
        "ordinary-differential-equations", "partial-differential-equations",
        "calculus-of-variations", "functional", "optimal-control",
        "hamiltonian", "geometric-mechanics", "distributional",
    ),
    "calculus-generalizations": (
        "fractional", "stochastic", "ito", "stratonovich", "malliavin",
        "discrete", "finite-difference", "q-calculus", "time-scales",
        "umbral", "nonstandard-analysis", "constructive-analysis",
        "interval", "automatic-differentiation", "nonsmooth", "subdifferential",
    ),
    "calculus-computation": (
        "numerical-differentiation", "numerical-integration", "spectral-methods",
        "finite-elements", "finite-volumes", "integral-transforms",
        "symbolic-calculus", "lambda-calculus", "relational-calculus",
    ),
    "algebra-foundations": (
        "elementary-algebra", "polynomial", "linear-algebra", "matrix",
        "boolean-algebra", "abstract-algebra", "group-theory", "ring-theory",
        "field-theory", "module-theory", "universal-algebra",
    ),
    "algebra-structures": (
        "commutative-algebra", "noncommutative-algebra", "lie-algebra",
        "tensor-algebra", "exterior-algebra", "clifford-algebra", "geometric-algebra",
        "homological-algebra", "tropical-algebra", "operator-algebra",
        "representation-theory", "category-theory", "algebraic-geometry",
        "computational-algebra", "hopf-algebra", "differential-algebra",
    ),
}

_TITLES = {
    "ito": "Itô calculus",
    "stratonovich": "Stratonovich calculus",
    "malliavin": "Malliavin calculus",
    "q-calculus": "q-calculus",
    "lambda-calculus": "Lambda calculus (computer science)",
    "relational-calculus": "Relational calculus (database logic)",
    "ordinary-differential-equations": "Ordinary differential equations",
    "partial-differential-equations": "Partial differential equations",
}

_IMPLEMENTED = {
    "differential", "multivariable", "finite-difference",
    "numerical-differentiation", "numerical-integration", "polynomial",
    "linear-algebra", "matrix", "ordinary-differential-equations",
}

TOPICS: tuple[Topic, ...] = tuple(
    Topic(slug, _TITLES.get(slug, slug.replace("-", " ").title()), family,
          "implemented" if slug in _IMPLEMENTED else "documented")
    for family, slugs in _DATA.items()
    for slug in slugs
)


def list_topics(family: str | None = None) -> tuple[Topic, ...]:
    """List topics, optionally filtered by exact family identifier."""
    return tuple(t for t in TOPICS if family is None or t.family == family)


def get_topic(slug: str) -> Topic:
    """Retrieve a topic by stable slug; raise KeyError if unknown."""
    for topic in TOPICS:
        if topic.slug == slug:
            return topic
    raise KeyError(slug)

