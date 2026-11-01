# Finite categories, functors, and natural transformations

`calc.categories` is an **executable finite example**, not a solver for general category theory. It uses names and complete finite tables. Every category is validated at construction; invalid examples fail immediately.

## Mathematics and conventions

A category consists of objects, typed arrows, an identity arrow on each object, and associative composition. For (f:A\to B) and (g:B\to C), this module stores the composite as `composition[("g", "f")]`, meaning (g\circ f). The table must specify the result for **every** composable pair and **no** noncomposable pairs. The category constructor verifies:

\[
f\circ 1_A=f=1_B\circ f,
\qquad h\circ(g\circ f)=(h\circ g)\circ f.
\]

A *covariant functor* (F:\mathcal C\to\mathcal D) assigns target objects and arrows, preserving their endpoints, identities, and composition. A *natural transformation* \(\eta:F\Rightarrow G\) assigns a component \(\eta_A:F(A)\to G(A)\) to each source object, satisfying the naturality square:

\[
G(f)\circ \eta_A=\eta_B\circ F(f)
\quad (f:A\to B).
\]

The code checks that equality for **every** source arrow. Names are opaque labels; the library does not infer equations or invent missing composites.

## Worked example: a chain poset as a thin category

Every finite partial order gives a *thin* category: one arrow \(a\to b\) exactly when \(a\le b\). Its composition exists precisely by transitivity. The relation must be explicitly reflexive and transitively closed.

```python
from calc.categories import poset_category, FiniteFunctor, NaturalTransformation

source = poset_category(
    ("a", "b", "c"),
    {("a", "a"), ("b", "b"), ("c", "c"),
     ("a", "b"), ("b", "c"), ("a", "c")},
)
target = poset_category(("x", "y"),
                        {("x", "x"), ("y", "y"), ("x", "y")})

def monotone_functor(mapping):
    arrows = {f.name: target.hom(mapping[f.source], mapping[f.target])[0]
              for f in source.arrows}
    return FiniteFunctor(source, target, mapping, arrows)

F = monotone_functor({"a": "x", "b": "x", "c": "y"})
G = monotone_functor({"a": "x", "b": "y", "c": "y"})
eta = NaturalTransformation(
    F, G, {"a": target.identity("x"),
           "b": target.hom("x", "y")[0],
           "c": target.identity("y")},
)
assert len(source.arrows) == 6
assert eta.components["b"] == target.hom("x", "y")[0]
```

The `monoid_category` helper likewise builds a one-object category from a finite associative multiplication table with identity. It checks the monoid laws through the category axioms. This captures simple algebraic examples, not all monoid theory.

## Scope and safeguards

- At most **32 objects** and **64 arrows**, with a complete composition table. Associativity checks have worst-case \(O(m^3)\) time for \(m\) arrows (at most \(64^3\) candidate triples).
- Iterable inputs are read only up to one past the relevant limit, so an infinite iterator is rejected without being exhausted. Mapping inputs are treated as finite supplied objects; callers should not supply side-effecting/dynamically mutating mappings.
- Constructors copy and expose tables as read-only views. All result objects are finite combinatorial structures; no numeric approximation is involved.
- No general co/limits, Yoneda computations, adjunctions, enriched/\(\infty\)-categories, set-theoretic universe management, symbolic inference, or automatic category enumeration.
- A passing finite table check proves the listed axioms **for that exact table**, not for an unlimited family of categories.

## Sources and reuse

- Emily Riehl, [*Category Theory in Context*](https://math.jhu.edu/~eriehl/context/) (author's book page; especially the chapters defining categories, functors, and natural transformations).
- David I. Spivak, [*Category Theory for Scientists*](https://ocw.mit.edu/courses/18-s996-category-theory-for-scientists-spring-2013/pages/open-textbook/), MIT OpenCourseWare, Section 4 (basic category theory and natural transformations).
- MIT OCW, [Algebraic Topology I: Lecture 3 — Categories, Functors, Natural Transformations](https://ocw.mit.edu/courses/18-905-algebraic-topology-i-fall-2016/2a313c70afa3c392a4d0cb36c6558b3e_MIT18_905F16_lec3.pdf).

These are **linked references**, not imported code or copied textbook content. Attribution does not imply the repository has a license to redistribute those texts.
