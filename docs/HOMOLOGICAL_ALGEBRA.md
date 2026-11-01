# Finite chain complexes and rational homology

This module computes **Betti numbers for supplied finite chain complexes of
rational vector spaces**. It is an entry point into homological algebra, not a
general-purpose algorithm for extracting topology from a point cloud,
triangulating a manifold, or computing integral homology with torsion.

## Mathematical model

A chain complex over the rational field \(\mathbb{Q}\) consists of
finite-dimensional vector spaces and boundary maps

\[
\cdots \xrightarrow{d_3} C_2 \xrightarrow{d_2} C_1
\xrightarrow{d_1} C_0 \rightarrow 0,
\qquad d_{n-1}\circ d_n = 0.
\]

Cycles are \(\ker d_n\); boundaries are \(\operatorname{im}d_{n+1}\).
Since every boundary is a cycle, homology is the quotient

\[
H_n = \ker d_n / \operatorname{im} d_{n+1},\qquad
\beta_n = \dim C_n-\operatorname{rank}d_n-\operatorname{rank}d_{n+1}.
\]

The code checks the chain identity *exactly* before accepting any complex.
It uses rational Gaussian elimination, so the ranks and Betti numbers have no
floating-point roundoff. Every matrix column corresponds to a basis element
of \(C_n\) and every row to a basis element of \(C_{n-1}\). Matrix orientations
matter: \(d_n\) is **rows = lower-degree generators, columns = current-degree
generators**.

The Euler–Poincaré identity provides a useful independent check:

\[
\sum_n(-1)^n\dim C_n=\sum_n(-1)^n\beta_n.
\]

## Executable examples

```python
from calc.homology import ChainComplex, rational_rank

# An oriented interval, edge e01 from v0 to v1:
# d1(e01) = v1-v0, so beta0=1 and beta1=0.
interval = ChainComplex((2, 1), {1: ((-1,), (1,))})
assert interval.betti_numbers() == (1, 0)

# A CW chain model for the circle: one vertex, one loop, zero boundary.
circle = ChainComplex((1, 1))
assert circle.betti_numbers() == (1, 1)

# An oriented triangle before filling its interior: beta1=1.
# Edges are ordered e01, e12, e02.
edges = ((-1, 0, -1), (1, -1, 0), (0, 1, 1))
unfilled = ChainComplex((3, 3), {1: edges})
assert unfilled.betti_numbers() == (1, 1)

# Adding the 2-simplex [012] gives d2([012])=e12-e02+e01.
filled = ChainComplex((3, 3, 1), {
    1: edges,
    2: ((1,), (1,), (-1,)),
})
assert filled.betti_numbers() == (1, 0, 0)
assert filled.euler_characteristic() == filled.homology_euler_characteristic() == 1

# Independent exact rank example.
assert rational_rank(((1, 2), (2, 4))) == 1
```

The circle example supplies the *cellular* boundary matrices directly.
It does not infer that those matrices describe a circle from geometric data.

## Scope and constraints

- Supports chain degrees \(0\) through at most \(7\) and each chain group
  dimension at most \(32\).
- Every supplied coefficient must be an `int` or `fractions.Fraction`;
  `float`, `bool`, strings, and complex numbers are rejected.
- Coefficient numerator and denominator each have a **2048-bit input cap**.
  Intermediate exact elimination can still grow; this is suitable for small
  educational examples, not large sparse homology computations.
- Omitted boundary matrices mean correctly shaped zero maps. For explicit
  maps, the shape must match `(dim C_(n-1), dim C_n)` even when dimensions
  are zero; an empty outer tuple represents a zero-row matrix.
- The `ChainComplex` object stores an immutable rational snapshot. Results
  are integer ranks/Betti numbers, not generators, bases, representatives,
  barcode persistence, or chain homotopies.
- Integer-coefficient homology, including torsion such as
  \(\mathbb Z/2\mathbb Z\), requires additional tools: rational Betti numbers
  alone cannot reveal torsion.

## References and provenance

- [MIT OpenCourseWare 18.905 Algebraic Topology I, Lecture Notes (2016)](https://ocw.mit.edu/courses/18-905-algebraic-topology-i-fall-2016/2dad3dd8354992d877fbeca4049817ee_MIT18_905F16_lecture_notes.pdf) — chain condition, cycles, boundaries and homology.
- [Allen Hatcher, *Algebraic Topology*, Chapter 2: Homology](https://pi.math.cornell.edu/~hatcher/AT/ATchapters.html) — underlying homology definitions, cellular and simplicial examples.

These sources are cited as mathematical references. No source code or textbook
content has been copied into the repository.
