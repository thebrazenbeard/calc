# Exact rational Clifford / geometric algebra

`calc.clifford` implements a small, transparent Clifford algebra `Cl(p,q)`
with exact rational coefficients. It is deliberately bounded: the goal is to
check and explore defining identities, not to replace an industrial geometric
algebra system.

## Convention

Given ordered generators `e_0,...,e_(n−1)`, specify a signature such as
`(+1,+1,−1)`. The geometric product is associative, distributive, and obeys

\[
e_i e_j+e_j e_i=2\eta_{ij},\qquad
\eta_{ij}=\operatorname{diag}(1,1,-1).
\]

An integer bitmask encodes a positively ordered blade; for example `0b101`
means `e_0 e_2`. The scalar unit has mask `0`. For a homogeneous grade-`r`
blade `A_r` and grade-`s` blade `B_s`, the exterior product retains grade
`r+s` from the geometric product. The **left contraction** retains grade
`s−r` only for `r≤s` (and is zero otherwise). Contraction conventions vary
across references, so the particular choice is part of this API.

Reversion multiplies grade `k` by `(-1)^{k(k−1)/2}`; grade involution by
`(-1)^k`. Clifford conjugation combines the two.

## Running examples

```python
from fractions import Fraction
from calc.clifford import CliffordAlgebra

cl = CliffordAlgebra((1, 1, -1))
e0, e1, e2 = (cl.generator(i) for i in range(3))

assert e0 * e0 == cl.scalar(1)
assert e2 * e2 == cl.scalar(-1)
assert e0 * e1 == -(e1 * e0)
assert e0.wedge(e1) == cl.basis(0b011)
assert (e0.wedge(e1)).reverse() == -e0.wedge(e1)

v = Fraction(1, 3) * e0 + 2 * e2
assert v.grade(1) == v
assert v.reversion_quadratic() == Fraction(1, 9) - 4
```

`CliffordAlgebra.multivector({0: a, 1: b, 3: c})` builds a sparse general
multivector with scalar, vector, and bivector contributions, respectively.
The `.as_dict()` method returns a fresh dictionary of exact `Fraction`
coefficients. Algebra objects and multivectors themselves are immutable.

## Scope and guarantees

- **Exactness:** only Python `int` or `fractions.Fraction` inputs; booleans,
  binary floats, NaNs, and complex inputs are rejected, with no silent
  rounding. Arithmetic is rational, not arbitrary-precision real analysis.
- **Metrics:** diagonal orthogonal *unit* signatures with each basis vector
  square `+1` or `−1`; no null directions, nonorthogonal metrics, or
  coordinate transformations.
- **Resources:** dimension `1..8`, hence at most 256 basis blades; each
  result coefficient's numerator and denominator are bounded to 4096 bits.
  This is not intended for large symbolic Clifford computations.
- **Quadratic form:** `A.reversion_quadratic()` returns the scalar part of
  `reverse(A) * A`. It equals the nonnegative sum of coefficient squares for
  all-positive signatures; it is **not a positive norm** for indefinite
  metrics. No inverse, square root, spinor representation, or exponential is
  implemented.
- **Meaning of "implemented":** tested selected operations, not the whole
  mathematical theory or a theorem-proving certificate.

## Verification

```sh
PYTHONPATH=src python -m unittest discover -s tests -p test_clifford.py -v
```

The tests exercise the defining signature, pairwise anticommutation,
association of every ordered triple of basis blades in `Cl(2,1)` (512
triples), graded wedge symmetry, left contraction, reversion, exact rational
operations, and input/resource limits. Exhaustive small-basis tests are not a
proof of every possible rational multivector expression.

## Independent mathematical references

- [Lundholm and Svensson, *Clifford algebra, geometric algebra, and applications*](https://www.math.lmu.de/~lundholm/clifford.pdf), lecture notes — multiplication, gradings, involutions, contractions.
- [Dorst and Mann, *Geometric algebra: a computational framework for geometrical applications*](https://staff.fnwi.uva.nl/l.dorst/clifford/dorst-mann-I.pdf), University of Amsterdam — geometric product and blade multiplication tables.
- [PyGAE `clifford` project](https://github.com/pygae/clifford) — further implementation ecosystem; this module copies no source code and does not require that dependency.
