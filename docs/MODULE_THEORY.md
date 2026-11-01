# Finite module theory over residue rings

This module computes exact, bounded examples of modules over
\(R=\mathbb Z/n\mathbb Z\). It models the free module \(M=R^d\), generated
submodules \(N\le M\), and the induced quotient module \(M/N\).
All entries are integer residue classes; there are no floating-point
approximations. This is a small executable model of the algebraic definitions,
**not** a general-purpose classification or presentation algorithm for modules.

## Module operations

For \(x,y\in M\), \(r\in R\), and coordinates \(i=1,\ldots,d\):

\[
(x+y)_i = (x_i+y_i)\bmod n,\qquad
(rx)_i=(r x_i)\bmod n.
\]

The zero vector and additive inverse are componentwise, too. These obey the
module axioms, including distributivity \(r(x+y)=rx+ry\),
\((r+s)x=rx+sx\), and \((rs)x=r(sx)\).

Because the scalar ring can have zero divisors, \(M\) need **not** be a
vector space. For instance, \(2\cdot(2)=0\) in \(\mathbb Z/4\mathbb Z\)
although neither factor is zero. Consequently, vector-space dimension/rank
theorems over fields must not be indiscriminately reused for these modules.

```python
from calc.modules import ResidueFreeModule

m = ResidueFreeModule(6, 2)         # (Z/6Z)^2, with 36 elements
assert m.add((5, 1), (4, 8)) == (3, 3)
assert m.scalar_multiply(4, (2, 3)) == (2, 0)
```

## Generated submodules

Given generators \(g_1,\ldots,g_k\), the generated submodule is their
finite residue-linear span:

\[
N=\Big\{\sum_{j=1}^{k}r_jg_j:\ r_j\in\mathbb Z/n\mathbb Z\Big\}.
\]

`generated_submodule` computes this with a bounded breadth-first closure
under addition by the generators. The finite group contains inverses because
each generator has finite additive order. It is exactly the smallest
submodule containing the generators, without requiring a basis or assuming
the ring is a field.

```python
n = m.generated_submodule((2, 0), (0, 3))
assert n.order == 6
assert n.index == 6
assert n.contains((4, 3))
assert not n.contains((1, 0))
```

The example has \(N=\{0,2,4\}\times\{0,3\}\), a six-element submodule
of the 36-element ambient module.

## Quotient modules

For \(N\le M\), a coset is \(x+N=\{x+n:n\in N\}\). The quotient carries
induced operations

\[
(x+N)+(y+N)=(x+y)+N,\qquad r(x+N)=(rx)+N.
\]

These are well-defined because adding any representative from \(N\) does
not change the coset. A quotient instance enumerates each disjoint coset
exactly once and assigns it a stable integer index. This index is merely an
implementation label; it is not a canonical module isomorphism.

```python
quotient = n.quotient()
assert quotient.order == 6
assert quotient.class_of((1, 0)) == quotient.class_of((3, 0))
zero = quotient.zero_class
first = quotient.class_of((1, 0))
assert quotient.add(first, quotient.neg(first)) == zero
```

By the finite additive coset theorem, \(|M/N|\cdot|N|=|M|\). The tests
verify this and the representative-independent quotient operations for
multiple moduli, including composite values.

## Resource and validation contract

- \(2\le n\le4096\), \(1\le d\le8\), and \(n^d\le4096\).
- At most 16 supplied submodule generators.
- Input vectors are read only through \(d+1\) items, so an infinite
  iterable is rejected without hanging.
- Coordinates and scalar multipliers must be Python integers; floats,
  booleans, malformed vectors, and ambiguous implicit coercions are rejected.
- The construction materializes full submodules and cosets in memory; it is
  intended for teaching, finite checks, and small combinatorial examples.
- All operations work over a single common modulus and fixed free-module
  rank. Arbitrary finitely presented modules, exact sequences, tensor
  products of modules, homological resolutions, invariant-factor
  classification and infinite modules remain outside this implementation.

## Verification

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_modules.py -v
```

The tests independently enumerate residue-linear spans by coefficient
tuples, compare the generated submodules with that oracle for 100 seeded
small instances, exhaustively check the module axioms on several full
spaces, and verify quotient class operations across all pairs of ambient
vectors for a small non-field scalar ring.

Tests establish the stated bounded implementation behavior, not the module
axioms for every possible abstract ring or module.

## Primary references

1. MIT OpenCourseWare, *Algebra II Student Notes*, Lecture 19, “Modules over
   a Ring”, especially the module over \(\mathbb Z/n\mathbb Z\) example and
   the quotient construction:
   https://ocw.mit.edu/courses/res-18-012-algebra-ii-student-notes-spring-2022/mit18_702s22_full_lec.pdf
2. Stanford University, Ben Lynn, *Submodules and Quotient Modules*, including
   the induced scalar action on a quotient:
   https://crypto.stanford.edu/pbc/notes/commalg/submodule.html

These are research/teaching references only. No upstream source code is
copied into this implementation, and linked materials retain their own
licenses.
