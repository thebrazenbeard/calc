# Exact tensor components and contractions

`calc.tensors` provides a small dependency-free **component algebra** over
the rational numbers. A `Tensor` is a dense, immutable, row-major array in a
*fixed basis*. It supports ordinary tensor products, exact index permutations,
explicit contractions, symmetric/alternating parts, and matrix products.

This is a bounded worked example of multilinear algebra, **not** a complete
tensor-calculus system. In particular, a rectangular array by itself does not
specify whether an index is covariant or contravariant, how it transforms
under a change of coordinates, or which metric/connection should be used.

## Mathematical contract

For component arrays (A_{i_1\dots i_r}) and (B_{j_1\dots j_s}), the outer
product has components

\[
(A\otimes B)_{i_1\dots i_r j_1\dots j_s}
   = A_{i_1\dots i_r}B_{j_1\dots j_s}.
\]

`permute_axes(order)` reorders axes explicitly. For a selected pair of axes
`a != b` of the same size, `contract(a,b)` sums the diagonal entries:

\[
C_{\text{remaining indices}} = \sum_{k=0}^{d-1}
    T_{\dots\,i_a=k\,\dots\,i_b=k\,\dots}.
\]

For rank two, this is the usual trace. Applying contraction to adjacent
matrix-product axes of (A\otimes B) gives ((AB)_{ij}=\sum_k A_{ik}B_{kj}).
Contraction here is an explicit array operation: mathematically, the natural
contraction pairs a vector-space factor with its **dual**. The API leaves that
variance decision to its caller instead of inventing a metric.

`symmetrize(a,b)` computes ((T+\tau_{ab}T)/2), and
`antisymmetrize(a,b)` computes ((T-\tau_{ab}T)/2). Their sum recovers (T).

## Executable example

```python
from fractions import Fraction
from calc.tensors import Tensor, identity_tensor

u = Tensor((2,), (Fraction(1, 2), 3))
v = Tensor((2,), (2, -1))
outer = u.tensor_product(v)
assert outer.shape == (2, 2)
assert outer[0, 1] == Fraction(-1, 2)
assert outer.contract(0, 1)[()] == Fraction(-2)

A = Tensor((2, 2), (1, 2, 3, 4))
B = Tensor((2, 2), (5, 6, 7, 8))
assert A @ B == Tensor((2, 2), (19, 22, 43, 50))
assert A @ B == A.tensor_product(B).contract(1, 2)
assert identity_tensor(2) @ A == A
assert A == A.symmetrize(0, 1) + A.antisymmetrize(0, 1)
```

Zero-dimensional *axes* are excluded. A **rank-zero** tensor is valid and has
shape `()` and exactly one scalar component, accessed with `scalar[()]`.
Integer and `Fraction` values are accepted; floats, complex coefficients and
booleans are rejected. Every output remains rational, subject to size fences.

## Resource and correctness boundaries

- At most **6 axes**, each dimension between **1 and 32**, and at most
  **4096 dense entries**. Outer products that exceed these bounds are refused
  *before* allocation. Iterable inputs are read with a finite sentinel so an
  unbounded generator cannot hang construction.
- Coefficient numerator and denominator bit lengths are each capped at
  **4096 bits**, including after operations.
- Arithmetic does not broadcast dimensions, infer Einstein repeated-index
  contractions, create a metric, or silently raise/lower indices.
- `matrix_multiply` accepts only rank-two tensors with compatible dimensions.
  It avoids allocating a potentially much larger four-axis outer product.
- A component identity is not evidence of tensorial coordinate behavior;
  this module does not implement change-of-basis laws or derivatives.
- Product and contraction routines are straightforward bounded loops, not
  optimized high-performance numerical kernels. Rational arithmetic can be
  substantially more expensive than binary floating point.

## Verification

`tests/test_tensors.py` checks 20 independent witnesses, including all six
axis orderings in a nontrivial three-axis array, explicit component-wise
contraction oracles, Kronecker-delta identities, product associativity and
bilinearity, symmetric/alternating decomposition, and resource rejection.

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_tensors.py -v
```

## Original references (linked, not copied)

- MIT OpenCourseWare, *18.101 Analysis II*, [Lecture 17: Multilinear
  Algebra](https://ocw.mit.edu/courses/18-101-analysis-ii-fall-2005/35663a54a97807582f11e26ceff7be1e_lecture17.pdf).
  Introduces multilinear tensor products.
- Brian Conrad, Stanford Math 210A, [Tensor algebras, tensor pairings, and
  duality](https://math.stanford.edu/~conrad/210APage/handouts/tensoralg.pdf).
  Explains the abstract construction and role of dual pairings.

No source code or long text passages are copied from these references.
