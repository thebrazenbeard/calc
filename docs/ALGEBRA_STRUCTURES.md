# Small, verifiable algebraic structures

This module is a **worked mathematical laboratory**, not a general-purpose
computer algebra system. It implements four sharply delimited structures and
exposes their assumptions so a reader can verify the axioms directly.

## 1. Finite groups

For a finite set \(G=\{0,\ldots,n-1\}\) and a binary operation \(\cdot\),
the `FiniteGroup` constructor accepts an \(n\times n\) Cayley table and checks:

1. Closure: every entry belongs to \(G\).
2. Associativity: \((a\cdot b)\cdot c=a\cdot(b\cdot c)\) for every triple.
3. A two-sided identity \(e\), with \(e\cdot a=a\cdot e=a\).
4. A two-sided inverse for every element.

It exposes `mul(a,b)`, `inverse(a)`, `element_order(a)`,
`cyclic_generators()`, `is_abelian`, and `identity`. Construction costs
\(O(n^3)\); the input limit is \(1\le n\le128\). That bound avoids implying
this cubic checker handles arbitrary large finite groups efficiently.

```python
from calc.structures import FiniteGroup

z4 = FiniteGroup([[(a + b) % 4 for b in range(4)] for a in range(4)])
assert z4.identity == 0
assert z4.inverse(1) == 3
assert z4.cyclic_generators() == (1, 3)
```

## 2. Prime fields

`PrimeField(p)` models \(\mathbb{F}_p=\mathbb{Z}/p\mathbb{Z}\) for a
**prime** \(2\le p<2^{64}\). A deterministic set of Miller–Rabin witnesses
is used for this *bounded* range; larger moduli are intentionally rejected.
Methods `add`, `sub`, `mul`, `div`, `inv`, `pow` accept integer residues,
normalizing them modulo \(p\). Inverting zero raises `ZeroDivisionError`.
Python's arbitrary-precision integers preserve exact modular operations.

```python
from calc.structures import PrimeField

f7 = PrimeField(7)
assert f7.inv(3) == 5          # 3 × 5 ≡ 1 (mod 7)
assert f7.div(4, 3) == 6     # 4 × 5 ≡ 6 (mod 7)
```

The multiplicative inverse rule is not extended to composite moduli, because
\(\mathbb{Z}/n\mathbb{Z}\) has nonzero zero-divisors when \(n\) is composite.

## 3. Boolean algebra

`boolean_join`, `boolean_meet`, `boolean_complement`, and `boolean_xor`
provide \(\lor\), \(\land\), \(\neg\), and exclusive-or for Python
`bool` operands. Both operands are always validated—even where an `or`/`and`
short-circuit could otherwise hide an invalid input. Python integers `0`/`1`
are rejected despite their truth values.

## 4. Tropical algebra (min-plus semiring)

Over \(\mathbb{R}\cup\{+\infty\}\), define:

\[
a\oplus b=\min(a,b),\qquad a\otimes b=a+b.
\]

Its additive identity is \(+\infty\) and its multiplicative identity is
ordinary zero. The implementation accepts Python `int`/`float`, and treats
\(+\infty\) as the absorbing element of multiplication. NaN and
\(-\infty\) are excluded. Floating-point overflow is signaled; approximate
floating-point results are **not** exact real-number algebra.

```python
from calc.structures import tropical_add, tropical_multiply

assert tropical_add(3, -2) == -2
assert tropical_multiply(3, -2) == 1
```

## Verification and limits

`tests/test_structures.py` verifies group axioms and a noncommutative \(S_3\)
composition table, rejects nonassociative/invalid tables, exhaustively checks
the field distributive law for small primes, tests De Morgan's laws, and tests
tropical distributivity on a finite grid including infinity. These tests are
strong examples; they are not formal proofs of arbitrary Python inputs.

Useful background references:

- MIT OpenCourseWare, [18.703 Modern Algebra](https://ocw.mit.edu/courses/18-703-modern-algebra-spring-2013/): groups, rings, fields.
- MIT OpenCourseWare, [18.06 Linear Algebra](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/): related exact linear structures.
- J. Richter-Gebert, B. Sturmfels and T. Theobald, [First Steps in Tropical Geometry](https://arxiv.org/abs/math/0306366) (2003): introduction to the min-plus tropical semiring. The research reference supports the broad topic, not this code's correctness.
- [Miller–Rabin deterministic witnesses (64-bit range)](https://miller-rabin.appspot.com/): reference implementation literature and witness results; the repository's own tests are distinct from external witness qualification.

No external code has been copied.
