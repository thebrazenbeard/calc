# Independent hostile review: Wave 11

**Review date:** 2026-11-01  
**Target:** `calc.rings`, `calc.representations`, `calc.affine_geometry`,
`calc.asymptotics`, and `calc.modules`.  
**Owned regression file:** `tests/test_hostile_wave11.py`.

This is an adversarial **model-authored executable review**, performed
separately from the implementation agents. It is not external peer review,
formal verification, or a proof of correctness for the entire mathematical
subjects named by the modules. No upstream source code was copied.

## Findings and independent evidence

| Surface | Falsifying challenge / independent oracle | Result |
| --- | --- | --- |
| `Z/nZ` ring ideals | For every modulus 2–23 and every single generator, enumerate all scalar multiples and annihilating residues directly; compare ideals and cardinalities | Pass |
| Ideal arithmetic | For modulus 2–14 and all pairs of generators, independently enumerate ideal sums, intersections, and products as finite sets | Pass |
| Chinese remainder theorem | Recover five preselected integers from congruences mod 4, 9, 5, 7; reject noncoprime moduli | Pass |
| Finite representations | Independently generate cyclic left-regular permutation matrices for orders 1–9 and compare every entry and character; check trivial/regular pairing, direct sum, intertwiners | Pass |
| Nontrivial characters | The rational sign representation of `C2` has self-pairing 1 and pairing 0 with the trivial representation; incompatible `C3` generator action must fail | Pass |
| Affine algebraic sets | Rational unit-circle grid points and exact Jacobian, membership in `(x²−1)` versus an unrelated polynomial, and a grid with no rational point on `x²−2` | Pass |
| Asymptotic analysis | 110 seeded exact monomial comparisons; independently classify big-O, little-o, Theta, and the quotient limit at both positive infinity and zero-plus | Pass |
| Finite modules | 96 generated-submodule cases independently enumerated by every scalar coefficient tuple; quotient coset partition, representative independence and action laws | Pass |
| Resource/type boundaries | Reject oversized module rank, booleans and nonintegral coordinates, nonintegral scalar action, and wrong vector arity | Pass |

The independent expected values are derived from elementary algebraic
identities, explicit small finite-set enumeration, and analytic polynomial
identities. They do not simply call another method of the same production
class to compute the oracle.

## Scope challenges and unresolved evidence

> A successful rational grid search cannot certify an entire algebraic
> variety, and an empty grid does **not** imply no algebraic points exist.

> Characters and matrix identities for finite rational group representations
> do **not** establish complete classification, irreducible decomposition, or
> complex/continuous representation theory.

> The exact asymptotic classifier uses leading terms of nonzero rational
> functions at selected approaches. It makes no claim about arbitrary
> oscillatory, transcendental, or multivariable limits.

> The residue-module quotient is an explicitly bounded finite partition,
> not a structure-theorem implementation for arbitrary modules over rings.

The review also does not establish a GitHub-hosted CI pass, real browser
interaction, a mathematical completeness guarantee, or external independent
assessment. Source references for definitions and conventions are linked in
the five subject documents, including MIT OCW, NIST DLMF, SageMath's official
polynomial ideal documentation, and Stanford lecture notes; links serve as
attribution, not permission to copy their contents.

## Reproduction

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_wave11.py -v
PYTHONPATH=src python -m unittest discover -s tests -q
```

At the review checkpoint on Python 3.12:

- **12/12 hostile review methods passed**.
- **915/915 repository tests passed** in the shared workspace.

**Verdict:** No confirmed unresolved defect within the above bounded
examples. That means the tested examples passed, not that the full named
mathematical theories have been implemented or verified.
