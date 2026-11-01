# Independent hostile review: Wave 8

**Review date:** 2026-11-01  
**Scope:** Chebyshev–Lobatto collocation, classical umbral calculus, finite categories/functors/natural transformations, rational homology, and bounded exact tensor arrays.  
**Independence:** The reviewer inspected the modules and wrote only `tests/test_hostile_wave8.py` and this report. Production modules and GitHub refs were not changed by the reviewer.

## Verdict

**19/19 new independent mathematical tests pass** on Python 3.12 in the shared workspace. No production defect was confirmed by the reviewed examples or bounded resource probes. This is **partial verification** of specific documented contracts, not a proof of all supported inputs or a certification of mathematical completeness. The source modules expose substantial scope limits.

## Evidence matrix

| Area | Adversarial / independent oracle | Observed result | Risk left unproved |
| --- | --- | --- | --- |
| Spectral collocation | Degree-7 polynomial interpolated on eight-degree Chebyshev–Lobatto grid; analytic derivative of polynomial; differentiation matrix maps `x` to `1` and constants to `0`, multiple sizes | Pass | Floating conditioning and error bounds at degree 128 and extreme data magnitudes; no PDE spectral convergence claim |
| Umbral calculus | Signed Stirling matrix inverse identity through degree 12; 72 seeded rational coefficient transformations; falling/rising factorial shift identities; third-order raw forward difference checked via independent binomial expansion; Newton interpolation on a negative fractional grid | Pass | Arbitrary Sheffer sequences, q-umbral operators and precision limits outside exact rational domain |
| Finite categories | Diamond-poset paths compose to the same arrow; monotone endofunctors and naturality squares; cyclic group of order three admits identity and inversion functors but no natural transformation between them; invalid unital magma rejected | Pass | Enumerating arbitrary categories, decidability beyond finite tables, or equivalence of differently presented categories |
| Rational homology | Four-edge square cycle has Betti `(1,1)`; after adding filling face `(1,0,0)`; 60 seeded rational chain complexes satisfy `d₁d₂=0`, vanish in homology and obey Euler–Poincaré; dependent rational row rank | Pass | Generating chain complexes from geometric topology; torsion and integral-coefficient homology; arbitrary sparse large complexes |
| Exact tensor algebra | Kronecker delta contraction recovers source matrix; inverse axis permutations across a nonsymmetric 2×3×4 shape; symmetric/antisymmetric decomposition and zero trace of skew part; exact matrix product associativity | Pass | Abstract covariant/contravariant type tracking, metric geometry, automatic differentiation and huge tensor contractions |
| Resource and invalid-input handling | Unbounded generators challenged in tensor, interpolation, barycentric and poset APIs; strict dimension/order maximums tested | Pass | Time/memory behavior of all values within their documented caps; high-bit rational intermediate growth |

## Reproduction

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_wave8.py -v
```

All 19 cases passed. The tests are deterministic, have no network or optional dependencies, and require only Python's standard library. Source examples also have their own separate test modules; this review intentionally introduces mathematically different witnesses instead of merely reasserting their fixture values.

## Source and licensing boundary

The mathematical identities tested here are classical and independently formulated. The repository's subject guides cite mathematical sources (including NIST DLMF and academic lecture notes), but a citation is not evidence that every upstream license permits copying. This review did not copy third-party source code, import external projects, or grant a license to the repository.

## Remaining uncertainty

- A finite sample of analytic identities and random exact matrices does not formally prove algebraic law compliance for every bounded input.
- The collocation module returns ordinary floating-point estimates. It does not promise certified enclosures, stable recovery from all dynamic ranges, or a full PDE spectral solver.
- The finite-category representation is explicitly tabulated and capped; it is not an implementation of arbitrary large categories.
- Homology over the rationals detects Betti ranks but **not integer torsion**. Topology must be supplied as chain groups and boundary matrices.
- Tensor components are arrays in a fixed basis; no index variance or manifold connection is inferred.

**Status:** 19 passing adversarial checks, no confirmed unresolved defect in this bounded review. Integration and publication remain the parent agent's responsibility.
