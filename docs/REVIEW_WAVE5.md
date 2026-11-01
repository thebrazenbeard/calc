# Wave 5 independent hostile verification

**Date:** 2026-11-01  
**Reviewer:** separate Wave 6 review agent, read-only on production  
**Scope:** `matrix_calc.py`, `grobner.py`, `time_scales.py`, `fractional.py`, `intervals.py`, packaged `topic_manifest.json`, `evidence.py`, and the evidence command.  
**Regression suite:** `tests/test_hostile_wave5.py`, independent of the implementation tests.  
**Result:** four classes of deterministically reproduced false-overflow defects found; parent patched production and all retained regression tests now pass. Other probes did not falsify the examined contracts. **This is not a proof of general mathematical correctness.**

## Reproduce the adversarial verification

From the repository root:

```bash
PYTHONPATH=src python -m unittest -v tests.test_hostile_wave5
```

The suite includes 24 named tests, hundreds of deterministic seeded interval
checks, 35 independently derived univariate ideal witnesses, exact rational
counterexamples and public-package contract checks. It uses Python's standard
library only. `unittest` reports each individual test separately; the seeded
checks each use `subTest` to expose a failing witness.

## Findings disclosed and fixed during review

| Severity | Finding and exact witness | Prior result | Post-fix regression |
| --- | --- | --- | --- |
| Moderate; high confidence | Finite-grid delta/nabla derivative of values `(-1e308,+1e308)` on times `(0,1e308)` equals exactly `2`. | Raised `ValueError` because the finite-value numerator overflowed before division. | `test_finite_difference_exactly_two_despite_subtraction_overflow` |
| Moderate; high confidence | Central matrix gradient/JVP of `f(X)=1e308*X[0][0]` at `X=((0,),)` with step `1` equals `1e308`, which is finite. | Raised `ValueError` on an overflowing `(+1e308)-(-1e308)` subtraction. | `test_gradient_of_linear_map_not_rejected_on_intermediate_overflow`; `test_jvp_of_linear_map_not_rejected_on_intermediate_overflow` |
| Moderate; high confidence | Delta integral on `(0,2,4)` of left samples `(1e308,-1e308)` equals zero; corresponding nabla right samples also sum to zero. | Rejected the individual `±2e308` products although the total was finite. | `test_delta_integral_cancellation_of_overflowed_terms`; `test_nabla_integral_cancellation_of_overflowed_terms` |
| Moderate; high confidence | `frobenius_inner(((1e308,1e308),), ((2.,-2.),))` equals zero, and the analogous trace derivative dot product also cancels. | `math.fsum` received `+inf` and `-inf` because products overflowed first. | `test_frobenius_cancellation_of_individually_overflowing_products`; `test_trace_gradient_cancellation_of_overflowed_dot_products` |

All four are *false positive rejection* defects, not evidence of an accepted
wrong finite answer. Parent production fixes use a normal binary64 fast path
and exact rational fallback in exceptional overflow cases. The regression
suite asserts the precise finite answers without expecting exceptions.

## Independent mathematical witnesses

1. **Matrix calculus:** Recomputed the shape and orientation of
   `∇_X tr(A X B)=(B A)^T` on a nonsquare, signed integer `A:2×3` and
   `B:4×2`; all 12 entries matched a hand-computed matrix. Checked an exact
   `Fraction` Frobenius inner product and exact `2X` gradient.
2. **Gröbner bases:** From `x²=y, xy=1`, directly derived the reduced lex
   basis `{x−y², y³−1}`; implementation matched, including leading monomials.
   Independently established gcd-based witnesses for 35 seeded univariate
   pairs sharing a linear factor, and tested rational coefficients. For the
   graded lex basis, all original generators and all Buchberger S-polynomials
   reduce to zero. Tested exponent, power and Buchberger pair budgets.
3. **Time scales:** Checked finite-grid discrete fundamental theorem on an
   irregular grid `(1,2,4,7)`: forward and backward telescoping weighted sums
   both reproduce `f(7)−f(1)=48` for `f(t)=t²`. The overflow regressions above
   explicitly probe valid finite output at extreme binary64 inputs.
4. **Fractional calculus:** Checked `D_C^{1/2} t` at `x=4` against
   `4/√π`, `D_{RL}^{1/2} 1` against `1/(2√π)` and ordinary first derivative
   of `t³` at `x=2` against `12`. Independent Grünwald–Letnikov constant
   sample with 8192 terms is close to `1/√π`, rather than Caputo's zero.
   Singular Riemann–Liouville endpoint is rejected, as specified.
5. **Intervals:** Compared addition, subtraction, all four multiplication
   endpoint products, and division (when zero is excluded) against exact
   `Fraction.from_float` rational extrema for 180 deterministic interval pairs,
   including very small and near-overflow values. Checked overflow enclosure,
   binary64 `0.1×0.1`, division-by-zero and disjoint intersections.

## Manifest, package and scholarly source audit

- Confirmed `importlib.resources.files("calc")` contains parseable JSON that
  equals `load_evidence()`. Confirmed `pyproject.toml` explicitly lists
  `topic_manifest.json` as `calc` package data. **This is a source/install
  contract check, not an independent wheel-build test.**
- Enumerated each manifest topic and verified unique slugs; valid linked
  source IDs in `docs/SOURCES.md`; accessible guide and test paths; declared
  implemented symbols can actually be imported; `example_tested` implies
  `implemented`, bindings and test paths; all topics explicitly say
  `formal_proof: false`.
- Verified the JSON CLI path by comparing
  `python -m calc evidence time-scales --json` directly against
  `get_evidence("time-scales")`, not a manually copied expected snapshot.
- The user-suggested [Jon Krohn ML-foundations repository](https://github.com/jonkrohn/ML-foundations)
  declares **MIT** on its GitHub repository listing; `docs/SOURCES.md` maps it
  to `S04` as **link-only**, with copying/attribution cautions. We did not
  conduct a file-by-file similarity or all-notebook copyright analysis.
- The official [MIT OpenCourseWare terms](https://ocw.mit.edu/pages/privacy-and-terms-of-use/)
  identify **CC BY-NC-SA 4.0** with attribution, noncommercial and share-alike
  obligations. The inventory treats MIT OCW as reference links rather than
  permissively licensed material to import. An [open time-scale paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11579191/)
  is marked open access in its publisher-hosted record, but that does not
  automatically clear unrelated linked papers. Original-source license checks
  must still precede any actual copied code/text/media.

## Residual risk and confidence

**High confidence:** the four false-overflow failures existed pre-fix and
the exact finite answers pass current production; finite-grid discrete
identities, lex elimination example, rational polynomial witnesses, and the
binary64 outward-interval containment checks.

**Medium confidence:** fractional closed forms in tested parameter regions
and grlex normal-form behavior. The latter includes an algorithmic property
(Buchberger zero reductions) but does not formally prove termination or
correctness for every ideal and budget combination.

**Not established:** broad symbolic algebra equivalence, interval functions
beyond four basic arithmetic operations, dense/general time scales, fractional
operator convergence over classes of nonsmooth functions, numerical error
bounds for finite differences, wheel contents on all build backends, every
external document's detailed license or independent CI of every topic branch.
The Gröbner exponent and coefficient/input budgets are educational limits, not
a comprehensive CPU/memory guarantee on adversarial polynomial evaluation.

**Disposition:** The named defects are fixed in production and covered by
strict regressions. No confirmed uncorrected mathematical falsifier was
observed within this scoped review. The rest of the mathematical atlas still
requires specialist audit and progressively stronger independent oracles.
