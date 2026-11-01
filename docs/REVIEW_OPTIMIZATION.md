# Independent hostile review: scalar optimization

Review subject: `src/calc/optimization.py` and its public functions `bisection`, `safeguarded_newton`, and `golden_section_minimize`, inspected locally on 2026-11-01. Reviewer: `/root/optimization_hostile` (native ChatGPT subagent). This is **internal independent source/test review**, not human review, formal proof, or cloud CI.

## Verdict

**PARTIAL — focused contracts pass after one confirmed repair.** Bounded scalar algorithms implement useful mathematical examples, with correctly stated continuity/unimodality assumptions and honest convergence exceptions. No global minimizer or certified root isolation is established. Current exact GitHub branch state and cloud CI were outside this review's readback.

## Confirmed finding and resolution

**MEDIUM; high confidence; fixed locally by the integration owner.** `golden_section_minimize` originally calculated `b-a` directly and used it to generate golden-section interior points. Even with individually finite endpoints, `[-1e308, 1e308]` has an unrepresentably large float width. The resulting arithmetic sent `-inf`, `inf`, and `nan` coordinates into the caller's objective and nonetheless returned `-1e308` for a constant function.

Direct pre-fix runtime reproduction:

```python
seen = []
f = lambda x: (seen.append(x), 0.0)[1]
golden_section_minimize(f, -1e308, 1e308, max_iter=30)
# Returned -1e308; sampled inputs included -inf, inf, nan.
```

The integrator added an upfront finite-width check to reject this input before objective evaluation. `tests/test_hostile_optimization.py::HostileMinimizationContracts.test_finite_endpoints_never_expose_nonfinite_query` now verifies no nonfinite query is allowed; the test passes on the repaired code. The test permits a controlled `ValueError` or a finite safe calculation, so it does not over-prescribe an implementation strategy.

## Requirement-to-evidence matrix

| Risk or requirement | Independent oracle or falsifier | Result |
| --- | --- | --- |
| Bracket sign and midpoint safety | `f(x)=x/1e308` across ±`1e308` has root zero | PASS |
| Discontinuity rejection | Nonzero sign jump cannot be a root; `f_tol` tighter than jump | PASS |
| Even-multiplicity root limitation | `(x-.5)^2` has no opposite-sign endpoints | Expected `ValueError` |
| Newton derivative fallback | Derivative raises `ZeroDivisionError` or returns zero at start; cubic/linear root known analytically | PASS |
| Nonfinite callback output | Interior `nan` after valid endpoints | Expected `ValueError` |
| Golden-section interior minimizer | Unique quadratic minimum at `-1.25` | PASS |
| Golden-section boundary minimum | Strictly decreasing linear function on finite interval | PASS |
| Extreme finite endpoints | Prior nonfinite objective queries from width overflow | Reproduced before fix; regression PASS after fix |
| Invalid boolean API arguments | Strict integer/real parameter contracts | PASS |

Focused execution: `PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_optimization.py -v` — **9 tests, exit 0** on Python 3.12 after the integrator's repair.

## Residual limitations and challenges

> A sign-changing bracket does not prove continuity. A discontinuity with jump smaller than `f_tol` can still masquerade as an approximate root. The docs disclose this, and the algorithm must not be advertised as certification.

> Golden-section search depends on unimodality. Finite samples do not prove unimodality or a global minimum, and flat roundoff-dominated objectives can shift the returned argmin.

- The function callbacks can have side effects or cost arbitrary time; they are normal user-provided Python callables, not a sandboxed formula interpreter.
- Very large/small scales and subnormal values can exhaust IEEE-754 resolution. `ConvergenceError` is a useful rejection rather than proof of nonexistence of a mathematical root.
- Tests are deterministic mathematical witnesses, not formal verification, exhaustive search, or a comparative high-precision package validation.
- The upstream SciPy [bisection reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html) explicitly states continuity/opposite-sign preconditions; its [bounded scalar optimization reference](https://docs.scipy.org/doc/scipy/reference/optimize.minimize_scalar-bounded.html) provides a separate production-method comparison. These are cited references; no upstream code was reused.
