# Wave 2: independent hostile mathematical and software review

**Review subject:** local `src/calc/{symbolic,autodiff,geometry,stochastic,structures}.py`, their own test modules and documentation, the 19 `docs/topics/*.md` lessons, and the wave-one shared core. The only GitHub baseline known to this reviewer is the reported wave-one `main` commit `81aebf196b2b742a3da7242e3357ea97c0baf4c8`; the wave-two files were inspected in the shared workspace. This is **internal subagent review**, not a second provider, a formally verified mathematics library, or GitHub CI qualification.

**Verdict: PARTIAL for the full 'every calculus/algebra' objective; PASS for the nine direct independent wave-two contract tests after two repaired defects.** Original claims remain restricted to their named algorithms. All changes to production code were made by the parent integrator, not the reviewer.

## Confirmed defects and disposition

1. **Medium severity, high confidence — nonfinite geometric metric output (fixed locally).** Previously, `geometry.metric_inner` silently returned `inf` for a valid finite metric and finite coordinate vectors whose result overflowed binary64: `metric_inner(lambda _: ((1.0,),), (0.0,), (1e308,), (1e308,))`. Returning a nonfinite scalar without rejecting or qualifying it compromises downstream numeric reasoning. A direct runtime probe observed `inf` before the fix; the parent then added output validation. The new `test_metric_inner_rejects_overflow_instead_of_returning_infinity` **passes after the fix**. The negative test itself was not executed against the original implementation. Exact repository publication of this fix requires separate GitHub readback.
2. **Medium severity, high confidence — non-advancing stochastic timestamps (fixed locally).** `euler_maruyama(drift=1,diffusion=0,t0=1e20,y0=0,dt=1e-9,steps=2,increments=(0,0))` previously returned times `(1e20,1e20,1e20)` while increasing state: the positive step is smaller than local floating-point time spacing. The reviewer added an expectation that the solver raise `ValueError`; it **failed first (1/9 contracts failing)**. The parent added per-step forward-time representability checks, and the focused suite then **passed 9/9**. This prevents silently fabricating time advancement.

## Independent mathematical witnesses

New tests in `tests/test_wave2_contracts.py` check:

| Mathematical invariant or threat | Independent witness | Result after patches |
| --- | --- | --- |
| Product, chain and logarithm derivatives | Analytic `d(x²sin x + ln x)/dx` at four positive points, singular `x=0` | PASS |
| Safe symbolic parsing | Reject import invocation, attribute access, comprehension, lambda, indexing, f-string, illegal arity, excessive exponent | PASS |
| Forward-mode multivariate AD | Analytic gradient of `exp(xy)+sin x` | PASS |
| Metric transformation | Affine Jacobian `A`, expected pullback `AᵀA=[[5,-1],[-1,10]]` | PASS |
| Exterior derivative | Identity `d(df)=0` for smooth polynomial-plus-trigonometric scalar | PASS (finite tolerance) |
| Numeric overflow | Finite metric/vector input whose inner product exceeds binary64 range must be rejected | PASS after fix |
| Itô discrete solver deterministic specialization | Zero diffusion, constant drift, supplied irrelevant Brownian samples | PASS |
| Time representability | Tiny positive time step at huge epoch must not silently produce duplicate timestamps | PASS after fix |
| Finite prime-field/tropical algebra | Exhaustive 11-residue distributivity; min-plus distributivity including `+∞` | PASS |

Fresh focused command: `PYTHONPATH=src python -m unittest discover -s tests -p test_wave2_contracts.py -v` — Python 3.12, 9 tests, exit 0 after both repairs. An initial focused run with the stochastic regression produced exit 1, with a `ValueError not raised` failure before the repair. This counterfactual sensitivity matters more than mere coverage count.

## Remaining hostile challenges

> **An indexed branch is not an implemented mathematical paradigm.** The catalogue has scores of terms, but most are research indices. A few examples of Itô SDEs, fractional power integrals, and exterior derivatives cannot support a claim of complete Itô, fractional, or exterior calculus. The package's `implemented` label is explained as *having a runnable example*; consumers must not read it as theory coverage.

> **Numerical checks are not formal error certificates.** Christoffel symbols and exterior derivatives use finite differences; Brownian paths use pseudorandom samples; standard quadrature can alias oscillations. Smoothness, grid resolution, floating-point spacing, and cancellation matter. Existing docs acknowledge most of these limits. No interval arithmetic or formal verification was executed.

> **The GitHub description promises separate topic branches.** The wave-two `docs/BRANCHES.md` is a branch design, not remote branch evidence. The 19 written guides have original worked examples and starting references, but actual branch heads, complete link availability, and CI status must be read back after publication. Do not claim 19 published branches merely from those local Markdown files.

> **Source presence does not prove source-derived correctness.** The reviewed subject guides contain link-out references and no apparent copied notebook code. The contributor source `jonkrohn/ML-foundations` was reported MIT-licensed by the wave-one researcher. This reviewer did not independently click and audit all 34 bibliography URLs or conduct a third-party copyright audit. No project license selection was observed.

### Residual risk register

- **Medium / high confidence:** Domain-sensitive symbolic simplification cannot serve as a full assumptions-aware CAS. In particular, polynomial conventions around `0**0`, real-valued logarithm domains, and cancellation should be specified for clients. The parser blocks execution syntax, but this is not a sandbox for arbitrary computation, nor has its resource usage been fuzzed at scale.
- **Medium / high confidence:** Finite-sample adaptive quadrature is vulnerable in principle to adversarial alias frequencies despite improved initial panels. Its own error estimate is heuristic, not a proof for arbitrary input functions.
- **Medium / medium confidence:** Floating-point overflow and underflow may still occur in other geometry and stochastic pathways; only the explicitly tested cases were corrected. Higher-dimensional conditioning and discontinuous metrics require focused stress tests.
- **Low / high confidence:** The subject-guide corpus remains an educational starting point, not an exhaustive proof corpus. A formal system (e.g., lambda calculus) should never be represented as a derivative operator merely because it is named a calculus.

**Evidence classes:** Direct runtime test for named functions; source inspection for interfaces, resource limits, docs and links; no independent CI, no deployed service, no proof assistant, no statistical convergence proof, and no full-source external license audit.
