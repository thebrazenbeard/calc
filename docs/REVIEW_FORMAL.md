# Independent hostile review: lambda and symbolic calculi

Review target: `src/calc/lambda_calc.py` and `src/calc/symbolic.py` in the local Wave 3 integration workspace, 2026-11-01. Reviewer: native specialist `/root/lambda_hostile` (internal, independent of implementation). Evidence: direct source inspection and fresh Python 3.12 `unittest` execution; **not** formal verification, network security testing, or a separate provider/model review.

## Material finding

**F-01 — Medium severity, high confidence: a manually constructed deep lambda AST escapes the advertised reduction resource error.**

`normalize(term, max_nodes=10)` first calls recursive `_count_nodes`. An `App` nested 1,200 levels deep raises Python's `RecursionError` before comparing the node count to `max_nodes`. The public `normalize` contract advertises an explicit node budget and `ReductionLimitExceeded` when a budget is crossed.

Minimal deterministic reproduction (no external inputs or side effects):

```python
from calc.lambda_calc import App, Var, normalize
t = Var("x")
for _ in range(1200):
    t = App(Var("f"), t)
normalize(t, max_nodes=10)  # observed: RecursionError
```

**Consequence:** a caller feeding programmatically assembled ASTs cannot reliably handle the documented bounded-evaluation exception; deep trees can exhaust Python's recursion stack instead of receiving a deterministic size refusal. This is a bounded robustness gap, *not* evidence that the source parser executes Python or that inputs parsed through the documented grammar bypass its nesting limit.

Smallest remediation: budget-aware, iterative node counting/preflight before recursive normalization, or a precise fail-closed recursion conversion to `ReductionLimitExceeded` followed by an internal depth guarantee. Keep the failing regression in `tests/test_hostile_formal.py` until repaired.

## Independent mathematical/security falsifiers

| Requirement or risk | Independent witness | Observed initial result |
| --- | --- | --- |
| Capture avoidance | Substitute free `y` into nested `λy.λz.x y`; check α-equivalence and free variables | Pass |
| Alpha-equivalence | Compare λ-term binding depths, not variable spellings | Pass |
| Nontrivial term computation | Church numeral addition of two and three normalizes to five | Pass |
| Normal-order laziness | Outer K-like term discards an Ω argument without evaluating it | Pass |
| Reversible presentation | Parse formatted nested abstraction/application terms | Pass |
| Reduction budget on manually deep AST | Build 1,200 App nodes; request `max_nodes=10` | **Fail (RecursionError)** |
| Symbolic domain preservation | `0*log(x)`, `0*(1/x)`, `(1/x)^0` at singular points | Pass |
| Symbolic calculus | Product and chain derivative of `exp(x)*sin(x²)` at five analytic points | Pass |
| Parser safety | Reject imports, attributes, comprehensions, lambda calls, subscripting | Pass |
| Real logarithm domain | Derivative of `log(x)^2` at `e`, reject negative argument | Pass |

## Explicit limits of this review

- All parser probes were local. No operating-system command or file access was requested through either parser.
- The lambda parser has explicit source/token/depth limits, reducing exposure from text input. Programmatically constructed trees use a different path and motivated F-01.
- Tests do not prove all capture-avoiding substitution cases or the absence of every parser denial-of-service input.
- Symbolic simplification deliberately avoids general algebraic cancellation; a derivative may algebraically extend a function's apparent domain and must not be taken as domain proof.

## Initial verification receipt

`PYTHONPATH=src python -m unittest discover -s tests -p 'test_hostile_formal.py' -v`

Initial result: **9 tests; 8 pass, 1 error** (`test_manual_deep_ast_respects_node_budget_before_python_recursion`). A further independently derived Church numeral test was added afterward (now 10 tests; 9 pass, 1 error until repaired). The failure is a deterministic implementation behavior, not a missing dependency. Follow-up verification must be recorded after any fix on the final integrated source; this document does not pre-approve such a fix.

## Local repair and fresh readback

The integrator subsequently changed `lambda_calc.py` to use an iterative, budget-aware node-count preflight, reject excessive tree depth (256) independently of a generous node budget, and convert residual recursion exhaustion to the domain-specific reduction exception. The reviewer reran the same regression with both `max_nodes=10` and `max_nodes=100000` on a 1,200-level AST. Both now raise `ReductionLimitExceeded` as expected.

Fresh verification: `PYTHONPATH=src python -m unittest discover -s tests -p 'test_hostile_formal.py' -q` — **10/10 passed, exit 0**, Python 3.12. F-01 is **resolved in the local shared workspace**; no claim is made here about remote publication or CI against a final GitHub SHA.
