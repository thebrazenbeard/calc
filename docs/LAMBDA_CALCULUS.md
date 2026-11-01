# Lambda calculus: variables, abstraction, and application

**Domain:** formal computation and mathematical logic. The λ-calculus has no real-valued derivative, numerical integral, or stochastic path implied by its name. It is indexed beside other mathematical *calculi* because it is a system of symbolic transformations.

## A minimal mathematical contract

Untyped terms follow the grammar

\[
M ::= x \mid (\lambda x.M) \mid (M\ N).
\]

An abstraction `λx.M` binds free occurrences of `x` in `M`; application associates left. β-reduction contracts `(λx.M) N` to `M[x := N]`, where substitution is **capture-avoiding**. α-equivalence disregards bound variable names: `λx.x` and `λy.y` denote the same term up to binder renaming.

The implementation in `src/calc/lambda_calc.py` provides immutable `Var`, `Abs`, `App` nodes, a restricted text `parse` function, `free_vars`, `substitute`, `alpha_equivalent`, `beta_step`, `normalize`, and `format_term`.

## Executable examples

```python
from calc.lambda_calc import parse, normalize, format_term, free_vars

identity = parse(r"(\x.x) y")
assert format_term(normalize(identity)) == "y"

constant = parse(r"((\x.\y.x) a) b")
assert format_term(normalize(constant)) == "a"

# A free y must remain free, rather than getting captured by λy.
captured_if_naive = parse(r"(\x.\y.x) y")
result = normalize(captured_if_naive)
assert free_vars(result) == frozenset({"y"})
```

Normal order reduces the *leftmost outermost* redex first. For example, `(λx.a) Ω` yields `a` without reducing its divergent argument `Ω = (λx.x x) (λx.x x)`. By contrast, normalizing `Ω` itself never reaches a normal form, so a finite execution must stop at the configured reduction-step limit.

## Syntax and resource limits

- Text accepts ASCII `\x.M` or Unicode `λx.M`; variable names match `[A-Za-z_][A-Za-z_0-9]*`, and applications are written by juxtaposition (`f x`). Parenthesize a lambda used as an argument: `f (λx.x)`.
- It **does not parse or execute Python**: arbitrary calls, attribute syntax, strings, imports, operators, macros, and programs are rejected.
- The source parser caps length, tokens, and nesting depth. `normalize` limits beta steps and AST node count, raising `ReductionLimitExceeded` rather than claiming termination.
- All reductions are **untyped**. Type checking, η-reduction, evaluation effects, lazy runtime implementation, optimization and observational equivalence are not provided.
- A finite step/node budget is a useful teaching guardrail, **not** a general proof of polynomial time or protection against every resource-exhaustion scenario. Terms can grow rapidly; do not expose this prototype as an unrestricted multi-tenant endpoint.

## Independent mathematical witnesses

1. `(λx.x) y → y` (identity).
2. `((λx.λy.x) a) b → a` (constant/K combinator).
3. `(λx.λy.x) y → λz.y` up to α-equivalence (naive substitution would capture the free `y`).
4. `(λx.a) ((λz.z z) (λz.z z)) → a` in one outermost reduction, despite the divergent argument.
5. `(λx.x x) (λx.x x)` never normalizes under β-reduction (the Ω term).

These examples are exercised in `tests/test_lambda_calc.py`, including input rejection, α-equivalence and shadowing. Code test results are separate from formal proof of any general theorem.

## Primary reference material

- Stanford Encyclopedia of Philosophy, [*The Lambda Calculus*](https://plato.stanford.edu/entries/lambda-calculus/), for syntax, β-conversion, substitutions, and reduction strategies.
- Cornell University CS 6110, [lecture on λ-calculus and capture-avoiding substitution](https://www.cs.cornell.edu/courses/cs6110/2023sp/lectures/lec02.pdf).
- [Stanford-hosted lambda calculus tutorial](https://theory.stanford.edu/~blynn/lambda/), including normal-order behavior and α-conversion.

These are links and references only; no upstream code or instructional prose was copied into this implementation.
