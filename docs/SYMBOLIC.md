# Symbolic calculus: a bounded demonstration

`calc.symbolic` implements a deliberately small, immutable expression tree and differentiation rules. It does **not** implement a computer algebra system, equation solver, proof assistant, or general symbolic integrator.

## Supported expressions

| Expression | Differentiation rule |
| --- | --- |
| Constant `c`, variable `x` | `d(c)/dx = 0`, `d(x)/dx = 1` |
| `f + g` | `f' + g'` |
| `f * g` | `f' g + f g'` |
| `f / g` | `(f'g - fg')/g²`, where `g ≠ 0` |
| `f ** n` | `n f^(n−1) f'` for bounded integer `n` |
| `sin(f)`, `cos(f)` | `cos(f)f'`, `−sin(f)f'` |
| `exp(f)`, `log(f)` | `exp(f)f'`, `f'/f` for `f > 0` in real arithmetic |

```python
from calc.symbolic import Var, sin, log, parse_expression

x = Var("x")
f = sin(x*x) + log(x+2)
derivative = f.diff("x")
print(derivative.eval({"x": 1.0}))

parsed = parse_expression("sin(x**2) + 3*x")
print(parsed.diff("x").eval({"x": 2.0}))
```

Arithmetic on exact integer and `fractions.Fraction` constants remains exact during local constant folding; transcendental evaluation and variable substitution use finite Python floats. Expressions with singularities, non-real logarithms, or non-finite evaluation fail with exceptions. Negative integer powers are supported, but arbitrary variable powers (such as `x**y`) are not.

## Parser security boundary

The optional `parse_expression` translates a narrow, allowlisted subset of Python expression syntax using `ast.parse` *without compiling or evaluating source as Python*. Inputs have a 2,048-character limit, a 256-node default limit, a 40-level default nesting limit, and integer exponents are limited to magnitude 64. These are resource guards, not a formal guarantee against all denial-of-service cases.

Only numeric literals, variables, addition/subtraction/multiplication/division, integer powers, parentheses, and single-argument calls to `sin`, `cos`, `exp`, or `log` are accepted. Function attributes, method calls, file operations, imports, comprehensions, strings, arrays, keywords, and other Python syntax are rejected. **Do not add general Python `eval`, `exec`, or arbitrary callable invocation to this parser.**

## Domain and simplification caveats

Simplification is local; it eliminates additive zero and multiplicative one and folds constant expressions. It does not assume `x/x == 1`, cancel zero times a potentially undefined expression, or prove more general algebraic identities. `0^0` on a literal or variable uses the explicit computational convention `1`, but complex bases such as `(1/x)^0` retain their input domain restrictions.

The formal derivative may have a larger apparent evaluation domain than the original expression. For example, differentiating `0*log(x)` can yield a formula evaluable outside `x > 0`, but this does **not** establish that the original function was differentiable there. Determine the original function's domain first.

## References and verification

- [MIT OpenCourseWare: Single Variable Calculus](https://ocw.mit.edu/courses/18-01sc-single-variable-calculus-fall-2010/) — derivative rules and examples.
- [Python `ast` documentation](https://docs.python.org/3/library/ast.html) — parser syntax-tree representation and cautions against untrusted complexity.
- [SymPy](https://docs.sympy.org/latest/tutorials/intro-tutorial/calculus.html) — a full-featured symbolic system; cited for comparison only, no code copied.

Tests exercise polynomial derivatives, chain/product/quotient rules, trig/exponential/log derivatives, independent variables, singularities, simplification identities, parser rejection, and resource limits. A passing finite test suite is not a proof of arbitrary symbolic equivalence.
