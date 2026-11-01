# Finite time-scale calculus

**Scope:** `calc.time_scales.FiniteTimeScale` implements the delta (forward) and nabla (backward) operators **only for an explicitly listed, finite, strictly increasing subset of real time points**. It does not implement the full theory on arbitrary closed subsets of the reals, continuous intervals, or accumulation points.

## Mathematical contract

For a time scale \(\mathbb T\), the forward and backward jumps are

\[
\sigma(t)=\inf\{s\in\mathbb T:s>t\},\qquad
\rho(t)=\sup\{s\in\mathbb T:s<t\}.
\]

At the included maximum/minimum, respectively, standard endpoint conventions set \(\sigma(\max\mathbb T)=\max\mathbb T\) and \(\rho(\min\mathbb T)=\min\mathbb T\). The graininess functions are \(\mu(t)=\sigma(t)-t\) and \(\nu(t)=t-\rho(t)\).

For a finite grid \(t_0<\cdots<t_n\), every \(t_i\), \(i<n\), is right-scattered, and every \(t_i\), \(i>0\), is left-scattered. The module computes

\[
f^\Delta(t_i)=\frac{f(t_{i+1})-f(t_i)}{t_{i+1}-t_i}\quad(i<n),\qquad
f^\nabla(t_i)=\frac{f(t_i)-f(t_{i-1})}{t_i-t_{i-1}}\quad(i>0).
\]

The finite delta integral is the **left-sampled** sum; the nabla integral is the **right-sampled** sum:

\[
\int_{t_i}^{t_j}f(t)\,\Delta t=\sum_{k=i}^{j-1}f(t_k)(t_{k+1}-t_k),\quad
\int_{t_i}^{t_j}f(t)\,\nabla t=\sum_{k=i+1}^{j}f(t_k)(t_k-t_{k-1}).
\]

Reversing integration bounds negates the result. Equal bounds yield zero without evaluating the callback. For finite time scales, the maximum has no delta derivative and the minimum has no nabla derivative: their jump conventions do **not** manufacture derivatives there.

## Executable examples

```python
from calc.time_scales import FiniteTimeScale

grid = FiniteTimeScale([0, 1, 2, 3])
assert grid.sigma(1) == 2
assert grid.rho(1) == 0
assert grid.delta_derivative(lambda t: t**2, 2) == 5
assert grid.nabla_derivative(lambda t: t**2, 2) == 3
assert grid.delta_integral(lambda t: t, 0, 3) == 3  # 0+1+2
assert grid.nabla_integral(lambda t: t, 0, 3) == 6  # 1+2+3
```

An irregular grid changes the denominators:

\(\mathbb T=\{0,\tfrac12,2\}\) gives \(D_\Delta t^2(\tfrac12)=(4-\tfrac14)/(2-\tfrac12)=\tfrac52\), rather than an ordinary derivative of \(1\).

On a finite geometric grid \(\{1,2,4,8\}\), \(\sigma(t)=2t\) for nonterminal points, and

\[
f(t)=t^2\quad\Longrightarrow\quad f^\Delta(t)=\frac{(2t)^2-t^2}{2t-t}=3t.
\]

This is the finite right-scattered counterpart of the Jackson q-derivative for \(q=2\), *not* an assertion that all time-scale operators are identical to q-calculus. The discrete fundamental theorem is witnessed by the telescoping sum \(3(1)+6(2)+12(4)=63=8^2-1^2\).

## API and limitations

- `FiniteTimeScale(points)` requires a nonempty ordered finite sequence, with no duplicates or nonfinite values. Coordinates and callback results use binary floating-point values; near-equal points that collapse to one float are rejected.
- `sigma`, `rho`, `mu`, `nu`, `is_right_scattered`, `is_left_scattered` expose jump/graininess structure.
- `delta_derivative` and `nabla_derivative` require callable real-valued functions and eligible points.
- `delta_integral` and `nabla_integral` require both bounds to be members of the grid.
- Finite arithmetic, overflow, endpoint and membership failures raise explicit exceptions. This is not arbitrary-precision, interval-certified or symbolic calculus.
- No right-dense/left-dense derivative limits, continuous integration, infinite time scales, dynamic equations, or mixed continuous/discrete accumulation sets are represented. Those require separate mathematical and computational contracts.

## Primary references

1. Martin Bohner and Allan Peterson, *Dynamic Equations on Time Scales: An Introduction with Applications* (Birkhäuser, 2001), [sample chapters hosted by Peterson](https://www.math.unl.edu/~apeterson1/sample_book.pdf) and [official institutional bibliography](https://scholarsmine.mst.edu/math_stat_facwork/474/).
2. Martin Bohner and Allan Peterson, [First and Second Order Linear Dynamic Equations on Time Scales](https://scholarsmine.mst.edu/math_stat_facwork/325/) (2001).
3. [Dynamics of Time Scales](https://www.math.ksu.edu/research/reu/results2012_files/taylor_monomials_paper.pdf), Kansas State University research paper, jump/graininess and delta/nabla preliminary definitions.

These are external references only. No textbook text, diagrams or code were copied into this implementation.
