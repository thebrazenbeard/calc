# Geometry track: finite-difference coordinate calculations

This is an executable introduction, **not** a complete differential-geometry or tensor-calculus system. All computations are local to a single coordinate chart, and derivatives are centered finite differences with heuristic accuracy.

## Contracts

`calc.geometry` takes ordinary Python callables. Coordinates and tangent vectors are sequences of finite numbers. The metric is a callable `metric(point)` returning a symmetric, nondegenerate \(n\times n\) matrix; it may be indefinite. Coordinate chart and metric must be sufficiently smooth near the point.

- `metric_inner(g, x, u, v)` computes \(g_{ij}(x)u^iv^j\).
- `coordinate_jacobian(F, x)` computes \(J^a{}_i=\partial_iF^a\) for a map whose output dimension can differ from its input dimension.
- `pullback_metric(F, g, x)` computes \(J^Tg(F(x))J\).
- `christoffel_symbols(g, x)` computes the Levi-Civita symbols
  \[
  \Gamma^i{}_{jk}=\tfrac12g^{i\ell}(\partial_jg_{\ell k}+\partial_kg_{\ell j}-\partial_\ell g_{jk}).
  \]
  Output indexing is `gamma[upper][lower1][lower2]` (zero-based).
- `exterior_derivative_0(f, x)` returns the components of \(df\).
- `exterior_derivative_1(w, x)` returns antisymmetric two-form components for increasing index pairs, \((d\omega)_{ij}=\partial_i\omega_j-\partial_j\omega_i\), where \(\omega=\sum\omega_i dx^i\).

The Christoffel symbols are **not tensors**; they transform with an inhomogeneous term. The Riemann curvature tensor, integration of forms over chains, manifold atlases, coordinate transitions, and symbolic simplification are not implemented here.

## Independent witnesses

Polar coordinates with \(r>0\) give \(ds^2=dr^2+r^2d\theta^2\), and hence

\[
\Gamma^r{}_{\theta\theta}=-r,\quad
\Gamma^\theta{}_{r\theta}=\Gamma^\theta{}_{\theta r}=1/r.
\]

The Euclidean metric has zero connection coefficients in Cartesian coordinates, but the same flat plane has nonzero coefficients in polar coordinates. The polar chart is singular at \(r=0\), correctly rejected by `christoffel_symbols`.

For the one-form \(\omega=-y\,dx+x\,dy\), the exact exterior derivative is \(d\omega=2\,dx\wedge dy\). For a smooth scalar \(f\), \(d(df)=0\) analytically; floating-point numerical approximations need tolerances.

## Numerical caveats

- Relative finite-difference step defaults to `1e-5` times the coordinate scale. It is tunable, not an error guarantee.
- Extremely ill-conditioned metrics or coordinates near singularities may be rejected or lose precision.
- Metric validation checks symmetry and nondegeneracy at the evaluation point, not global geometry.
- No numerical method can infer smoothness or correct topology from a bare callable.

## Sources

These are references, not copied source material:

1. [MIT OpenCourseWare 18.950: Differential Geometry](https://ocw.mit.edu/courses/18-950-differential-geometry-fall-2008/), especially Chapter 2 notes on connection and Christoffel symbols.
2. [University of California, Irvine: Differential Geometry lecture notes](https://lu.math.uci.edu/pdfs/lecture_notes/Math-240BC/dgbook.pdf), section on differential forms and the exterior derivative.
3. [MIT 8.962: General Relativity lecture notes](https://web.mit.edu/8.962/www/lecnotes/8_962TA-lec-all.pdf), coordinate transformations and Christoffel-symbol calculations.
