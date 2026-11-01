# Independent hostile review: complex calculus

**Scope:** `src/calc/complex_calc.py` at the Wave 3 local integration surface on 2026-11-01. This is independent mathematical testing by the `complex_hostile` native specialist, not a formal proof or a cloud CI verification. The reviewer did not alter production code.

## Requirement/risk matrix

| Requirement or risk | Independent witness | Evidence class |
| --- | --- | --- |
| Directional complex derivative | `d(exp z)/dz = exp z` across four distinct directions | Analytic invariant + executed unittest |
| Nonholomorphic diagnostic | `d(conj z)` differs by direction; CR residual is `(2,0)` | Analytic invariant + executed unittest |
| Path parameterization | `∫(2+z) dz = 5/2` along `[0,1]` under two affine schedules | Exact antiderivative + executed unittest |
| Contour orientation | `∮ dz/z = 2πi`; reverse path negates integral | Cauchy's integral formula + executed unittest |
| Winding domain | Winding requires a closed loop; common translation of path and center cannot change this | Geometric/topological invariant + executed negative unittest |
| Finite differences | Perturbations must move the represented point | IEEE-754 representability + executed negative unittest |
| Simpson parameter step | Distinct endpoints must produce nonzero, advancing quadrature nodes | IEEE-754 representability + executed negative unittest |

## Reproduced defects (pre-fix)

1. **Medium, high confidence — coordinate-scale-dependent open-path acceptance.** `winding_number(lambda t: 1e12+1000+0.01*t, lambda t: 0.01, center=1e12, start=0, stop=1)` returned `0` even though start and end differ by roughly `0.01`. The same segment without the large common translation raised `ValueError`. Cause: closure threshold used `max(1,|path_start|,|path_end|,|center|)`; translation changes the threshold, violating path-closure invariance. A robust closure check should use local geometry, ideally a conservative absolute tolerance plus a relative threshold based on path displacement/size around the center, with representability checks.
   The same scale bug also **rejects a valid shifted unit circle**: `winding_number(lambda t: 1e12+exp(1j*t), lambda t: 1j*exp(1j*t), center=1e12)` raises `ValueError` claiming the path meets its center, while the translated circle around zero correctly has winding `1`. The center-distance threshold must also use local radius and precision, not absolute world coordinates.
2. **Medium, high confidence — unrepresentable complex differentiation step.** `complex_derivative(lambda z: z*z, 1, step=5e-324)` returned `0j` where the analytic derivative is `2`. `1±h` round to the same float. Reject a non-progressive point/step rather than report a false derivative.
3. **Low/medium, high confidence — unrepresentable contour mesh spacing.** `contour_integral(lambda z: 1, lambda t: complex(t), lambda t: 1, 0, 5e-324, panels=2)` returned zero because the derived panel width underflowed to zero; the exact integral equals the minimum positive subnormal. Reject non-progressive quadrature grids rather than hide a precision loss.

## Limitations that no finite sample can close

- Local CR checks cannot establish global holomorphicity or identify an unsampled singularity.
- Simpson quadrature can alias arbitrarily oscillatory integrands. More panels reduce some errors but do not give a universal error certificate.
- Numerical winding tests are heuristic. Homotopy/topological statements require independent regularity and path nonintersection assumptions.
- Caller-supplied Python callables execute with the caller's authority; the module is not a sandbox for untrusted functions.

**Remediation boundary:** Parent integrator owns production fixes and final integrated verification. Test failures should remain visible until corresponding mathematical invariants hold, not be weakened to achieve a green suite.

## Post-repair verification

The parent integrator subsequently amended `complex_calc.py` to reject unrepresentable spatial and quadrature steps, use a coordinate-translation-stable closure threshold, and check center proximity using local radius plus floating-point ULP resolution. The tests above were **not relaxed**.

Fresh local command:

```text
PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_complex.py -v
Ran 8 tests — OK (exit 0), Python 3.12
```

This qualifies only the tested examples and negative input classes. Continuous holomorphicity and exact winding remain unproven by the numerical algorithm.
