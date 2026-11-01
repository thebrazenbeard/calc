# Adversarial review — bounded 1D PDE solvers

**Reviewer:** native subagent `/root/pde_hostile`; **target:** local
`src/calc/pde.py` and `docs/PDE.md` on 2026-11-01, after the third
implementation wave. This is independent model review plus executable
Python 3.12 checks, **not** a formal proof or a GitHub Actions result.
The parent integrator owns Git commits and publication. Exact publication
commit must be recorded separately after GitHub readback.

## Material findings and repair evidence

### 1. Heat FTCS falsely overflows at stable CFL equality

**Severity:** medium; **confidence:** high (direct runtime reproduction).

Original expression evaluated `values[j-1] - 2*values[j] + values[j+1]`
before applying its factor. With `initial=[0, 1e308, 0]`, `dx=1`,
`diffusivity=1`, `dt=.5`, `steps=1`, the scheme has `r=1/2` and the exact
interior output is `0`. The original code instead raised `ArithmeticError`
because its intermediate `2e308` overflowed.

**Repair observed:** production now computes the stable convex combination
`(1-2r) u_j + r u_{j-1} + r u_{j+1}` using compensated summation. The
corresponding new regression passed after the parent integrator's fix.

### 2. Wave scheme falsely overflows at Courant equality

**Severity:** medium; **confidence:** high (direct runtime reproduction).

With `initial=[0,1e308,0]`, zero velocity, `wave_speed=dx=dt=1`,
`steps=2`, the Courant-one discrete scheme should give the finite center
sequence `(1e308, 0, -1e308)`. The original first-step stencil raised
`ArithmeticError` from its `-2*center` intermediate.

**Repair observed:** the parent regrouped the first and subsequent
recurrences to avoid the cancellable multiplication and introduced finite
result validation. The independent regression passed after this change.

These failures were established by *direct pre-fix runtime probes*.
The new unittest cases were run against the corrected code; no claim is made
that they were executed as unittest failures before repair.

## Independent mathematical test matrix

| Required invariant or risk | Independent oracle | Result |
| --- | --- | --- |
| Heat CFL `r=1/2` finite pulse | Convex-neighbor formula | PASS after repair |
| Wave Courant `ν=1` finite pulse | Exact three-level discrete recurrence | PASS after repair |
| Heat boundedness | Maximum principle for convex stencil | PASS |
| Heat spatial convergence | Analytic `exp(-π²t) sin(πx)` mode; grid refinement 11→21→41 | PASS |
| Dimensional scaling | Same `κdt/dx²` yields same heat history | PASS |
| Wave dimensional scaling | Same `c dt/dx` and zero initial velocity yields same wave history | PASS |
| Wave initial-velocity term | Analytic sine-mode solution `sin(πt)/π` | PASS |
| Fixed endpoints | Boundary invariance under all steps | PASS |
| Immutable input/history | Tuple output and no input mutation | PASS |
| Invalid negative speed and zero time step | Rejection | PASS |

Executed from `calc/`:

```sh
PYTHONPATH=src python -m unittest discover -s tests -p 'test_hostile_pde.py' -v
```

Result after repairs: **10 tests, 0 failures, exit code 0**.

## Unresolved limits

1. The tests support the declared **bounded 1D explicit schemes** only. No
   general PDE solver, convergence proof for nonsmooth inputs, certified
   global error bound, or high-dimensional/stiff-system capability is claimed.
2. A CFL limit bounds the linear discretization's stability; it does not
   guarantee useful accuracy on a coarse grid. Near singular scales and
   floating-point extremes remain numerical rather than formally verified.
3. Documentation links to MIT OCW were present, but their remote availability
   and current copyright terms were not independently rechecked in this pass.
4. Tests ran in the local Python runtime. GitHub publication, CI status and
   branch qualification are separate evidence and were not audited here.

**Verdict:** The two directly reproduced defects were repaired and the
ten independent focused contracts pass. **Specification: PASS within the
tested bounded solver contract. Mathematical generality and external release
qualification: PARTIAL / not established.**
