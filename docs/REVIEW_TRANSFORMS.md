# Independent hostile review: discrete transforms

**Scope:** `src/calc/transforms.py`, `tests/test_transforms.py`, `docs/TRANSFORMS.md`; additional independent contracts in `tests/test_hostile_transforms.py`. Review performed locally on Python 3.12 on 2026-11-01. This is a **MODEL_REVIEW + TEST**, not formal verification or proof for arbitrary floating-point inputs. No source files were modified by this reviewer.

## Material findings

### F1 — zero z-transform coefficient triggers false overflow (confirmed; fixed)

- **Severity:** medium. **Confidence:** high (direct runtime reproduction before parent patch).
- **Reproduction before correction:** `z_transform_finite((1,0,0),1e-300)` raised `ValueError` because the algorithm evaluated `z**(-2)` before multiplying by a zero coefficient.
- **Independent mathematical result:** the finite sum is exactly 1 for every nonzero z. Rejecting it is wrong even within binary floating-point arithmetic.
- **Correction:** the parent integrator added a zero-coefficient guard. `test_trailing_zeros_must_not_force_spurious_overflow` was executed afterward and passed.

### F2 — nonzero weighted term overflows before scaling (confirmed; fixed)

- **Severity:** medium. **Confidence:** high (direct runtime reproduction and independently calculated finite answer).
- **Reproduction:** `z_transform_finite((1,0,1e-300),1e-200)` raised `ValueError` on `z**(-2)`.
- **Independent mathematical result:** `1 + 10^-300 * (10^-200)^-2 ≈ 10^100` is representable as a float. It is only the *unscaled intermediate* `z**(-2)` that overflows, so this is a numerical algorithm defect rather than an unavoidable range limit for the requested value.
- **Correction:** the parent integrator added a guarded fallback that repeatedly divides the coefficient by `z` if the direct inverse-power path overflows. Intermediate values are checked for finiteness; the path does not claim arbitrary precision.
- **Regression:** `test_small_coefficients_can_cancel_exponent_range` initially failed; it now passes for positive-real and imaginary `z`, including the expected negative complex-phase sign.

## Independent contract coverage

| Requirement | Analytic or independent oracle | State at review |
| --- | --- | --- |
| DFT/IDFT normalization | N=1, odd-length roundtrip, frequency-bin shift and modulation theorems | Passed |
| Complex-sample energy | Parseval identity under forward-unnormalized convention | Passed |
| Cancellation stability | DC sum `1e16 + 1 - 1e16` equals 1 | Passed |
| Linear/circular convolution | Impulse identity, commutativity, associativity, zero-padding equivalence | Passed |
| Finite z-transform | Negative real and imaginary z, zero coefficients, representable weighted terms | Passed after F1/F2 repairs |
| Input/boundary rejection | Nonfinite, empty, Boolean/string, true overflow | Passed |

The tests use explicit mathematical identities, small hand-derived examples, and a independent scale argument rather than reproducing the implementation's loops.

## Residual limits (not implementation defects under stated scope)

- Direct DFT is `O(N²)`, not FFT. Input length is capped at 4096 but there is no performance service-level objective or benchmark in this review.
- Arbitrary finite precision cannot guarantee correctly rounded Fourier coefficients for ill-conditioned cancellation, underflow, or large magnitudes. Tests do not provide a uniform error bound.
- Finite z-transform is a finite sum, not an infinite transform or a region-of-convergence determination. None of the continuous Fourier/Laplace/Mellin theory is implemented here.
- The transform module is dependency-free but the repository currently has no chosen project license; the external NumPy/SciPy documentation is linked as a convention reference and no third-party implementation is vendored.
- GitHub CI, exact post-publication repository head, and cloud runtime execution are outside this local review and remain unverified by this reviewer.

**Verdict:** locally tested transform contracts **PASS** after two reproduced defects were repaired. Overall release qualification remains **PARTIAL** pending final exact-commit integration tests, GitHub CI readback, and further numerical conditioning coverage. Fresh focused command: `PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_transforms.py -v` — **11/11 passed**, exit 0.
