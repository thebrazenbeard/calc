# Differential calculus

**Definition.** A real function \(f\) is differentiable at an interior point \(a\) if the finite limit

\[f'(a)=\lim_{h\to0}\frac{f(a+h)-f(a)}{h}\]

exists. The derivative is local linear response, not merely a slope sampled at two distant points. Differentiability implies continuity, but the converse is false.

**Worked example.** For \(f(x)=x^3\), expand \((a+h)^3-a^3=3a^2h+3ah^2+h^3\). Dividing by \(h\) and taking the limit gives \(f'(a)=3a^2\); at \(a=2\), the answer is \(12\). The original expansion is a mathematical witness independent of a numerical differencing algorithm.

**Failure boundary.** The absolute-value function has left derivative \(-1\) and right derivative \(+1\) at zero, so its ordinary derivative there does not exist. Forward finite differences are not an exact symbolic derivative and can suffer cancellation when \(h\) is too small.

**Repository relation.** `calc.differential.derivative` approximates a derivative. Exact polynomial differentiation lives separately in `calc.algebra`; analytic limits and error bounds remain theoretical work. Continue to [integral calculus](integral.md), [multivariable calculus](multivariable.md), and [automatic/symbolic method inventory](../CATALOG.md).

**Primary starting source:** [MIT OCW 18.01SC: Single Variable Calculus](https://ocw.mit.edu/courses/18-01sc-single-variable-calculus-fall-2010/). Referenced, not reproduced.

