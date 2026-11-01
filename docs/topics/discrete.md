# Discrete calculus

**Definition.** Forward difference with spacing \(h\ne0\) is \(\Delta_h f(x)=f(x+h)-f(x)\). Division by \(h\) yields a difference quotient that may approach a derivative as \(h\to0\), but at fixed spacing it is a distinct operator. Summation plays the role of antidifference.

**Worked example.** \(\Delta_1(x^2)=(x+1)^2-x^2=2x+1\), not \(2x\). A second difference is \(\Delta_1^2(x^2)=2\). Summing the first difference telescopes: \(\sum_{k=0}^{n-1}(2k+1)=n^2\).

**Failure boundary.** Ordinary product rules have shifted terms: \(\Delta(fg)(x)=f(x)\Delta g(x)+g(x+h)\Delta f(x)\). Reusing the unshifted differential Leibniz rule gives wrong identities. Differences also depend on the step and boundary conventions.

**Repository relation.** `calc.advanced` contains selected normalized finite-difference examples; finite-difference approximations to smooth derivatives also occur in `calc.differential`. See [q-calculus](q-calculus.md), [time-scale calculus](time-scales.md) and [differential calculus](differential.md).

**Primary starting source:** [SymPy official finite-difference documentation](https://docs.sympy.org/latest/explanation/special_topics/finite_diff_derivatives.html).

