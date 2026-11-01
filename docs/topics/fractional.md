# Fractional calculus

**Definition.** Fractional calculus studies derivatives/integrals of noninteger order, but no single operator covers all contexts. For \(\alpha>0\), the left Riemann–Liouville fractional integral from zero is

\[(I^\alpha_{0+}f)(t)=\frac1{\Gamma(\alpha)}\int_0^t(t-s)^{\alpha-1}f(s)\,ds.\]

**Worked example.** For \(f(t)=1\), a change of variables gives \(I^\alpha 1=t^\alpha/\Gamma(\alpha+1)\). At \(\alpha=1/2\), \(t=1\), this is \(2/\sqrt\pi\). More generally, \(I^\alpha t^\beta=\Gamma(\beta+1)t^{\alpha+\beta}/\Gamma(\alpha+\beta+1)\) for the stated integrable domain.

**Failure boundary.** Riemann–Liouville and Caputo derivatives treat initial constants differently: for \(0<\alpha<1\), the Caputo derivative of a constant vanishes, whereas the Riemann–Liouville derivative generally does not. A model or program must name the operator, lower terminal and function space.

**Repository relation.** `calc.advanced` covers a narrow analytic power-function fractional integral, not arbitrary fractional differential equations. See [integral calculus](integral.md) and [stochastic calculus](stochastic.md).

**Primary starting source:** [Review of fractional derivative definitions, Journal of Computational Physics (2019)](https://doi.org/10.1016/j.jcp.2019.03.008).

