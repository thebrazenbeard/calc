# Integral calculus

**Definition.** The Riemann integral of a bounded function on a closed interval is the common limit of its tagged sums as partition mesh tends to zero, when that limit exists independently of tags. Lebesgue integration uses a different construction and applies under other hypotheses.

**Worked example.** For \(f(x)=x^2\) on \([0,1]\), the antiderivative \(x^3/3\) and the fundamental theorem of calculus give

\[\int_0^1 x^2\,dx=[x^3/3]_0^1=1/3.\]

The independent exact value is the test oracle for a numerical quadrature routine.

**Failure boundary.** A function with a nonintegrable singularity cannot be made integrable by asking the numerical sampler to skip the singular point. For example, the improper integral \(int_0^1 x^{-1}\,dx\) diverges. The fundamental theorem requires suitable continuity/absolute-continuity assumptions; numerical sums need error analysis.

**Repository relation.** `calc.numerical.trapezoid`, `simpson`, and `adaptive_simpson` approximate definite integrals; they do not establish general integrability and their estimated tolerances are not certificates for arbitrary functions. See [differential calculus](differential.md), [fractional calculus](fractional.md), and the [measure-theory catalogue](../CATALOG.md).

**Primary starting source:** [MIT OCW 18.01SC](https://ocw.mit.edu/courses/18-01sc-single-variable-calculus-fall-2010/).

