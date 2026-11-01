# Dedicated branch: stochastic

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Stochastic calculus

**Definition.** Itô calculus integrates adapted processes with respect to semimartingales, in particular Brownian motion. Unlike ordinary smooth-path calculus, Brownian increments have nontrivial quadratic variation, giving an additional second-derivative term in Itô's formula.

**Worked example.** For standard Brownian motion \(W_t\), Itô's rule applied to \(f(w)=w^2\) yields

\[d(W_t^2)=2W_t\,dW_t+dt.\]

Taking expectations gives \(E[W_t^2]=t\), because the Itô integral is mean-zero under the appropriate integrability conditions. In Stratonovich notation the ordinary-looking chain rule holds, but its stochastic integral differs from the Itô one.

**Failure boundary.** A single Monte Carlo path or seeded pseudorandom simulation cannot prove a stochastic identity, pathwise convergence or distributional convergence. Noise coefficient evaluation depends on whether a numerical method uses Itô or Stratonovich semantics; a silent switch can introduce drift bias.

**Repository relation.** Subject-specific stochastic simulation can illustrate increment scaling and SDE discretization; it is not Malliavin calculus or a rigorous stochastic integral implementation. See [differential calculus](differential.md), [fractional calculus](fractional.md), and the [advanced index](../CATALOG.md).

**Primary starting source:** [MIT OCW Topics in Mathematics with Applications in Finance, Lectures 24–25](https://ocw.mit.edu/courses/18-642-topics-in-mathematics-with-applications-in-finance-fall-2024/).

