# Time-scale calculus

**Definition.** A time scale \(\mathbb T\) is a nonempty closed subset of the real line. Its forward jump \(\sigma(t)=\inf\{s\in\mathbb T:s>t\}\), with the usual endpoint convention, and graininess \(\mu(t)=\sigma(t)-t\) govern the delta derivative. On \(\mathbb T=\mathbb R\), delta differentiation reduces to ordinary differentiation; on \(\mathbb T=\mathbb Z\), it is a forward difference.

**Worked example.** On integer time with \(f(t)=t^2\), the delta derivative is \(f^\Delta(t)=f(t+1)-f(t)=2t+1\). On real time, the same \(f\) has \(f'(t)=2t\). Both are manifestations of one formalism with different graininess.

**Failure boundary.** Endpoint definitions and left/right scattered points matter. A generic closed set can have mixed dense and isolated points; replacing its delta derivative with one fixed finite step fails. Domain conventions must be fixed before a robust program is possible.

**Repository relation.** A catalogue entry and conceptual bridge, **not** a full time-scale solver. See [discrete calculus](discrete.md) and [q-calculus](q-calculus.md).

**Primary starting source:** [Time-scale calculus applied to pharmacokinetics, PMC (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11579191/). This is an application paper, not a complete axiomatic source.

