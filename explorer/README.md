# Calculus field notes: browser explorer

A self-contained interactive visualization of selected introductory calculus
ideas. Open `explorer/index.html` directly in a modern browser, or serve the
repository root with:

```bash
python -m http.server 8000
```

Then visit `http://localhost:8000/explorer/`.

The interface works offline. It uses ordinary HTML, CSS, SVG, and JavaScript;
it has no Python runtime, CDN, third-party scripts, analytics, dynamic code
execution, or network dependency. No mathematical input is compiled or
evaluated: the only selectable functions are explicitly defined in `app.js`.

## What can be explored

- Four explicit differentiable functions: `x²`, `sin(x)`, `x³ − 3x`, and
  `1 / (1 + x²)`.
- A movable tangent point (pointer drag, keyboard-accessible slider, or arrow
  keys while the graph is focused).
- Adjustable horizontal viewing radius from 2 to 8.
- Signed area under the function between two user-chosen in-view bounds.
- Separate analytic reference values and numerical estimates, never presented
  as formal error certificates.

The tangent is drawn from the selected function's known analytic derivative.
The derivative estimate is the symmetric difference

\[\frac{f(x+h)-f(x-h)}{2h},\qquad h=10^{-5}\max(1,|x|).\]

The integral estimate uses composite Simpson quadrature with 600 subintervals.
The reference integral is computed from a known antiderivative; the area fill
visually distinguishes positive and negative contributions. Because all
arithmetic is binary64, smooth-function examples are approximations rather than
certified numerical bounds.

The explorer illustrates scalar functions only. It does not parse symbolic
expressions, prove limits, solve arbitrary integrals, or execute Python modules.
The Python algorithms and evidence catalogue remain in `src/calc/`.

## Verification

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_explorer_assets.py -v
```

These standard-library tests inspect the distributed asset contract and basic
offline/security invariants. They are **not** a substitute for running the
interface in a browser and testing keyboard, pointer, focus, and narrow-screen
behavior. A deployer can host this static directory under GitHub Pages after
enabling Pages in the repository settings; this project does not silently
publish it.

## Accessibility and design

- Semantic labels, real HTML controls, live numerical outputs, keyboard and
  SVG arrow-key navigation, and visible focus indicators.
- Mobile layouts down to 360 CSS px and reduced-motion support.
- Graph has a dynamically updated accessible description and text readouts;
  the graph itself is supplemental rather than the sole delivery of values.

The local assets and interface text are original to this repository. No
third-party code or artwork is bundled in the explorer.
