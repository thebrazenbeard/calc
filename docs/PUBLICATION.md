# Publication and verification receipts

This record distinguishes **GitHub writes/readbacks** from local tests and CI. It is not a claim that all mathematical paradigms have been implemented, or that cloud CI passed.

## Checkpoint history (2026-11-01)

| Subject | Verified exact SHA | Evidence |
| --- | --- | --- |
| Initial bootstrap | `9f2e8611c988bbea32dab8ac5eaca716bc784f83` | GitHub create-file response and commit fetch |
| Wave 1 foundation on `main` | `81aebf196b2b742a3da7242e3357ea97c0baf4c8` | GitHub ref update and remote README/algebra file readback; 71/71 local tests |
| Wave 2 on `main` | `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9` | GitHub ref update; 135/135 local tests and README sample |
| Wave 3 on `main` | `43d8aad9a0dfaf2bcec7fe8f0d0e25a45d5f6a87` | GitHub ref update and remote PDE/lambda/source receipt file readback; 213/213 local tests; Python wheel built and imports passed |

Wave 4 review, including five module-focused independent hostile checks and a research-gap audit, raised the local suite to **261/261 passing tests** after repairs. The new tests cover overflow and precision failures in complex calculus, discrete transforms, PDEs, optimization and formal lambda-term reduction. Wave 4's exact remote head must be read back separately after publication.

The standard-library tests were executed with Python 3.12 in a local isolated workspace. At the Wave 2 check, GitHub's combined status API returned an empty `statuses` list. That does **not** establish a passing cloud CI run. Branch workflows and test execution on GitHub were not independently read back.

## Curated subject branches

Each branch started from the exact Wave 2 `main` commit above and received its own `TOPIC.md` with an original subject guide. GitHub content readback checked that file's slug on all 19 branches; no failures were observed. These branches are **research/teaching entry points**, not independently complete solvers.

| Branch | Exact `TOPIC.md` commit SHA |
| --- | --- |
| `topic/abstract-algebra` | `cd095fb1d41e730e61f608ee8fdc7b878e4cc2dc` |
| `topic/boolean-algebra` | `138aab7b6328325a395d2e231838a3cd2200e030` |
| `topic/calculus-of-variations` | `9d4dcbf529e35fe9b7ad252c283bce9632e80ae7` |
| `topic/complex-analysis` | `78567ff7b00d37af925e414cbd20e858a0123f7a` |
| `topic/differential-geometry` | `f5bbacbd529dc16ee91dc120dacf519ab4803a6f` |
| `topic/differential` | `e5ffb0e08cd7a15286b591b317139f7dadc63f91` |
| `topic/discrete` | `aa98bf72c596181e309ac649af68d8542a389625` |
| `topic/exterior` | `830c006b62a52abe6ab98e9c2587f7170ef7691c` |
| `topic/fractional` | `53c95cf502692c9ab1425b8545599b3da62ed932` |
| `topic/integral` | `08d5f55d99791b9900175ef11f955238e7673f95` |
| `topic/lambda-calculus` | `28c566ab91f37e492d4334b7f8b958b831a3f465` |
| `topic/lie-algebra` | `e5cc223cb9d8f7d1cc10b4243a89edcf9920c2de` |
| `topic/linear-algebra` | `1909546ebff085eb19e839abcb6800641a5914be` |
| `topic/multivariable` | `d33d81a7d1d80704f345f229f8f7f5c3e863a542` |
| `topic/q-calculus` | `af54bd79565f1fc0c870cdba3c157dedae810039` |
| `topic/stochastic` | `e1e671b5ebd926ceeade9686dcbd1c4cbc4f7ba1` |
| `topic/time-scales` | `86a1e9a7874e28a892d32eb5c6a5fc3414dd0e15` |
| `topic/tropical-algebra` | `d60ca0b9854fca3e6cb828187b948f7f833712f0` |
| `topic/vector` | `e6402b12abe743f1cf01b2c54369345d3889ef76` |

## Known boundaries

- The catalogue is open-ended: an indexed name is not an implementation or a proof.
- The topic branches are intentionally pinned to a known Wave 2 foundation; later `main` fixes do not automatically propagate into their history.
- No repository license has been chosen. The 34 linked external research sources have independent terms.
- The recorded tests do not formally verify the numerical or stochastic algorithms for all possible inputs.
