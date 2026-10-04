# de Bruijn–Newman constant of the Davenport–Heilbronn function: code notes

- **Certified runs** use python-flint 0.8.0 (Arb ball arithmetic).
- **Exploratory runs** use mpmath.
- Outputs stored in `runs/` are named in the table. Rows marked "not stored" record an exploratory observation whose output is not kept in this repository, and no claim depends on them.

## The lower bound and exploration

| file                    | what                                                                       | result                                                    |
| ----------------------- | -------------------------------------------------------------------------- | --------------------------------------------------------- |
| dh\_explore.py           | combination and functional-equation checks; the first off-line zero       | ρ = 0.808517 + 85.699348i (`runs/dh_explore.txt`)         |
| dh\_track\_zero85\_gh.py   | heat-flow track of the first off-line zero (Gauss–Hermite)                 | axis reached for t ∈ (0.176, 0.178) (`runs/dh_track_zero85_gh.txt`) |
| dh\_track\_zero85.py      | the same track by adaptive quadrature                                      | not stored                                                |
| dh\_track\_zero240.py     | track of the zero from 0.86953 + 240.40467i                                | axis reached for t ∈ (0.2300, 0.2310) (`runs/dh_track_zero240.txt`) |
| cert\_dh.py              | the H\_t evaluator; also a segment-ball winding test                        | the winding test failed (enclosures too wide); not used   |
| cert\_dh2.py             | Rouché certificate (Arb)                                                   | used by dh\_rouche\_run.py                                  |
| dh\_rouche\_run.py        | Rouché at t = 0.22 and 0.23                                                | **Λ\_DH > 0.23** (claims #1; `certificates/03-davenport-heilbronn/rouche-t023.json`) |
| dh\_sigma0\_upper.py      | zero-free half-plane (Arb bisection)                                       | σ₀ = 1.3951362, **Λ\_DH ≤ 1.60254** (claims #2)            |
| dh\_head\_torus\_upper.py  | an attempt to improve #2 with an exact head bound                          | no gain; not stored                                       |
| dh\_track\_zero857.py, dh\_track\_zero531\_783.py, dh\_track\_zero2442.py, dh\_refine\_1000\_3000.py | tracks of higher off-line zeros | not stored                                                |

## The upper bound (conditional; papers/03-davenport-heilbronn-lambda/claims.md, #4)

The certificate is `certificates/03-davenport-heilbronn/lambda-058.json` (t0 = 0.5, y0 = 0.4, X = 5e4), with the rectangle runs in `rect_runs_050/`.

| file                    | what                                                                         | result                                           |
| ----------------------- | ---------------------------------------------------------------------------- | ------------------------------------------------ |
| rs\_identity\_dh.py       | Theorem A′ (the RS integral identity) for the DH coefficients, numerically   | relative residual ≤ 3.7e-22 (`runs/rs_identity_dh.txt`) |
| c0\_census.py            | c₀ over the fundamental discriminants with \|D\| ≤ 400 (proof-rs-general.md, Lemma G3) | `runs/c0_census.txt`                   |
| dh\_family\_dim.py        | dimension of the real self-dual space for 3 ≤ q ≤ 19                         | `runs/dh_family_dim.txt`                         |
| dh\_final\_region\_scan.py, dh\_next\_scan.py | float feasibility scans for the final region                | not stored; used only to choose (t0, y0, X)      |
| dh\_upper.py             | part (a): the final region x ≥ X, y ∈ [y0, Ytop], verified arithmetic        | in lambda-058.json                               |
| dh\_rect\_scan.py, dh\_direct\_timing.py, dh\_trap\_check.py | reconnaissance and timing of the evaluators  | not stored                                       |
| dh\_rect.py              | part (b): G = H/B\_t avoids (−∞, 0] on ∂R, 1-D Chebyshev certificates         | `rect_runs_050/`                                 |
| dh\_rect\_summary.py      | assembles (a) + (b) into lambda-058.json                                     | —                                                |

**Why Theorem A′ holds for DH (a proof, not only the numerical check).**

- **What the proof uses.** papers/supporting/riemann-siegel-identity/proof-rs-general.md uses the coefficients only through:
  - periodicity mod q;
  - parity (a₋ₙ = (−1)^𝔞 aₙ);
  - a₀ = 0 and Σ\_{n mod q} aₙ = 0;
  - real values;
  - the Gauss-sum self-duality of Lemma G4, aₙ = τ⁻¹ Σ\_m a\_m e^{2πimn/q} with τ = i^𝔞√q.
- **DH has the first five properties.** It is 5-periodic and odd (𝔞 = 1), with a₀ = 0 and 1 + κ − κ − 1 = 0.
- **Self-duality holds exactly.**
  - Put u = √(10 − 2√5) and v = √(10 + 2√5), so uv = 4√5 and κ = (u − 2)/(√5 − 1).
  - The cases n = 1, 2 reduce to 2u = (√5 − 1)v and 2v = (√5 + 1)u. Squared, both sides equal 40 ∓ 8√5; since both sides are positive, the squaring is reversible.
  - The cases n = 3, 4 follow by oddness.

**Part (b) method.**

- **The reduction.** B\_t(z) = exp(lMt((1 + y − ix)/2, t)) is analytic and nonvanishing for y > −3. If G = H\_t/B\_t avoids (−∞, 0] on ∂R, then G has winding number 0 around R, so H\_t has no zeros in R.
- **The pieces.**
  - Bottom edge: direct for x ≤ 704, RS for x ≥ 700.
  - Top edge: direct for x ≤ 404, RS for x ≥ 400.
  - Left side x = 0: direct.
  - Right side x = X = 5e4: part (a), where |G − 1| < 1, with the error terms evaluated from x = X.
- **RS pieces.**
  - The error e\_C comes from Lemma I with c = N + ½ fixed per N-block and δ localised to fifths of the block, with K\_h from Lemma H (cellwise interval bound; 1.179).
  - The Lemma I tail uses the Gaussian e^{−(π/2q)u²} and needs (4/q)πU/4 > a. That condition fails at c = 10.5, and near x = 640 the bottom-edge RS error exceeds the margin. Hence RS starts at x = 700 (bottom) and x = 400 (top).
  - The error e\_A + e\_B is ε\_max times the block sums.
- **Direct pieces.**
  - The trapezoid rule on the Gaussian ξ-average, with step 0.6 and strip half-width 3.
  - The discretisation error uses the Trefethen–Weideman bound 2M/(e^{2πa/h} − 1), with M ≤ e^{a²/4t}·sup ∫ e^{−w²/4t}|ξ|.
  - The tail uses `cert_dh.DH.tail_bound` (from W − h).
- **ξ evaluation.**
  - Values are taken at 1 − s for Re s < ½, which avoids the Γ poles.
  - Balls near s = 1, where the Hurwitz terms have poles that cancel, are bounded by maximum modulus.
  - Near x = 0 the grid offset is chosen to keep every node away from s = 0 and s = 1.
- **Guards.** A nan margin is treated as a failure (`not (lb > 0)`).
