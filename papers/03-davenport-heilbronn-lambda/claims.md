# Claims and evidence — the de Bruijn–Newman constant of the Davenport–Heilbronn function

| # | Claim (precise) | Rigor | Digits | Novelty | Certificate | Lean statement | Witness | Notes |
| --- | --------------- | ----- | ------ | ------- | ----------- | -------------- | ------- | ----- |
| 1 | For the Davenport–Heilbronn function f (coefficients (1, κ, −κ, −1, 0) mod 5), H\_{0.23}(z) = (4π·0.23)^{−1/2}∫e^{−w²/0.92}ξ((1 + w + iz)/2)dw has exactly one zero within 10⁻⁵ of a point with Im z ≈ −0.0349, hence a non-real zero. Therefore Λ\_DH > 0.23 | certificate |  | no de Bruijn–Newman or heat-flow study of this function found (searches September–October 2026) | certificates/03-davenport-heilbronn/rouche-t023.json |  |  | Arb Rouché certificate (cert\_dh2.py). |
| 2 | For the same f, every zero of H₀ lies in \|Im z\| ≤ 2σ₀ − 1 with σ₀ = 1.3951362 (the zero-free half-plane from 1 − Σ\_{n≥2}\|a\_n\|n^{−σ} > 0, plus ξ(s) = ξ(1 − s)). Hence Λ\_DH ≤ 1.60254 by de Bruijn's strip theorem | certificate |  |  | certificates/03-davenport-heilbronn/upper-bound.json |  |  | Arb bisection (dh\_sigma0\_upper.py). |
| 3 | (Withdrawn.) A bound Λ\_DH ≤ 0.68 at t₀ = 0.6 is not claimed: t₀ = 0.6 lies outside the range t ≤ ½ of Prop 6.1^χ. Superseded by #4 | withdrawn |  |  |  |  |  | Kept only so that the numbering of #4 is stable. |
| 4 | For the same f, H\_{0.5} has no zeros with \|Im z\| ≥ 0.4. Hence Λ\_DH ≤ 0.5 + 0.4²/2 = 0.58 by de Bruijn's strip theorem, **conditional on the open items O1–O4 of papers/supporting/riemann-siegel-identity/derivation.md** | conjectural (computation: certificate) |  |  | certificates/03-davenport-heilbronn/lambda-058.json |  |  | Parameters (t₀, y₀, X) = (0.5, 0.4, 5·10⁴); runs in code/03-davenport-heilbronn/rect\_runs\_050/. |

## Reasoning

### #1

- **Why a non-real zero bounds Λ.** H\_t is real entire of order 1, since ξ(s) = ξ(1 − s) and the coefficients are real. For t ≥ Λ it lies in the Laguerre–Pólya class, and the backward heat operator preserves that class without any positivity of Φ. So a non-real zero at t₁ gives Λ > t₁.
- **Why the zero exists.**
  - Rouché compares H with its linearisation p on a circle of radius r. |h₀| and d come from Arb point values.
  - |H − p| is bounded by Cauchy estimates from a rigorous modulus majorant on a disk of radius 0.05 (M ≤ 248 at t = 0.23 and M ≤ 215 at t = 0.22).
  - The margins are 9.9·10⁻⁶ < 1.71·10⁻⁵ at t = 0.23, and 1.37·10⁻⁴ < 2.82·10⁻⁴ at t = 0.22.
- **Numerics.** The zero from 0.86953 + 240.40467i reaches the axis for t between 0.2300 and 0.2310 (`code/03-davenport-heilbronn/runs/dh_track_zero240.txt`), so the bound is near-sharp for this zero.
- **Ceiling.** f is not in the Selberg class; this says nothing about ζ.

### #2

- **The zero-free half-plane.** For σ > 1, |f| ≥ 1 − Σ\_{n≥2}|a\_n|n^{−σ}. The sum is decreasing in σ and certified < 1 at σ₀. By symmetry the zeros lie in 1 − σ₀ ≤ σ ≤ σ₀, i.e. |Im z| ≤ 2σ₀ − 1.
- **The bound.** De Bruijn: if every zero of H₀ has |Im z| ≤ Δ, then H\_t is real-rooted for t ≥ Δ²/2.

### #4

- **The argument** (main.md §4.2): the strip of #2 at t₀ = 0.5, then a final region for x ≥ 5·10⁴, a winding certificate on the boundary of [0, 5·10⁴] × [0.4, Ytop], and symmetry.
- **The certificate.** The final region has minimum lower bound 0.0604 over 482 blocks, against an error bound of 0.0272 (evaluated from x = X), and a tail lower bound ≥ 0.95.
- **The rectangle.** It has 50,359 boundary pieces with worst margin 0.1025:
  - direct evaluation on the bottom edge for x ≤ 704, on the top edge for x ≤ 404, and on the side x = 0;
  - the RS approximation beyond.
  - `dh_rect_summary.py` checks the coverage, every ok flag and the common (t₀, y₀, X).
- **Why conditional.** See main.md §4.3: Prop 6.1^χ (O1), Lemma T for T′ < 10⁶ (O2), Lemma M's constant and its range y ≤ 1 (O3), and Lemma I for q = 5 (O4).
- **What limits it.** The float feasibility scan (dh\_next\_scan.py) indicates that smaller (t₀, y₀) need a larger final-region X, and the rectangle cost grows with X.
- **Ceiling.** It concerns the Davenport–Heilbronn function only. Unconditionally, Λ\_DH > 0.23; conditionally, Λ\_DH ∈ (0.23, 0.58].
