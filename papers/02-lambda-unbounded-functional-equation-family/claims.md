# Claims and evidence — Λ across the q = 7 Theorem A′ family

| # | Claim (precise) | Rigor | Digits | Novelty | Certificate | Lean statement | Witness | Notes |
| --- | --------------- | ----- | ------ | ------- | ----------- | -------------- | ------- | ----- |
| 1 | For θ = 1/10 and a\_θ = cos θ·(·/7) + sin θ·b (b the exact DFT-self-dual odd partner of main.md §1), H\_{0.30} for f\_θ has exactly one zero within 2·10⁻⁵ of about 102.085 − 0.1925i (exact centre in the certificate file), hence a non-real zero. Therefore Λ(f\_{1/10}) > 0.30 | certificate |  | no de Bruijn–Newman study of such families found (searches September–October 2026); the family itself goes back to Heilbronn 1938 | certificates/02-functional-equation-family/rouche-theta1\_10-t030.json |  |  | Arb Rouché certificate (x7\_cert.py). |
| 2 | In the q = 7 odd Theorem A′ family, Λ is upper semicontinuous at the L-function: limsup\_{θ→0} Λ(f\_θ) ≤ Λ(L(s, χ₋₇)). Hence, under GRH(χ₋₇), Λ(f\_θ) → 0 | conjectural |  |  |  |  |  | proof-semicontinuity.md. Steps 0 and 2–4 are complete. Lemma 1′ (uniform far-field dominance) is sketched. It rests on the effective approximation of papers/supporting/riemann-siegel-identity/derivation.md, whose Prop 6.1^χ is not written out (derivation.md, open items O1 and O5). |
| 3 | For θ = 704/1000 (next to θ\* = arctan(−1/b(1)), where a\_θ(1) = 0) in the same family, H\_{24.5} for f\_θ has exactly one zero within 10⁻⁶ of about 4.4719 − 0.9283i, hence a non-real zero. Therefore Λ(f\_{704/1000}) > 24.5 | certificate |  |  | certificates/02-functional-equation-family/rouche-theta704\_1000-t2450.json |  |  | Arb Rouché certificate (x7\_cert.py, high-precision mode, Stirling tail bound and ξ-reflection, W = 120, R = 0.6). |
| 4 | The rate as θ → 0 for the zero near height 51.14. In the same family, H\_t for f\_θ has a non-real zero (Rouché, one zero within 10⁻⁵ of the centre in each certificate file) at (θ, t) = (1/20, 0.12), (3/100, 0.06) and (1/50, 0.02). Therefore Λ(f\_{1/20}) > 0.12, Λ(f\_{3/100}) > 0.06 and Λ(f\_{1/50}) > 0.02 | certificate |  |  | certificates/02-functional-equation-family/rouche-theta1\_20-t012.json, rouche-theta3\_100-t006.json, rouche-theta1\_50-t002.json |  |  | Arb Rouché certificates (x7\_cert.py, R = 0.05, r = 10⁻⁵, δ = 10⁻⁶). Each t lies within one tracking step (0.02) of this zero's numerical collision time. |
| 5 | Growth toward θ\*: for θ = 70406/100000 in the same family (a\_θ(1) ≈ 1.85·10⁻⁶), H\_{30} for f\_θ has exactly one zero within 10⁻⁶ of about 4.4275 − 1.7507i, hence a non-real zero. Therefore Λ(f\_{70406/100000}) > 30 | certificate |  |  | certificates/02-functional-equation-family/rouche-theta70406\_100000-t3000.json |  |  | Arb Rouché certificate (x7\_cert.py, high-precision mode, Stirling tail and ξ-reflection, W = 150, R = 1). |
| 6 | Growth toward θ\*, third point: for θ = 70406119530951995496/10²⁰ in the same family (a\_θ(1) ≈ 1.0·10⁻¹⁰), H\_{44} for f\_θ has exactly one zero within 10⁻⁶ of about 6.1011 − 1.0381i, hence a non-real zero. Therefore Λ(f\_θ) > 44 | certificate |  |  | certificates/02-functional-equation-family/rouche-theta70406119530951995496\_100000000000000000000-t4400.json |  |  | Arb Rouché certificate (x7\_cert.py, high-precision mode, Stirling tail and ξ-reflection, W = 220, R = 0.9, one Arb Newton step). |
| 7 | Theorem U: in the q = 7 odd DFT-self-dual family, Λ(f\_θ) → +∞ as θ → θ\* (a\_θ(1) → 0, a\_θ(2) ≠ 0). So sup\_θ Λ(f\_θ) = +∞, although every f\_θ has the functional equation of L(s, χ₋₇) | conjectural |  | no de Bruijn–Newman study of such families found (searches September–October 2026) |  |  |  | Written proof, proof-unbounded.md (a Laplace lemma plus Rouché against the two-term model). It was reviewed adversarially by a second, independent agent of the research system, and the gaps found were fixed. Not formalised. |

## Reasoning

### #1

- **Why a non-real zero bounds Λ.**
  - a\_θ is real, 7-periodic, odd and exactly DFT-self-dual with eigenvalue i√7 (main.md §1). So ξ(s) = ξ(1 − s), with root number 1 (Theorem A′, proof-rs-general.md, Remark 5).
  - Hence H\_t is real entire and even of order 1. For t ≥ Λ it lies in the Laguerre–Pólya class, and the backward heat flow preserves that class.
  - So a non-real zero at t₁ gives Λ > t₁.
- **The zero.**
  - Rouché compares H with its linearisation on |z − c| = 2·10⁻⁵. |h₀| ≤ 1.3·10⁻¹⁴ and |d| ≥ 3.138 come from Arb point values.
  - The modulus majorant is M ≤ 198 on radius 0.05, and it gives |H − p| ≤ 3.17·10⁻⁵ < 6.28·10⁻⁵ ≤ |p|.
  - The coefficients are Arb balls containing the exact values, so the certificate covers the exact f\_{1/10}.
- **Numerics.** The heat-flow track (`runs/x7_track_theta0.1.txt`) reaches the axis between t = 0.30 and 0.32, so the bound is near-sharp for this zero.
- **Ceiling.** f\_θ is not in the Selberg class.

### #2

- **The proof** (proof-semicontinuity.md):
  - (0) the heat flow is linear in the coefficients;
  - (1) uniform far-field dominance leaves no zeros with |x| ≥ X and y₀ ≤ |y| ≤ Y (Lemma 1′, sketched);
  - (2) for t₁ > Λ₀ the zeros of H\_{t₁}^0 are real and simple (Csordas–Smith–Varga 1994, Cor. 1, plus a self-contained backward-heat Hermite argument);
  - (3) Hurwitz on a compact box keeps them real for small θ;
  - (4) de Bruijn's strip lemma gives Λ(f\_θ) ≤ t₁ + y₀²/2.
- **Why the grade is conjectural.** Lemma 1′ is sketched, not written out. A self-contained qualitative large-x asymptotic would suffice. Nothing is formalised.
- **Consequences.**
  - Under GRH(χ₋₇), Λ₀ = 0 (Λ₀ ≥ 0 by Dobner), so Λ(f\_θ) → 0.
  - A uniform positive lower bound for Λ(f\_θ), θ ≠ 0, would disprove GRH(χ₋₇).

### #3

- **Why θ\* matters.** At θ\* the leading coefficient a\_θ(1) vanishes, so the first-term dominance that keeps zeros near the line fails.
- **The deep zero.** Near θ\*, f\_θ ≈ a(1) + a(2)2^{−s}, with a(1) = 9.454·10⁻⁵ at θ = 704/1000. So f\_θ has zeros near the two-term balance σ₀ = log₂|a(2)/a(1)| = 12.863, repeating in t with period 2π/log 2. `deep_zeros.py` refines the lowest to (12.8674, 4.5492), i.e. Im z = −24.735 (`runs/deep_zeros_0.704.txt`).
- **The certificate.** High-precision mode at t = 24.5 gives |h₀| ≤ 2.6·10⁻¹¹, |d| ≥ 0.513, M ≤ 608 on R = 0.6, and |H − p| ≤ 1.7·10⁻⁹ < 5.1·10⁻⁷ ≤ |p| on r = 10⁻⁶. Since r < |Im c|, the zero is non-real.
- **Why the radius is small.** The modulus majorant ∫Gaussian·sup|ξ| ignores the cancellation that makes H small, so Rouché needs r < |d|R²/M.
- **Two features of the machinery,** both rigorous:
  - (i) **A Stirling tail bound (TAIL = 2).** For |w| ≥ W, |ξ| ≤ norm·C\_a·ζ(2)(7/π)^x Γ(x), with x = (|w| + B)/4. The log-integrand is then concave once x > 5t/24, and the code asserts g′(W) < 0. At t = 24.5 this bound exceeds the numerically integrated tail by a factor of about 1.4·10⁵ (`runs/tail_check_t2450.txt`).
  - (ii) **Reflection.** ξ is evaluated at 1 − s when Re s < ½. This is exact for the self-dual a, and the ball coefficients contain the exact a.
- **Ceiling.** Non-multiplicative functions; nothing about ζ.

### #4

- **Why each certificate bounds Λ.** As in #1.
- **The zero.** The off-line zero near height 51.14 at θ = 0.1, 0.05, 0.03, 0.02 has σ − ½ = 0.428, 0.271, 0.181, 0.115 (`runs/x7_track_theta*.txt`).
- **The collision times.** The heat-flow tracks reach the axis in t ∈ (0.30, 0.32), (0.12, 0.14), (0.06, 0.08) and (0.02, 0.04). The certificate times 0.30, 0.12, 0.06 and 0.02 are the last tracking step before the axis.
- **Interpretation (numeric).** This fits an on-line close pair of zeros of L(s, χ₋₇) near t ≈ 51.14 that collides at θ\_c ≈ 0.0132 and leaves the line (the pair-collision law, main.md §4.1). It is consistent with the semicontinuity of #2.
- **Ceiling.** These are lower bounds from one zero at low height.

### #5 and #6

- **The zeros.** The two-term balance gives σ₀ = 18.541 and 32.713. `deep_zeros.py` refines the lowest zeros to (18.5413, 4.5341) and (32.7135, 4.5324), i.e. Im z = −36.083 and −64.427 (`runs/deep_zeros_*.txt`).
- **Choosing θ.** θ\* = arctan(−1/b(1)) = 0.704061195… was computed in Arb (`theta_star.py`). The rational θ of #6 gives a\_θ(1) = 1.00000·10⁻¹⁰ (Arb enclosure; `runs/theta_star.txt`).
- **The certificates.**
  - At t = 30: |h₀| ≤ 4.0·10⁻⁹, |d| ≥ 0.934, M ≤ 5285 on R = 1, and |H − p| ≤ 5.3·10⁻⁹ < 9.3·10⁻⁷ ≤ |p| on r = 10⁻⁶.
  - At t = 44, after one Arb Newton step: |h₀| ≤ 3.8·10⁻⁸, |d| ≥ 5.658, M ≤ 7594 on R = 0.9, and |H − p| ≤ 9.4·10⁻⁹ < 5.6·10⁻⁶ ≤ |p| on r = 10⁻⁶.
- **The rate (heuristic).**
  - The heat flow moves the two-term balance toward the line at the initial rate Re φ′(σ₀) − ½ log 2 ≈ ½ log σ₀ (main.md §5.2).
  - This suggests Λ(f\_θ) ≍ σ₀/log σ₀ ≍ log(1/|a\_θ(1)|)/log log(1/|a\_θ(1)|). The rate is a conjecture.
- **Ceiling.** Certified points and a heuristic rate. Unboundedness is proved separately (Theorem U, #7).

### #7

- **The idea.** Near θ\*, f\_θ ≈ a(1) + a(2)2^{−s}, whose zeros sit at σ₀ = log₂|a(2)/a(1)| → ∞, at bounded height.
- **Step 1 (the Laplace lemma).**
  - For fixed t, the heat flow acts on the n = 1, 2 terms G(s)n^{−s} as multiplication by exp((t/4)(φ′(s) − log n)²)(1 + o(1)), uniformly near the deep zero. This holds because φ″ = O(1/σ₀) and φ‴ = O(1/σ₀²) there.
  - The n ≥ 3 terms, the reflected region and the far tails are negligible: there one needs only upper bounds, dominated by (2/3)^{σ₀} and e^{−w²/4t}.
- **Step 2 (the model).** So H\_t, normalised, is uniformly close on a small circle to the model 1 + sgn·K·2^{−u}. The model has a simple zero at Re s = σ₀ − (t/2)(Re φ′ − ½ log 2) ≈ σ₀ − (t/4) log σ₀.
- **Step 3 (Rouché).** Rouché puts a zero of H\_t there. Its real part tends to ∞, so it is non-real, and hence Λ(f\_θ) > t, by de Bruijn's monotonicity of real-rootedness in t.
- **Numerical support.** Lemma L's drift prediction matches the heat-flow tracks to 0.1% at the three deep zeros (`runs/drift_check.txt`).
- **What would change this.** A flaw in the uniform Stirling estimates at bounded height, which are standard.
- **Ceiling.** Non-multiplicative functions; nothing about ζ. It shows that the functional equation alone puts no upper bound on Λ.
