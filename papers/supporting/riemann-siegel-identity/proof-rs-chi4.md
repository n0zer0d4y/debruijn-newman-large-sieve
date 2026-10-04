# The Riemann–Siegel integral identity for L(s, χ₋₄): a complete proof

This is a written proof. The numeric spot checks (§9) are numeric-experiment grade; nothing here has been formalised.

## Notation

- ω := e^{iπ/4}, so ω² = i, ω̄ = e^{−iπ/4}, iω = −ω̄ and iω̄ = ω.
- χ = χ₋₄, and L(s) = Σ\_{n≥1} χ(n) n^{−s}.
- **Branch.** log is the principal branch on ℂ ∖ (−∞, 0], and w^{a} := exp(a log w). Note log w̄ = conj(log w) off the cut.
- **The lines.** For real c, let ℓ\_c := {c + uω : u ∈ ℝ} and ℓ̃\_c := {c + uω̄ : u ∈ ℝ}.
  - ∫\_{c↙} is the integral over ℓ\_c traversed from u = +∞ to u = −∞.
  - ∫\_{c↘} is the integral over ℓ̃\_c traversed from u = −∞ to u = +∞.
- **The functions.**
  - Ψ(w) := ω e^{iπw²/4} / (4 cos(πw/2)), and Ψ\*(w) := conj Ψ(w̄).
  - F(s) := ∫\_{c↙} w^{−s} Ψ(w) dw.
  - F̃(s) := ∫\_{c↘} w^{s−1} Ψ\*(w) dw.
  - X(s) := (π/4)^{s−1/2} Γ((2 − s)/2) / Γ((s + 1)/2).
- **Two Mellin transforms.**
  - A(u) := ∫₀^∞ ξ^{u−1} e^{−iξ²/π} / (2 cosh ξ) dξ, for Re u > 0.
  - B(s) := ∫₀^∞ x^{s−1} e^{ix²/π} / (2 cosh x) dx, for Re s > 0.

  Both converge absolutely, because |e^{±ix²/π}| = 1 and 1/cosh decays exponentially. Note B(s) = conj A(s̄).

**Theorem.** For every c ∈ (0, 1) and every s ∈ ℂ,

L(s, χ₋₄) = F(s) − X(s) F̃(s)

as meromorphic functions of s. In particular F̃(2k) = 0 for every integer k ≥ 1, which cancels the poles of X.

**Strategy.**

- The Mordell identity gives a closed form for one entire function Φ (Lemma 3).
- On ℓ\_c the quantity Re(ωw) = c/√2 is constant and positive. So w^{−s} has a Laplace representation along the whole line, and F becomes a Mellin transform of Φ (Prop. 4).
- Φ is its own cosine transform (Lemma 6), which gives a second Mellin form of F (Prop. 7).
- Reflection (Lemma 8) and a Γ-identity (Lemma 10) then close the argument in the strip 0 < σ < 1, and analytic continuation finishes it.

The geometric-series-plus-Hankel route suggested earlier is not used. On ℓ\_c the series Σ χ(n) e^{±iπnw/2} does not converge uniformly at the real crossing point w = c, and the Laplace trick avoids that problem.

## 0. Preliminaries

**(P1) Decay along the lines.**

- For w = a + uω with a, u real:
  - w² = a² + 2auω + iu², so Im(w²) = √2·au + u² and |e^{iπw²/4}| = exp(−(π/4)(u² + √2 au)).
  - |cos(πw/2)| ≥ |sinh(π Im w/2)| = sinh(π|u|/(2√2)).
  - If a is not an odd integer, |cos(πw/2)| is bounded below on the line.
- Hence |Ψ(w)| ≤ C\_a exp(−πu²/4 + O(|u|)) on ℓ\_a. This is locally uniform for a in any compact set avoiding the odd integers, and the O(|u|) is uniform too. The same holds for Ψ\* on ℓ̃\_a, since |Ψ\*(w)| = |Ψ(w̄)|.

**(P2) Branch and entireness.** Let c > 0 and w ∈ ℓ\_c.

- Im w = u/√2 vanishes only at w = c > 0, so ℓ\_c misses (−∞, 0].
- arg w ∈ (−3π/4, π/4) and |w| ≥ c/√2, the distance from 0 to ℓ\_c.
- Moreover ωw = cω + iu, so **Re(ωw) = c/√2 > 0**, arg(ωw) = arg w + π/4 ∈ (−π/2, π/2), and log(ωw) = log w + iπ/4.
- Hence |w^{−s}| = |w|^{−σ} e^{t·arg w} ≤ |w|^{−σ} e^{3π|t|/4}.
- With (P1), the integral defining F converges absolutely and locally uniformly in s, so **F is entire**.
- The same argument on ℓ̃\_c (arg w ∈ (−π/4, 3π/4)) shows **F̃ is entire**.

**(P3) Independence of c.** Take 0 < c₁ < c₂ < 1. The closed strip between ℓ\_{c₁} and ℓ\_{c₂} meets ℝ only in [c₁, c₂] ⊂ (0, 1). So it contains no pole of Ψ (the poles are the odd integers) and no point of the cut. Apply Cauchy's theorem on the parallelogram cut off at u = ±R: the two short sides have bounded length and the integrand is O(e^{−πR²/8}) on them, by (P1). Letting R → ∞ shows the integrals over ℓ\_{c₁} and ℓ\_{c₂} agree. The same holds for F̃ and for Φ below.

**(P4) Moving past the pole at w = 1.**

- The residue. Near w = 1, cos(πw/2) = −(π/2)(w − 1) + O((w − 1)³) and e^{iπ/4} = ω. So Ψ(w) = ω²/(−2π(w − 1)) + O(1) = 1/(2πi(w − 1)) + O(1), i.e. Res\_{w=1} Ψ = 1/(2πi).
- The claim. Let E be entire with |E(w)| ≤ α e^{β|w|}, and let 0 < c < 1 < c′ < 3. Then

  ∫\_{c↙} Ψ E dw − ∫\_{c′↙} Ψ E dw = 2πi · Res\_{w=1}(ΨE) = E(1).

- Proof. ℓ\_c is traversed in direction ω̄·(−1) = e^{5πi/4}, i.e. down and to the left, and its left-hand normal points towards larger Re w. So the strip between the two lines lies to the left of ℓ\_c as traversed, and the positively oriented boundary is ℓ\_c (as traversed) followed by ℓ\_{c′} reversed. The only pole in the strip is w = 1. The truncation argument of (P3) applies, because the Gaussian of (P1) dominates e^{β|w|}.

## 1. The Mordell identity and the function Φ

**Lemma 1 (Mordell identity).** Let Q(w) := (ω/2) e^{iπ(w² + 2w)/4}. Then

Ψ(w) + Ψ(w + 2) = Q(w).

*Proof.* There are three facts: e^{iπ(w+2)²/4} = e^{iπw²/4} e^{iπw} e^{iπ} = −e^{iπw} e^{iπw²/4}; cos(π(w + 2)/2) = −cos(πw/2); and 1 + e^{iπw} = 2 e^{iπw/2} cos(πw/2). Hence

Ψ(w) + Ψ(w + 2) = ω e^{iπw²/4}(1 + e^{iπw})/(4 cos(πw/2)) = (ω/2) e^{iπw²/4} e^{iπw/2} = Q(w). ∎

**Lemma 2 (Gaussian line integral).** For every b ∈ ℂ and every real c,

∫\_{c↙} e^{iπ(w+b)²/4} dw = −2ω.

*Proof.*

- Put β := (c + b)ω̄. Then c + b + uω = ω(u + β), so (c + b + uω)² = i(u + β)² and the integrand is e^{−π(u+β)²/4}.
- Moving the real Gaussian's contour by the complex shift β is justified by Cauchy's theorem and Gaussian decay, so ∫\_ℝ e^{−π(u+β)²/4} du = 2.
- With dw = ω du and the orientation u: +∞ → −∞, the integral is −2ω. ∎

**Lemma 3 (closed form of Φ).** For z ∈ ℂ let

Φ(z) := ∫\_{c↙} Ψ(w) e^{iπzw/2} dw, 0 ≤ c < 1.

Then Φ is entire in z and independent of c ∈ [0, 1), and

2 cos(πz/2) · Φ(z) = 1 − ω e^{−iπz²/4}.

*Proof.*

- **Entire and independent of c.** On ℓ\_c, |e^{iπzw/2}| ≤ e^{(π|z|/2)|w|}, which (P1) absorbs. So Φ is entire and, by (P3), independent of c. The value c = 0 is allowed here: there is no w^{−s}, so no cut, and there is no pole in [0, c].
- **Integrate the Mordell identity.** Multiply Lemma 1 by e^{iπzw/2} and integrate over ℓ\_c. All three integrals converge by (P1).
  - The shifted term. The substitution v = w + 2 maps ℓ\_c onto ℓ\_{c+2} with the same orientation, so

    ∫\_{c↙} Ψ(w + 2) e^{iπzw/2} dw = e^{−iπz} ∫\_{(c+2)↙} Ψ(v) e^{iπzv/2} dv = e^{−iπz} (Φ(z) − e^{iπz/2}).

    The last step is (P4) with E(v) = e^{iπzv/2}, since c + 2 ∈ (2, 3) for c ∈ (0, 1).

  - The Gaussian term. Q(w) e^{iπzw/2} = (ω/2) e^{−iπ(1+z)²/4} e^{iπ(w+1+z)²/4}, so by Lemma 2 its integral is (ω/2) e^{−iπ(1+z)²/4}(−2ω) = −i e^{−iπ(1+z)²/4}.

- **Solve for Φ.** The identity is Φ(z)(1 + e^{−iπz}) − e^{−iπz/2} = −i e^{−iπ(1+z)²/4}. Multiply by e^{iπz/2} and note iπz/2 − iπ(1 + z)²/4 = −iπ/4 − iπz²/4 and i e^{−iπ/4} = ω:

  2 cos(πz/2) Φ(z) = 1 − i e^{−iπ/4} e^{−iπz²/4} = 1 − ω e^{−iπz²/4}. ∎

*Remarks.*

- At odd z = n both sides vanish, because e^{−iπn²/4} = e^{−iπ/4} = ω̄. This is consistent with Φ being entire.
- Φ is even in z.
- Setting z = 2ix/π with x real gives

  **h(x) := Φ(2ix/π) = (1 − ω e^{ix²/π}) / (2 cosh x)**,

  because −iπz²/4 = ix²/π and cos(ix) = cosh x.

## 2. F as a Mellin transform

**Proposition 4.** For σ = Re s > 0,

Γ(s) F(s) = ∫₀^∞ x^{s−1} h(x) dx.

*Proof.*

(i) **Laplace representation.** For Re a > 0 and σ > 0, a^{−s} = Γ(s)^{−1} ∫₀^∞ y^{s−1} e^{−ay} dy with principal branches: this holds for a > 0 and extends by analytic continuation in a. By (P2), Re(ωw) = c/√2 > 0 on ℓ\_c and w^{−s} = e^{iπs/4}(ωw)^{−s}. Hence

w^{−s} = e^{iπs/4} Γ(s)^{−1} ∫₀^∞ y^{s−1} e^{−ωwy} dy (w ∈ ℓ\_c).

(ii) **Fubini.** |Ψ(w) y^{s−1} e^{−ωwy}| = |Ψ(w)| y^{σ−1} e^{−cy/√2} is integrable over ℓ\_c × (0, ∞), by (P1) and σ > 0. Swapping the integrals, with −ωwy = iπζw/2 for ζ(y) = 2iωy/π,

Γ(s) F(s) = e^{iπs/4} ∫₀^∞ y^{s−1} Φ(2iωy/π) dy.

(iii) **Rotating the ray.**

- For z = 2iωy/π: z² = −4iy²/π², so −iπz²/4 = −y²/π; and πz/2 = iωy, so cos(πz/2) = cosh(ωy). By Lemma 3,

  Φ(2iωy/π) = (1 − ω e^{−y²/π})/(2 cosh(ωy)).

  This is analytic in y on the closed sector S = {−π/4 ≤ arg y ≤ 0} minus 0.

- For y = Re^{iφ} ∈ S:
  - arg(ωy) ∈ [0, π/4], so Re(ωy) ≥ R/√2 and |cosh(ωy)| ≥ sinh(R/√2);
  - |e^{−y²/π}| = e^{−R² cos 2φ/π} ≤ 1;
  - |y^{s−1}| ≤ R^{σ−1} e^{π|t|/4}.
- So the integral over the arc |y| = R, −π/4 ≤ φ ≤ 0, tends to 0 as R → ∞. The arc |y| = r contributes O(r^σ) → 0 as r → 0. Cauchy's theorem therefore moves the ray arg y = 0 to arg y = −π/4.
- Write y = ω̄x with x > 0. Then (ω̄x)^{s−1} dy = ω̄^{s} x^{s−1} dx, with the principal branch at arg = −π/4, and e^{iπs/4} ω̄^{s} = 1. Also 2iω(ω̄x)/π = 2ix/π. So

  Γ(s) F(s) = ∫₀^∞ x^{s−1} Φ(2ix/π) dx = ∫₀^∞ x^{s−1} h(x) dx. ∎

**Corollary 5.** For σ > 0,

Γ(s) F(s) = Γ(s) L(s) − ω B(s).

*Proof.* Split h(x) = 1/(2 cosh x) − ω e^{ix²/π}/(2 cosh x). It remains to show Γ(s)L(s) = ∫₀^∞ x^{s−1}/(2 cosh x) dx for σ > 0.

- For σ > 1: 1/(2 cosh x) = e^{−x}/(1 + e^{−2x}) = Σ\_{k≥0} (−1)^k e^{−(2k+1)x} = Σ\_n χ(n) e^{−nx}. Integrating termwise is legitimate, because Σ\_n n^{−σ} Γ(σ) < ∞.
- Both sides are analytic for σ > 0, since the integral converges locally uniformly there, so the identity extends to σ > 0. ∎

## 3. h is a cosine transform

**Lemma 6.** For real x,

h(x) = (2ω̄/π) ∫₀^∞ e^{−iξ²/π} cos(2xξ/π) / (2 cosh ξ) dξ.

*Proof.*

- Let m(x) := ∫\_ℝ e^{−iξ²/π} e^{2ixξ/π}/(2 cosh ξ) dξ, which converges absolutely. Taking even parts, m(x) = 2 ∫₀^∞ e^{−iξ²/π} cos(2xξ/π)/(2 cosh ξ) dξ, so it suffices to prove m(x) = πω h(x).
- **Rotate the ξ-line onto ω̄ℝ.**
  - The two sectors swept are {−π/4 ≤ arg ξ ≤ 0} and {3π/4 ≤ arg ξ ≤ π}. They contain no pole of 1/cosh; the poles lie on iℝ.
  - On |ξ| = R in the first sector, write ξ = Re^{−iθ} with 0 ≤ θ ≤ π/4. The modulus of the integrand is at most

    exp(−R² sin 2θ/π + 2|x| R sin θ/π) · 2e^{−R cos θ}/(1 − e^{−2R cos θ}).

  - Using sin 2θ ≥ 4θ/π and sin θ ≤ θ, the exponent is at most θR(2|x|/π − 4R/π²) − R/√2 ≤ −R/√2 once R ≥ π|x|/2. So the arc contribution tends to 0; the second sector is symmetric.
  - With ξ = ω̄r, and −iξ²/π = −r²/π:

    m(x) = ω̄ ∫\_ℝ e^{−r²/π} e^{2ixω̄r/π}/(2 cosh(ω̄r)) dr.

- **Compute Φ on ℓ₀.** On ℓ₀ (c = 0 is allowed by Lemma 3), w = uω with u: +∞ → −∞.
  - Ψ(uω) = ω e^{−πu²/4}/(4 cosh(πuω̄/2)), because cos(πuω/2) = cosh(iπuω/2) = cosh(−πuω̄/2).
  - Hence Φ(z) = −(ω²/2) ∫\_ℝ e^{−πu²/4} e^{iπzuω/2}/(2 cosh(πuω̄/2)) du.
  - With r = πu/2:

    Φ(z) = −(i/π) ∫\_ℝ e^{−r²/π} e^{izωr}/(2 cosh(ω̄r)) dr.

- **Match the two.** At z = −2ix/π we get izω = (2x/π)ω = 2ixω̄/π, since iω̄ = ω. So m(x) = ω̄·(iπ)Φ(−2ix/π) = iπω̄ Φ(2ix/π), using that Φ is even. Since iω̄ = ω, this is πω h(x). ∎

## 4. A second Mellin form of F

**Proposition 7.** For 0 < σ < 1,

F(s) = (2ω̄/π) (π/2)^{s} cos(πs/2) · A(1 − s).

*Proof.*

- By Prop. 4 and dominated convergence (h is bounded near 0 and O(e^{−x}) at ∞), Γ(s)F(s) = lim\_{ε↓0} ∫₀^∞ x^{s−1} e^{−εx} h(x) dx.
- **Insert Lemma 6.** For ε > 0 the double integral is absolutely convergent: bound it by x^{σ−1} e^{−εx}(2/π)/(2 cosh ξ). So Fubini gives

  ∫₀^∞ x^{s−1} e^{−εx} h(x) dx = (2ω̄/π) ∫₀^∞ (e^{−iξ²/π}/(2 cosh ξ)) I\_ε(2ξ/π) dξ,

  where, for a > 0, I\_ε(a) := ∫₀^∞ x^{s−1} e^{−εx} cos(ax) dx = (Γ(s)/2)[(ε − ia)^{−s} + (ε + ia)^{−s}]. This follows from the Laplace formula of Prop. 4(i), since Re(ε ± ia) = ε > 0.

- **The limit ε → 0.** Since σ > 0 and |ε ± ia| ≥ a, we have |(ε ± ia)^{−s}| ≤ a^{−σ} e^{π|t|/2}. Also (∓ia)^{−s} = a^{−s} e^{±iπs/2}, so I\_ε(a) → Γ(s) a^{−s} cos(πs/2).
- **Domination in ξ.** The majorant |Γ(s)| e^{π|t|/2} (2ξ/π)^{−σ}/(2 cosh ξ) is integrable exactly when σ < 1. Dominated convergence therefore gives

  Γ(s) F(s) = (2ω̄/π) Γ(s) cos(πs/2) ∫₀^∞ (2ξ/π)^{−s} e^{−iξ²/π}/(2 cosh ξ) dξ = (2ω̄/π)(π/2)^{s} Γ(s) cos(πs/2) A(1 − s).

- Divide by Γ(s), which is finite and nonzero in the strip. ∎

## 5. Reflection

**Lemma 8.** For all s ∈ ℂ, F̃(s) = −conj F(1 − s̄).

*Proof.*

- For w = c + uω̄ ∈ ℓ̃\_c: w̄ = c + uω ∈ ℓ\_c, Ψ\*(w) = conj Ψ(w̄), and w^{s−1} = conj(w̄^{s̄−1}), since log w̄ = conj log w off the cut. Hence

  F̃(s) = ∫\_{−∞}^{∞} conj((c + uω)^{s̄−1} Ψ(c + uω)) ω̄ du = conj ∫\_{−∞}^{∞} (c + uω)^{s̄−1} Ψ(c + uω) ω du.

- By definition, with the orientation +∞ → −∞, F(1 − s̄) = −∫\_{−∞}^{∞} (c + uω)^{s̄−1} Ψ(c + uω) ω du. ∎

**Corollary 9.** For 0 < σ < 1,

F̃(s) = −(2ω/π)(π/2)^{1−s} sin(πs/2) · B(s).

*Proof.* 1 − s̄ also lies in the strip. Apply Prop. 7 at 1 − s̄ and conjugate:

- conj A(s̄) = B(s);
- conj (π/2)^{1−s̄} = (π/2)^{1−s};
- conj cos(π(1 − s̄)/2) = cos(π(1 − s)/2) = sin(πs/2);
- conj ω̄ = ω. ∎

## 6. A Γ-identity

**Lemma 10.** As meromorphic functions,

X(s) · (2/π)(π/2)^{1−s} sin(πs/2) · Γ(s) = 1.

*Proof.*

- Legendre duplication gives Γ(s) = 2^{s−1} π^{−1/2} Γ(s/2) Γ((s + 1)/2). Reflection gives Γ(s/2) Γ(1 − s/2) = π/sin(πs/2), and Γ((2 − s)/2) = Γ(1 − s/2). So

  X(s)Γ(s) = (π/4)^{s−1/2} 2^{s−1} π^{1/2} / sin(πs/2).

- Multiplying by (2/π)(π/2)^{1−s} sin(πs/2) gives

  π^{s−1/2} 2^{1−2s} · 2^{s−1} · π^{1/2} · 2π^{−1} · π^{1−s} 2^{s−1} = π^0 2^0 = 1. ∎

## 7. Proof of the Theorem

**In the strip 0 < σ < 1.** By Corollaries 5 and 9,

F − X F̃ = L − ωB/Γ(s) + X(s)(2ω/π)(π/2)^{1−s} sin(πs/2) B(s)
= L + (ωB(s)/Γ(s)) · [X(s)(2/π)(π/2)^{1−s} sin(πs/2) Γ(s) − 1] = L,

by Lemma 10.

**Everywhere.**

- F and F̃ are entire (P2), and L is entire.
- X is meromorphic, with simple poles exactly at s = 2, 4, 6, …: these come from Γ((2 − s)/2), while 1/Γ((s + 1)/2) is entire.
- So F − XF̃ is meromorphic on ℂ, and it equals L on a nonempty open set. By the identity theorem, L = F − XF̃ on ℂ minus {2, 4, …}.
- Since L is entire, X F̃ is regular at s = 2k. X has a genuine simple pole there, with (π/4)^{2k−1/2} ≠ 0 and Γ(k + 1/2) finite and nonzero. Hence F̃(2k) = 0 for k ≥ 1.
- Finally, by (P3), everything is independent of c ∈ (0, 1). ∎

## 8. Remarks

1. **The Riemann–Siegel form** (derivation.md §2.4) follows at once. Apply (P4) at every odd n ≤ N (the same argument with E(w) = w^{−s}, the strip ℓ\_{c}…ℓ\_{c\_N} staying off the cut), and the analogue for F̃, where Res Ψ\* = −χ(n)/(2πi). This gives

   L = Σ\_{n≤N} χ(n) n^{−s} + X Σ\_{n≤N} χ(n) n^{s−1} + R\_N.

2. **Where the Mordell identity is used:** only in Lemma 3. Lemma 6 reuses Lemma 3 on the line through 0.
3. **Literature.** The identity should be the q = 4 case of Siegel's Riemann–Siegel formula for L(s, χ) (Ann. of Math. 44, 1943, doi:10.2307/1968761). That text was not available, so the correspondence is unverified. Titchmarsh §2.10 was not consulted either. The proof above is self-contained.
4. **The general real primitive χ, and every self-dual periodic sequence,** is proved in proof-rs-general.md (Theorem A′), with a symmetrised kernel. The new input there is the closed form of Φ\_K (Lemma 3′), obtained from K's own quasi-periodicity, and the self-duality lemma (G4), where the Gauss sum enters.

## 9. Numeric spot checks (interactive mpmath; illustrations, not reproduced by a stored script)

All values are relative errors |LHS − RHS|/|RHS|. The Mellin integrals were computed in the rotated, Gaussian-damped forms. B(s) = ω^{s} ∫₀^∞ y^{s−1} e^{−y²/π}/(2 cosh(ωy)) dy and A(u) = ω̄^{u} ∫₀^∞ …/(2 cosh(ω̄y)) dy, by the rotation of Prop. 4(iii); the head ∫₀^{1/4} was done by an exact Taylor series in y.

| Identity                                             | Points                                                     | Rel. error                                                 |
| ---------------------------------------------------- | ---------------------------------------------------------- | ---------------------------------------------------------- |
| Lemma 3 (Φ closed form), c = 0.5 and c = 0.1         | z = 0.3; 2.2 − 1.1i; 0.8276i                               | 4.9e−34; 1.1e−33; 1.2e−33 (the same for both c)            |
| Lemma 3 at dps 50                                    | z = −1.7 + 0.6i                                            | 3.6e−34 (dps 32 gave 3.0e−17: cancellation, not a failure) |
| Lemma 2 (Gaussian line integral, as used in Lemma 3) | z = 0.7 + 0.4i                                             | 3.9e−33                                                    |
| (P4) shift/residue: Φ\_c − Φ\_{c+2} = e^{iπz/2}        | z = 0.7 + 0.4i                                             | 8.1e−34                                                    |
| Lemma 6 (m(x) = πω h(x), rotated form)               | x = 0.7; 3                                                 | 2.2e−35; 2.0e−35                                           |
| Lemma 6, direct oscillatory ξ-integral vs rotated    | x = 0.7; 3                                                 | 1.6e−33; 5.9e−33                                           |
| Corollary 5 (Γ F = Γ L − ωB)                         | s = 0.3 + 5i; 0.7 − 3i                                     | 2.5e−37; 8.6e−37                                           |
| Corollary 5                                          | s = 2.5 + i                                                | 7.8e−34                                                    |
| Proposition 7                                        | s = 0.8 + 2i; 0.4 + 7i                                     | 1.3e−37; 7.9e−37                                           |
| Proposition 7                                        | s = 0.2 − i                                                | 4.3e−30                                                    |
| Corollary 9                                          | s = 0.4 + 7i; 0.6 − 2i                                     | 7.5e−37; 9.3e−38                                           |
| Lemma 8 (F̃(s) = −conj F(1 − s̄))                      | s = 0.4 + 7i; 0.6 − 2i                                     | 0.0 (structural: the same quadrature nodes)                |
| Lemma 10 (Γ-identity)                                | s = 0.37 + 4.1i; −2.3 + 0.7i; 5.5 − 3i                     | 1.7e−33; 3.4e−33; 4.1e−33                                  |
| Γ(s)L(s) = ∫ x^{s−1}/(2 cosh x) dx                   | s = 0.3 + 5i                                               | 2.1e−9 (naive oscillatory quadrature only)                 |
| Theorem, whole identity (a second, independent check) | s = ½ + 20i; 2.5 + 3i; −1.5 + 7i; 0.3 − 15i; 0.5 + 0.5i; 3 | 5.7e−31; 3.3e−32; 4.9e−31; 2.1e−23; 1.4e−31; 1.1e−31       |
| F̃(2k) = 0 (same run)                                 | k = 1, 2, 3                                                | \|F̃\| ≈ 5e−35, 1e−35, 2e−36                                |

In an earlier pass, the Prop. 4 check with naive oscillatory quadrature on [0, ∞) agreed only to 1e−6 to 1e−8. That limit was the quadrature at x → 0 with large |t|, cured by the rotated form; it is not a failure of the identity.
