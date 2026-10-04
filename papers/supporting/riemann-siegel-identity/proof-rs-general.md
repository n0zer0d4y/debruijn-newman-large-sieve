# The Riemann–Siegel integral identity for every real primitive χ: proof and a correction

This extends `proof-rs-chi4.md`, cited below as [χ₋₄]. It is a written proof, not formalised. The numeric spot checks of §9 are illustrations: they were computed interactively with mpmath, and only the c₀ census (`code/03-davenport-heilbronn/c0_census.py`) and the Davenport–Heilbronn instance of Theorem A′ (`code/03-davenport-heilbronn/rs_identity_dh.py`) are reproduced by stored scripts.

## Summary

- **The naive identity is false in general.** With the kernel K defined below, the dual kernel conj K(−w̄) and the sign (−1)^𝔞,
  - L = F\_K + (−1)^𝔞 X F̃\_K fails for D = 8. The relative error is about 8.66 at s = 1.7 − 2i and 1.3·10⁻¹⁰ at s = ½ + 7i (interactive check). The defect decays like e^{−πt} as t → +∞, which is why a check at s = ½ + 15i alone cannot see it.
  - The defect is exactly the explicit term Δ(s) of Corollary B. It vanishes precisely when c₀ = 0 (defined below).
- **Theorem A (the correct general identity).** Replace K by the symmetrised kernel 𝒦 := K + (c₀/2)e^{iπw²/q}. It has the same residues as K and parity 𝒦(−w) = (−1)^{𝔞+1}𝒦(w). Then L = F\_𝒦 + (−1)^𝔞 X F̃\_𝒦 holds for all s.
  - Here F̃\_𝒦 uses conj 𝒦(−w̄), which equals (−1)^{𝔞+1} conj 𝒦(w̄). So for 𝒦 the "naive" and the "−w̄" dual kernels agree up to the sign (−1)^{𝔞+1}.
- **When c₀ ≠ 0.** c₀ = 0 automatically for odd χ, and for q odd. Among the 242 fundamental discriminants with |D| ≤ 400, c₀ ≠ 0 exactly for **D = 8, 40, 104, 136, 232, 296, 328**. In particular χ₋₄ and D = 5, 12, 24, −8, −3 are unaffected, since there K = 𝒦 (census: `code/03-davenport-heilbronn/c0_census.py`).
- **The Gauss sum** enters through the self-duality (Lemma G4) conj(c\_j e^{−iπj²/(4q)}) = ω i^{−𝔞} c\_j. That lemma uses τ(χ) = i^𝔞√q and the Landsberg–Schaar reciprocity.

## Notation

- χ is a real primitive character mod q ≥ 3, with parity 𝔞 ∈ {0, 1}, so χ(−1) = (−1)^𝔞.
- ω = e^{iπ/4}; G(w) := e^{iπw²/q}; ζ := e^{iπw/q}.
- The coefficients: c\_j := (1/2q) Σ\_{n mod 2q} χ(n) e^{−iπn²/q} e^{−iπjn/q} for 0 ≤ j < 2q. This is the unique solution of e^{iπn²/q} Σ\_j c\_j e^{iπjn/q} = χ(n), n mod 2q, since it is a DFT on ℤ/2q.
- The kernel: K(w) := G(w) Σ\_j c\_j ζ^j/(ζ^{2q} − 1).
- c₀ is the j = 0 coefficient, and c′\_j := c\_j e^{−iπj²/(4q)}.
- **Lines and branch** as in [χ₋₄]:
  - ℓ\_c = {c + uω} is traversed u: +∞ → −∞ ("c↙"); ℓ̃\_c = {c + uω̄} is traversed u: −∞ → +∞ ("c↘").
  - Take 0 < c < 1, with the principal branch of w^{±s}. On ℓ\_c, Re(ωw) = c/√2, arg w ∈ (−3π/4, π/4), and log(ωw) = log w + iπ/4.
- **For a kernel k**, F\_k(s) := ∫\_{c↙} w^{−s} k(w) dw and F̃\_k(s) := ∫\_{c↘} w^{s−1} conj k(−w̄) dw.
- X(s) := (π/q)^{s−1/2} Γ((1 − s + 𝔞)/2)/Γ((s + 𝔞)/2). The functional equation is L(s) = X(s) L(1 − s), with root number 1 for real primitive χ.
- **Transforms.**
  - θ₀(x) := Σ\_{n≥1} χ(n) e^{−nx}. Then ∫₀^∞ x^{s−1} θ₀ dx = Γ(s) L(s) for σ > 0, by the argument of [χ₋₄] Cor. 5.
  - For σ > 0, 𝔅\_K(s) := ∫₀^∞ x^{s−1} g\_K(x) dx, where g\_K is defined in Lemma 3′.
  - 𝔅♯\_K(u) := conj 𝔅\_K(ū).

**Theorem A.** For all s ∈ ℂ, with 𝒦 = K + (c₀/2)G,

L(s, χ) = ∫\_{c↙} w^{−s} 𝒦(w) dw + (−1)^𝔞 X(s) ∫\_{c↘} w^{s−1} conj 𝒦(−w̄) dw.

**Corollary B.** For the kernel K,

F\_K(s) + (−1)^𝔞 X(s) F̃\_K(s) = L(s, χ) + Δ(s),
Δ(s) := (q c₀ ω/(4π)) (2π/q)^s (4π/q)^{(1−s)/2} Γ((1 − s)/2) e^{3πis/4}.

Here Δ ≡ 0 if and only if c₀ = 0, and c₀ = 0 whenever 𝔞 = 1 or q is odd (Lemma G3).

## 1. Algebraic lemmas on K

**(G1) Partial fractions.**

K(w) = (1/2q) Σ\_{n mod 2q} χ(n) e^{iπ(w² − n²)/q}/(e^{iπ(w−n)/q} − 1).

*Proof.* Substitute c\_j and sum the geometric series Σ\_{j<2q} (ζ e^{−iπn/q})^j = (ζ^{2q} − 1)/(ζ e^{−iπn/q} − 1). Each summand is invariant under n ↦ n + 2q.

*Consequences.*

- The poles of K lie at the integers, with Res\_{w=m} K = χ(m)/(2πi); near m ≡ n (mod 2q), e^{iπ(w−n)/q} − 1 ≈ (iπ/q)(w − m).
- K is regular at w = 0, because χ(0) = 0.
- On ℓ\_a, |G(w)| = exp(−(π/q)(u² + √2 au)). The sum Σ\_j c\_j ζ^j/(ζ^{2q} − 1) is bounded on ℓ\_a when a ∉ ℤ, because every exponent j ≤ 2q − 1 is below 2q.
- Hence (P1)–(P3) of [χ₋₄] hold verbatim for K: Gaussian decay e^{−πu²/q + O(|u|)}, F\_K and F̃\_K entire, and independence of c. (P4) holds for moving ℓ\_c past the integers m, each contributing χ(m) E(m).

**(G2) Quasi-periodicity.** K(w + q) = e^{2πiw} K(w).

*Proof.* In (G1) put n = m + q: χ(m + q) = χ(m), (w + q)² − (m + q)² = w² − m² + 2q(w − m), and e^{2πi(w−m)} = e^{2πiw}. ∎

Hence K(w) − K(w + q) = −G(w) Σ\_j c\_j ζ^j.

**(G3) Reflection in w.**

K(−w) = −(−1)^𝔞 (K(w) + c₀ G(w)), with c₀ = (1/2q) Σ\_{n mod 2q} χ(n) e^{−iπn²/q}.

Moreover c₀ = 0 if 𝔞 = 1 or q is odd.

*Proof.*

- In (G1), K(−w) is a sum of terms χ(n) e^{iπ(w² − n²)/q}/(e^{−iπ(w+n)/q} − 1). Replace n by −n, using χ(−n) = (−1)^𝔞 χ(n).
- Use 1/(e^{−iu} − 1) = −1 − 1/(e^{iu} − 1).
- For c₀: n ↦ −n gives c₀ = (−1)^𝔞 c₀. For q odd, n ↦ n + q multiplies the summand by e^{−iπ(2n + q)} = −1. ∎

So 𝒦 := K + (c₀/2)G satisfies 𝒦(−w) = (−1)^{𝔞+1} 𝒦(w), and it has the same residues as K (G is entire).

**(G4) Self-duality: where the Gauss sum enters.** For every j,

conj(c\_j e^{−iπj²/(4q)}) = ω i^{−𝔞} c\_j.

*Proof.*

- The left side is (1/2q) S\_j, with S\_j := Σ\_{n mod 2q} χ(n) e^{iπ(n + j/2)²/q}.
- Write χ(n) = τ^{−1} Σ\_{m mod q} χ(m) e^{2πimn/q}. This holds for all n, since χ is primitive and real, and τ = τ(χ) = i^𝔞 √q.
- Complete the square: (n + j/2)² + 2mn = (n + m + j/2)² − m² − mj. Hence

  S\_j = τ^{−1} Γ\_j Σ\_{m mod q} χ(m) e^{−iπ(m² + mj)/q}, with Γ\_j := Σ\_{n mod 2q} e^{iπ(n + j/2)²/q}.

- **Γ\_j by Landsberg–Schaar.** If q + j is odd, then n ↦ n + q flips the sign of the summand, so Γ\_j = 0. If q + j is even, the summand has period q. The Landsberg–Schaar identity Σ\_{n<q} e^{πi(n² + jn)/q} = √q e^{πi(q − j²)/(4q)} (with (a, b, c) = (1, j, q) and ac + b even) then gives Γ\_j = 2√q ω.
- **The m-sum is q c\_j.** For q + j even it equals ½ Σ\_{m mod 2q} χ(m) e^{−iπ(m² + mj)/q} = q c\_j. For q + j odd the same pairing shows c\_j = 0.
- Hence S\_j = 2q√q ω c\_j/τ in all cases, and the claim follows from √q/τ = i^{−𝔞}. ∎

## 2. The closed form of Φ\_K (the analogue of Lemma 3)

**Lemma 3′.** Let Φ\_K(z) := ∫\_{c↙} K(w) e^{2πizw/q} dw, with 0 ≤ c < 1. It is entire in z, independent of c, and

Φ\_K(z) = [P(z) − ω√q e^{2πiz} C(z)] / (1 − e^{2πiz}),

where P(z) := Σ\_{0<n<q} χ(n) e^{2πizn/q} and C(z) := Σ\_j c\_j e^{−iπ(j + 2z)²/(4q)}.

*Proof.*

- **The Gaussian line integral.** Exactly as in [χ₋₄] Lemma 2, ∫\_{c↙} e^{iπ(w+b)²/q} dw = −ω√q, because the integrand is e^{−π(u+β)²/q} after c + b + uω = ω(u + β).
- **Integrate (G2).** Multiply K(w) − K(w + q) = −G(w) Σ\_j c\_j ζ^j by e^{2πizw/q} and integrate over ℓ\_c.
  - The substitution v = w + q together with (P4) gives ∫\_{c↙} K(w + q) e^{2πizw/q} dw = e^{−2πiz} (Φ\_K(z) − Σ\_{n=1}^{q} χ(n) e^{2πizn/q}), where χ(q) = 0. The poles crossed are n = 1, …, q.
  - The right side is −Σ\_j c\_j ∫ e^{(iπ/q)(w² + (j + 2z)w)} dw = ω√q Σ\_j c\_j e^{−iπ(j + 2z)²/(4q)}.
- So (1 − e^{−2πiz}) Φ\_K = −e^{−2πiz} P + ω√q C. Multiply by −e^{2πiz}. ∎

**The specialisation z = iqx/(2π), x > 0.** Here e^{2πizn/q} = e^{−nx}, e^{2πiz} = e^{−qx}, and e^{−iπ(j + 2z)²/(4q)} = e^{−iπj²/(4q)} e^{jx/2} e^{iqx²/(4π)}. Hence

h\_K(x) := Φ\_K(iqx/(2π)) = θ₀(x) + g\_K(x), g\_K(x) := −ω√q e^{iqx²/(4π)} Σ\_j c′\_j e^{jx/2}/(e^{qx} − 1).

- θ₀ = P/(1 − e^{−qx}) is regular at 0, since Σ\_{n mod q} χ(n) = 0. As h\_K is entire, g\_K is regular at 0 too.
- g\_K(x) = O(e^{−x/2}) as x → +∞, since j ≤ 2q − 1.
- (G4) rewrites g\_K in terms of K itself:

  **conj g\_K(ξ) = −i^{−𝔞} √q · K(−iqξ/(2π)) for all real ξ.**

  To see this, note K(−iqξ/(2π)) = e^{−iqξ²/(4π)} Σ\_j c\_j e^{jξ/2}/(e^{qξ} − 1), then conjugate g\_K and apply (G4).

**The pure Gaussian.** For G itself, the Gaussian line integral gives Φ\_G(z) := ∫\_{c↙} G(w) e^{2πizw/q} dw = −ω√q e^{−iπz²/q}.

## 3. F\_K and F\_G as Mellin transforms (the analogue of Prop. 4)

**Proposition 4′.** For σ > 0:

Γ(s) F\_K(s) = ∫₀^∞ x^{s−1} h\_K(x) dx = Γ(s) L(s) + 𝔅\_K(s),

F\_G(s) = −(ω√q/2) e^{iπs/4} Γ(s/2) (4π/q)^{s/2}/Γ(s).

The second formula holds for all s by continuation, since Γ(s/2)/Γ(s) = 2^{1−s}√π/Γ((s+1)/2) is entire.

*Proof.*

- **The Laplace step** of [χ₋₄] Prop. 4(i)(ii) is verbatim, with e^{−ωwy} = e^{2πizw/q} for z = iqωy/(2π). It gives Γ(s) F\_k(s) = e^{iπs/4} ∫₀^∞ y^{s−1} Φ\_k(iqωy/(2π)) dy for k = K, G.
- **For G, no rotation is needed.** Φ\_G(iqωy/(2π)) = −ω√q e^{−qy²/(4π)} is a real Gaussian in y, and ∫₀^∞ y^{s−1} e^{−qy²/(4π)} dy = ½ Γ(s/2)(4π/q)^{s/2}.
- **For K, rotate** the ray from arg y = 0 to arg y = −π/4, as in [χ₋₄] Prop. 4(iii). By Lemma 3′, on the sector −π/4 ≤ arg y ≤ 0:
  - |Φ\_K(iqωy/(2π))| ≤ C e^{−Re(ωy)/2} for |y| ≥ 1. The P-part is Σ χ(n) e^{−nωy}/(1 − e^{−qωy}). In the Gaussian part, |e^{−iπz²/q}| = e^{−q Re(y²)/(4π)} ≤ 1 and |e^{2πiz} e^{−iπjz/q}| = e^{(j/2 − q) Re(ωy)} ≤ e^{−Re(ωy)/2}.
  - The denominator 1 − e^{−qωy} ≠ 0 for y ≠ 0 in the sector.
  - So the arcs vanish, as before. On arg y = −π/4, y = ω̄x gives e^{iπs/4} ω̄^s = 1 and z = iqx/(2π).
- Finally θ₀ contributes Γ(s)L(s). ∎

## 4. h\_K as a cosine or sine transform (the analogue of Lemma 6)

**Lemma 6′.** For real x > 0, with λ := qx/(2π):

- (𝔞 = 1, so c₀ = 0) h\_K(x) = −(√q/π) ∫₀^∞ conj g\_K(ξ) cos(λξ) dξ;
- (𝔞 = 0) h\_K(x) = −(√q/π) ∫₀^∞ conj g\_K(ξ) sin(λξ) dξ + (iqc₀/(2π)) E(λ), where

  E(λ) := ∫₀^∞ e^{−iqξ²/(4π)} e^{−iλξ} dξ = ω̄ ∫₀^∞ e^{−qr²/(4π)} e^{−iλω̄r} dr.

  The second form is absolutely convergent: rotate ξ = ω̄r on the sector −π/4 ≤ arg ξ ≤ 0.

*Proof.*

- **Rotate ℓ₀ onto iℝ.** By Lemma 3′, h\_K(x) = ∫\_{0↙} K(w) e^{−xw} dw. Rotate the upper half-ray (arg π/4) onto arg π/2 and the lower half-ray (arg 5π/4) onto arg 3π/2, separately. The two sectors contain no pole.
  - **Arcs of radius R.** |G| = e^{−(π/q)R² sin 2φ} and |e^{−xw}| = e^{−xR cos φ}. In the upper sector cos φ ≥ 0. In the lower sector sin 2φ ≥ 0, and near φ = 3π/2 the Gaussian e^{−2πR²δ/q} beats e^{xRδ} once R > xq/(2π), as in [χ₋₄] Lemma 6. The kernel quotient is bounded. So the arc contributions are O(1/R).
  - **The case c₀ ≠ 0.** Then K(iy) → −c₀ e^{−iπy²/q} as y → +∞, a pure chirp. On the upper sector split K = [K + c₀G] − c₀G. The bracket is O(|G| e^{−π Im w/q}), so it gives an absolutely convergent ray integral. For −c₀G e^{−xw} the arc estimate above shows the imaginary-ray integral exists as an improper integral and equals the rotated one.
- **Parametrise iℝ** by w = −iqξ/(2π), with ξ: −∞ → +∞ (w from +i∞ to −i∞, matching the orientation of ℓ₀). Using (G4),

  h\_K(x) = (i^{𝔞+1}√q/(2π)) ∫\_ℝ conj g\_K(ξ) e^{iλξ} dξ,

  an improper integral at ξ → −∞ when c₀ ≠ 0.

- **Fold onto (0, ∞).** By (G3) and (G4), conj g\_K(−ξ) = −(−1)^𝔞 conj g\_K(ξ) + (−1)^𝔞 i^{−𝔞} √q c₀ e^{−iqξ²/(4π)}.
  - For 𝔞 = 1, conj g\_K is even, which gives the cosine form (i² = −1).
  - For 𝔞 = 0, the odd part gives 2i ∫₀^∞ conj g\_K sin, and the chirp gives √q c₀ E(λ).
- The integrals over (0, ∞) of conj g\_K are absolutely convergent, since g\_K = O(e^{−ξ/2}). The E-integral equals its rotated form by the same arc estimate. ∎

For q = 4 this is exactly [χ₋₄] Lemma 6.

## 5. The dual Mellin form (the analogue of Prop. 7)

**Proposition 7′.** For 0 < σ < 1, with cs\_𝔞 = sin (𝔞 = 0) or cos (𝔞 = 1):

∫₀^∞ x^{s−1} h\_K(x) dx = −(√q/π)(2π/q)^s Γ(s) cs\_𝔞(πs/2) 𝔅♯\_K(1 − s) + [𝔞 = 0] (iqc₀/(2π)) (2π/q)^s 𝓜E(s),

𝓜E(s) = (ω̄/2) e^{−iπs/4} Γ(s) Γ((1 − s)/2) (4π/q)^{(1−s)/2}.

*Proof.*

- **The transform term.** This is the e^{−εx}-regularised Fubini argument and dominated convergence of [χ₋₄] Prop. 7, with ∫₀^∞ x^{s−1} e^{−εx} sin(ax) dx = (Γ(s)/(2i))[(ε − ia)^{−s} − (ε + ia)^{−s}] → Γ(s) a^{−s} sin(πs/2). The majorant ξ^{−σ}|g\_K(ξ)| is integrable for σ < 1.
- **𝓜E.** Use the rotated form of E. Fubini holds, because ∫∫ λ^{σ−1} e^{−λr/√2} e^{−qr²/(4π)} dλ dr = Γ(σ)2^{σ/2} ∫ r^{−σ} e^{−qr²/(4π)} dr < ∞ for 0 < σ < 1. Then ∫₀^∞ λ^{s−1} e^{−iω̄rλ} dλ = Γ(s)(ωr)^{−s} (iω̄ = ω).
- **The scaling.** 𝓜[E(qx/(2π))](s) = (2π/q)^s 𝓜E(s). ∎

## 6. Reflection and the Γ-identity (the analogues of Lemma 8, Cor. 9 and Lemma 10)

**Lemma 8′.** For any kernel k, F̃\_k(s) = −conj F\_{k⁻}(1 − s̄), where k⁻(v) := k(−v). The proof is that of [χ₋₄] Lemma 8, with k⁻(w̄) in place of Ψ(w̄).

With (G3), for all s:

F̃\_K(s) = (−1)^𝔞 [F♯\_K(1 − s) + conj(c₀) F♯\_G(1 − s)], F̃\_G(s) = −F♯\_G(1 − s),

where F♯(u) := conj F(ū).

**Corollary 9′.** For 0 < σ < 1, F♯\_K(1 − s) = L(1 − s) + 𝔅♯\_K(1 − s)/Γ(1 − s). This is Prop. 4′ at 1 − s̄, conjugated, using that χ is real.

**Lemma 10′ (Γ-identity with parity 𝔞).**

X(s) Γ(s)/Γ(1 − s) = (√q/π) (2π/q)^s Γ(s) cs\_𝔞(πs/2).

Equivalently, X(s)/Γ(1 − s) = (√q/π)(2π/q)^s cs\_𝔞(πs/2).

*Proof.* Use duplication for Γ(s) and Γ(1 − s), together with Γ(s) sin(πs/2) = 2^{s−1}√π Γ((s+1)/2)/Γ(1 − s/2) (𝔞 = 0), or Γ(s) cos(πs/2) = 2^{s−1}√π Γ(s/2)/Γ((1 − s)/2) (𝔞 = 1). The powers collect to (π/q)^{s−1/2} 2^s π^{−1/2} = (√q/π)(2π/q)^s. For q = 4, 𝔞 = 1 this is [χ₋₄] Lemma 10. ∎

## 7. Proof of Corollary B and Theorem A

**Assembling.** Fix s with 0 < σ < 1. Since (−1)^𝔞 (−1)^𝔞 = 1,

F\_K + (−1)^𝔞 X F̃\_K = [L + 𝔅\_K/Γ(s)] + X[L(1 − s) + 𝔅♯\_K(1 − s)/Γ(1 − s) + conj(c₀) F♯\_G(1 − s)].

- XL(1 − s) = L, and 𝔅\_K/Γ(s) = (∫ x^{s−1} h\_K)/Γ(s) − L.
- By Prop. 7′ and Lemma 10′, the 𝔅♯\_K-terms cancel exactly. What remains is

  F\_K + (−1)^𝔞 X F̃\_K − L = [𝔞 = 0] { (iqc₀/(2π)) (2π/q)^s 𝓜E(s)/Γ(s) + X conj(c₀) F♯\_G(1 − s) }.

**Evaluating the remainder.**

- By (G4) at j = 0, conj c₀ = ω c₀ when 𝔞 = 0.
- By Prop. 4′, F♯\_G(1 − s) = −(ω̄√q/2) e^{−iπ(1−s)/4} Γ((1 − s)/2)(4π/q)^{(1−s)/2}/Γ(1 − s).
- Insert these and Lemma 10′ (𝔞 = 0). The braces become (q c₀/(4π))(2π/q)^s (4π/q)^{(1−s)/2} Γ((1 − s)/2) · [ω e^{−iπs/4} − 2 sin(πs/2) ω̄ e^{iπs/4}], and the bracket equals ω e^{3iπs/4}. This is Δ(s).
- **Corollary B** follows on the strip, and on all of ℂ by the identity theorem. All terms are meromorphic, L and F\_K, F̃\_K are entire, and Δ is meromorphic.

**Theorem A.**

- F\_𝒦 = F\_K + (c₀/2) F\_G and F̃\_𝒦 = F̃\_K + (conj c₀/2) F̃\_G = F̃\_K − (conj c₀/2) F♯\_G(1 − s).
- So F\_𝒦 + (−1)^𝔞 X F̃\_𝒦 − L = Δ + [𝔞 = 0]{(c₀/2) F\_G(s) − (X/2) conj(c₀) F♯\_G(1 − s)}. For 𝔞 = 1 the c₀-terms vanish.
- For 𝔞 = 0, the same substitutions give

  Δ − (X/2) ω c₀ F♯\_G(1 − s) = (q c₀/(4π))(2π/q)^s (4π/q)^{(1−s)/2} Γ((1 − s)/2) · ω e^{iπs/4} cos(πs/2),

  while (c₀/2) F\_G(s) = −(c₀ ω√q/4) e^{iπs/4} Γ(s/2)(4π/q)^{s/2}/Γ(s).

- Their sum vanishes if and only if Γ((1 − s)/2) cos(πs/2) = (√π/2) 2^s Γ(s/2)/Γ(s) = π/Γ((1 + s)/2), which is the reflection formula. ∎

**Corollary (poles of X).** The poles of X are at s ≡ 1 + 𝔞 (mod 2), s ≥ 1 + 𝔞. Since L is entire, F̃\_𝒦 vanishes at each of them.

## 8. Remarks

1. **The Riemann–Siegel form.** Shifting ℓ\_c past the integers n ≤ N gives L = Σ\_{n≤N} χ(n) n^{−s} + X Σ\_{n≤N} χ(n) n^{s−1} + R\_N(s) with the kernel 𝒦, as used in derivation.md §2.4. For 𝒦 the dual residues are compatible with the sign (−1)^𝔞.
2. **Why a single check on the line can miss the defect.** At s = σ + it, |Δ(s)| ≍ |t|^{−σ/2} e^{−π|t|/4} e^{−3πt/4}. As t → +∞ this is about e^{−πt}, which is invisible at t = 15. As t → −∞ it grows like e^{π|t|/2}. Unlike the true identity, the K-form is not symmetric under s ↦ s̄, so it fails badly at s = ½ − 15i. Test identities at low |t|, at both signs of t, and off the line.
3. **For the Dirichlet-character case.** χ₋₄ and the even character mod 5 have c₀ = 0, so for them K = 𝒦.
4. **Novelty.** Siegel's 1943 formula (Ann. of Math. 44) was not read, so whether it uses the symmetric kernel is unverified.
5. **Beyond characters.**
   - *What the proof uses.* It never uses multiplicativity, primitivity or the character property as such. It uses a: ℤ → ℝ only through:
     - (i) q-periodicity;
     - (ii) parity, a(−n) = (−1)^𝔞 a(n);
     - (iii) a(0) = 0;
     - (iv) Σ\_{n mod q} a(n) = 0;
     - (v) real values;
     - (vi) the DFT self-duality Σ\_{m mod q} a(m) e^{2πimn/q} = i^𝔞 √q · a(n) for all n. This is the only input of Lemma G4, which in §1 enters through τ(χ) = i^𝔞√q.
   - *Theorem A′.* Theorem A, Corollary B and Remark 1 hold verbatim for every such a, with L(s) := Σ a(n) n^{−s}. The functional equation L(s) = X(s)L(1−s) then follows from the identity itself.
   - *Example.* The Davenport–Heilbronn coefficients (1, κ, −κ, −1, 0) mod 5 satisfy (vi) exactly: with u = √(10−2√5) and v = √(10+2√5), it reduces to 2u = (√5−1)v and 2v = (√5+1)u, which hold because both sides are positive and their squares are equal (40 ∓ 8√5). The Polymath-type machinery then gives Λ\_DH ≤ 0.58, conditional on the open items of derivation.md (papers/03-davenport-heilbronn-lambda, claim #4).
   - *The general family.* Consider the real (i^𝔞√q)-eigenspace of the DFT on sequences of parity 𝔞 with a(0) = 0 and Σ a = 0.
     - For 3 ≤ q ≤ 19 its dimension is ⌊(q+1)/4⌋ for 𝔞 = 1 and ⌊q/4⌋ for 𝔞 = 0 (numeric-experiment; `code/03-davenport-heilbronn/dh_family_dim.py`, SVD null space; output in `runs/dh_family_dim.txt`).
     - For q = 5 each parity gives a one-dimensional space: the Legendre symbol (even) and the Davenport–Heilbronn vector (odd).
     - From q = 7 (odd parity, dimension 2) on, the space contains combinations of several characters, such as the circle of papers/02-lambda-unbounded-functional-equation-family.
     - This is a family of non-multiplicative functions with a functional equation, to which the Λ machinery of this paper applies unchanged.

## 9. Numeric spot checks (interactive mpmath, dps 25; relative errors; not reproduced by a stored script except the census)

| Check                                                             | D = 5 (𝔞 = 0, c₀ = 0) | D = 8 (𝔞 = 0, c₀ = ½e^{−iπ/8}) | D = −3 (𝔞 = 1)   |
| ----------------------------------------------------------------- | --------------------- | ------------------------------ | ---------------- |
| (G2) K(w + q) = e^{2πiw}K(w)                                      | 9.7e−25               | 9.2e−25                        | 3.1e−25          |
| (G3) reflection                                                   | 2.3e−24               | 1.6e−25                        | 2.0e−25          |
| (G4) max\_j of the difference                                      | 4.9e−25               | 2.3e−25                        | 1.7e−25          |
| Lemma 3′, z = 0.3 + 0.2i and −1.1 + 0.4i                          | 3.7e−24, 9.3e−25      | 7.6e−25, 3.2e−25               | 1.9e−25, 1.3e−25 |
| h\_K = θ₀ + g\_K vs the line integral of Φ\_K(iqx/2π), x = 0.7 and 2 | 1.5e−24, 2.9e−25      | 2.5e−25, 1.0e−24               | 3.0e−25, 1.2e−25 |
| Lemma 6′ (sin/cos transform plus E), x = 0.7 and 2                | 2.8e−25, 3.1e−25      | 3.5e−25, 6.4e−25               | 1.0e−23, 1.4e−23 |
| F\_G closed form, s = 0.4 + 3i and 2.2 − i                         | 5.2e−26, 1.6e−26      | 3.6e−26, 2.7e−27               | 4.8e−26, 1.2e−26 |
| Γ-identity (Lemma 10′), same s                                    | 2.8e−26, 5.8e−26      | 3.2e−26, 4.0e−28               | 1.3e−26, 1.0e−25 |
| The naive identity with K, s = ½ + 7i and 1.7 − 2i                | 2.1e−24, 3.7e−24      | **1.3e−10, 8.66 (fails)**      | 3.2e−25, 2.8e−24 |

For D = 8 in detail, at s = 1.7 − 2i, 0.3 + i and ½ + 7i:

- the defect F\_K + X F̃\_K − L equals the predicted Δ(s) to 1.2e−24, 1.7e−25 and 1.5e−25 (relative to |L|). At s = 1.7 − 2i the defect is −8.8240491 + 4.2490888i;
- Theorem A with the symmetric kernel 𝒦 has errors 1.9e−25, 6.6e−26 and 1.7e−25.

For D = 12, 24 and −8 (all with c₀ = 0 to about 1e−25), the naive identity holds to 1e−23 to 1e−24 at the same points.

A c₀ census over the 242 fundamental discriminants with 0 < |D| ≤ 400 (`code/03-davenport-heilbronn/c0_census.py`, output in `runs/c0_census.txt`) gives c₀ ≠ 0 exactly for D = 8, 40, 104, 136, 232, 296 and 328. In every one of these cases arg c₀ ≡ −π/8 (mod π), consistent with conj c₀ = ω c₀.
