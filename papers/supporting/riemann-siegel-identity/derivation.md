# The Polymath-15 effective approximation for self-dual periodic Dirichlet series

This document collects the analytic inputs behind the upper bound of papers/03-davenport-heilbronn-lambda (claim #4) and behind Lemma 1′ of papers/02-lambda-unbounded-functional-equation-family/proof-semicontinuity.md. It transfers the effective approximation of D.H.J. Polymath [P19, Thm 1.3] from ζ to Dirichlet series with real, periodic, DFT-self-dual coefficients, the class of proof-rs-general.md Remark 5.

**Status.** Some lemmas are written in full and some are sketches. The open items are listed in §4. Every certificate that uses this document is therefore graded conditional.

## 1. Setting

- **Coefficients.** a: ℤ → ℝ is q-periodic, with parity 𝔞 ∈ {0, 1} (a(−n) = (−1)^𝔞 a(n)), a(0) = 0, Σ\_{n mod q} a(n) = 0, and the self-duality Σ\_{m mod q} a(m) e^{2πimn/q} = i^𝔞 √q · a(n). Put L(s) = Σ\_{n≥1} a(n) n^{−s}.
  - Real primitive characters are examples. So is the Davenport–Heilbronn function (q = 5, 𝔞 = 1), whose self-duality is proved exactly in papers/03-davenport-heilbronn-lambda/claims.md (#4).
- **Functional equation.** ξ(s) := (q/π)^{(s+𝔞)/2} Γ((s + 𝔞)/2) L(s) satisfies ξ(s) = ξ(1 − s). This follows from Theorem A′ (proof-rs-general.md).
- **Heat flow.** H\_0(z) := ξ((1 + iz)/2) is real on ℝ, even, and entire of order 1. H\_t(z) := (4πt)^{−1/2} ∫ e^{−w²/(4t)} ξ((1 + w + iz)/2) dw solves ∂\_t H = −∂\_z² H.
- **de Bruijn's strip theorem** [dB50]: if every zero of H\_0 lies in |Im z| ≤ Δ, then every zero of H\_t lies in |Im z| ≤ (max(Δ² − 2t, 0))^{1/2}. Hence a strip of half-width y₀ at time t₀ gives Λ ≤ t₀ + y₀²/2.
- Nothing below uses positivity of Φ, an Euler product, or a numerical verification of GRH. The certificate of papers/03-davenport-heilbronn-lambda is a fixed-time argument (H\_{t₀} has no zeros with |Im z| ≥ y₀), not Polymath's barrier criterion.

## 2. The effective Riemann–Siegel approximation

### 2.1 The normalising factor B\_t

Following [P19, (6)–(12)], with w = (s + 𝔞)/2:

- M\_0(s) := (q/π)^{w} √(2π) exp((w − ½) Log w − w), the Stirling model of (q/π)^{w}Γ(w). There is no pole, so ζ's factor s(s − 1)/2 is absent.
- α(s) := (log M\_0)′(s) = ½ log(q/π) + ½ Log w − 1/(4w) = ½ log(q/π) + ½ Log((s + 𝔞)/2) − 1/(2(s + 𝔞)).
- M\_t(s) := exp((t/4) α(s)²) M\_0(s), and B\_t(x + iy) := M\_t((1 + y − ix)/2).

The code uses exactly these formulas (`cert_final.Char.am`, `lM0`, `lMt`).

### 2.2 Main terms and error terms

Put s₊ = (1 + y − ix)/2 and s₋ = (1 − y + ix)/2. Then:

- b\_n = exp((t/4) log² n);
- s\* = s₊ + (t/2) α(s₊);
- κ = (t/2)(α(s₋) − α((1 + y + ix)/2));
- γ = M\_t(s₋)/M\_t(s₊);
- N = ⌊(q x/(4π) + q t/16)^{1/2}⌋.

The main term is

f\_t(x + iy) = Σ\_{n≤N} a(n) b\_n n^{−s\*} + γ Σ\_{n≤N} a(n) n^{y} b\_n n^{−conj(s\*)−κ}.

For real x the second exponent equals w := s₋ + (t/2)α(s₋), so f\_t is analytic in x in the form Σ a(n) b\_n n^{−s\*} + γ Σ a(n) b\_n n^{−w}. The conjugated s\* in the B-sum is what reproduces the Riemann–Siegel dual term n^{s−1} at t = 0.

The approximation to be justified is

H\_t(x + iy)/B\_t(x + iy) = f\_t(x + iy) + O\_≤(e\_A + e\_B + e\_C),

where:

- e\_A and e\_B come from the Dirichlet-term estimate (§2.3);
- e\_C comes from the Riemann–Siegel remainder (§2.5).

Assembly: H\_t is the Gaussian average of ξ. ξ is expanded by the Riemann–Siegel identity of §2.4 into the terms n ≤ N, their reflections, and two remainder integrals. Each term is then averaged as in [P19, §4].

### 2.3 The Dirichlet-term estimate (Prop 6.1^χ)

**Statement.** Let T ≥ 10, 0 < t ≤ ½ and σ real. The Gaussian average r\_{t,n}(σ + iT) of M\_0(s) n^{−s} (the n-th Dirichlet term with the Γ-factor) satisfies

r\_{t,n}(σ + iT) = M\_t(σ + iT) b\_n n^{−(σ+iT)−(t/2)α(σ+iT)} (1 + O\_≤(ε\_{t,n})),

ε\_{t,n} := exp( [(t²/8)|α(σ + iT) − log n|² + t/4 + 1/6] / (T − 1.25) ) − 1.

**Argument (sketch; open item O1).** This is [P19, Prop 6.1] with three changes:

- (a) Im α(σ + iT) = ½ arg((s + 𝔞)/2) + T/(2|s + 𝔞|²) > 0. So the shifted points stay in Im ≥ T, and Polymath's −0.15 correction (their (46)) is not needed.
- (b) α′(s) = 1/(2(s + 𝔞)) + 1/(2(s + 𝔞)²), so |α′| ≤ 1/(2T) + 1/(2T²) ≤ 1/(2(T − 1)) for Im s ≥ T. This replaces their (43).
- (c) (1 − t/(2(T − 1)))^{−1/2} ≤ exp(t/(4(T − 1.25))) for t ≤ ½. This replaces their (47).

The Stirling remainder for Γ((s + 𝔞)/2) is exp(O\_≤(1/(6(T − 0.66)))), by Boyd's bound [P19, Lemma 5.1(v)] applied with |Im((s + 𝔞)/2)| = T/2. It is absorbed into the denominator T − 1.25. The rest of Polymath's proof is claimed to transfer line by line, **but it has not been written out here**.

**Consequence.** e\_A + e\_B ≤ ε\_max · (Σ\_{n≤N} |a(n)| b\_n n^{−Re s\*} + |γ| Σ\_{n≤N} |a(n)| b\_n n^{−Re w}), with ε\_max the largest ε\_{t,n}. With |α| + log N ≤ 2 log N + 3, the code uses

ε\_max ≤ exp( [(t²/8)(2 log N + 3)² + t/4 + 1/6] / (x/2 − 1.25) ) − 1,

which decreases in x.

### 2.4 The Riemann–Siegel integral identity

proof-rs-general.md proves (Theorem A′) that, for every a as in §1,

L(s) = ∫\_{c↙} w^{−s} 𝒦(w) dw + (−1)^𝔞 X(s) ∫\_{c↘} w^{s−1} conj 𝒦(−w̄) dw,

where:

- 𝒦 = K + (c₀/2)e^{iπw²/q} is the symmetrised kernel, with c₀ = 0 when 𝔞 = 1 or q is odd;
- X(s) = (π/q)^{s−½}Γ((1 − s + 𝔞)/2)/Γ((s + 𝔞)/2);
- K(w) = e^{iπw²/q} Σ\_{j<2q} c\_j e^{πijw/q}/(e^{2πiw} − 1) has residue a(n)/(2πi) at every integer n.

proof-rs-chi4.md treats the case q = 4 separately, with the kernel Ψ(w) = ω e^{iπw²/4}/(4 cos(πw/2)).

Shifting the line past the integers n ≤ N (the residue theorem, given the Gaussian decay along the lines) gives

L(s) = Σ\_{n≤N} a(n) n^{−s} + (−1)^𝔞 X(s) Σ\_{n≤N} a(n) n^{s−1} + R\_N(s),

with R\_N the two remainder integrals on lines through c = N + ½.

For the Davenport–Heilbronn coefficients, code/03-davenport-heilbronn/rs\_identity\_dh.py checks the identity numerically.

### 2.5 The remainder: Lemma I and the explicit e\_C

**Notation.**

- s = σ + iT with T > 0.
- w\_R := (qT/(2π))^{1/2}, the real saddle of w^{−s} e^{iπw²/q}.
- c := N + ½ (a half-integer), δ := w\_R − c ∈ [−½, ½], λ := (w\_R/c)², so T = 2πc²λ/q.
- The line is w = c + u e^{iπ/4}, and J(σ, T, c) := ∫\_ℝ |w^{−s} e^{iπw²/q}| du along it.
- Q\_s := 4/q.

**Step 1 (the kernel).** On the line, |𝒦(w)| ≤ K\_h |e^{iπw²/q}|, with K\_h from Lemma H (§3). So |∫\_{c↙} w^{−s}𝒦(w) dw| ≤ K\_h J(σ, T, c).

**Step 2 (the exponent).** |e^{iπw²/q}| = e^{−(π/q)(√2cu + u²)}. With z = u e^{iπ/4}/c and y = u/(√2c), log|w^{−s}e^{iπw²/q}| = −σ log c − σ log|1 + z| + (2πc²/q) g(y), where g(y) := λ arg(1 + y + iy) − y − y².

**Step 3 (global bound).** For λ ≤ 1.5, g(y) ≤ (λ − 1)y − y²/2 for all real y. The cases:

- For y ≥ 0, use arg ≤ min(π/4, y).
- For y = −r < 0, put φ(r) = atan2(r, 1 − r):
  - r ≥ 2λ: use φ ≥ 0;
  - r ≤ 3/(2λ): use φ(r) ≥ atan r ≥ r − r³/3;
  - 1 ≤ r ≤ 2λ: use φ ≥ π/2 > λ/2 ≥ max\_r(r − r²/(2λ)).

Hence (2πc²/q) g(y) ≤ ℓ₀u − Q\_s π u²/8, with ℓ₀ := Q\_s π c(λ − 1)/(2√2). Also |σ log|1 + z|| ≤ |σ| m(u), with m(u) := max(½ log 2, log(1 + |u|/c)) ≤ |u|/c + 0.35.

For half-integer c ≥ 10.5 and |δ| ≤ ½, λ ≤ (1 + 1/21)² < 1.1.

**Step 4 (main region |u| ≤ ρc, ρ = ¼).** Write log(1 + z) = z − z²/2 + r(z), with |r(z)| ≤ |z|³/(3(1 − ρ)), Re z² = 0 and Im z² = 2y². Then

−σ log|1 + z| + (2πc²/q) g(y) = ℓu − κu² + ϑ·A|u|³, |ϑ| ≤ 1,

where:

- ℓ := ℓ₀ − σ/(√2c);
- κ := Q\_s π(1 + λ)/4;
- A := Q\_s πλ/(6c(1 − ρ)) + |σ|/(3c³(1 − ρ)).

**Lemma I.** Let c ≥ 10.5 be a half-integer, |δ| ≤ ½, ρ = ¼, and let σ be real with κ′ := κ − Aρc > κ/2. Then J(σ, T, c) ≤ c^{−σ}[G₀ + E + Tail], with:

- G₀ := (π/κ)^{1/2} e^{ℓ²/(4κ)};
- E ≤ A e^{ℓ²/(4κ′)} · 4(1/κ′² + |ℓ/(2κ′)|³ (π/κ′)^{1/2}), obtained by substituting u = v + ℓ/(2κ′) in A∫|u|³e^{ℓu − κ′u²}du, using e^P ≤ 1 + Pe^P and P ≤ Aρc·u² on the main region;
- Tail := ∫\_{|u|>ρc} e^{|σ|m(u) + ℓ₀u − Q\_s πu²/8} du ≤ 2e^{0.35|σ|} e^{aU − Q\_s πU²/8}/(Q\_s πU/4 − a), with U = c/4 and a = |σ|/c + |ℓ₀|. This requires Q\_s πU/4 > a; the exponent is concave with slope a − Q\_s πu/4 < 0 on u ≥ U.

The case q = 4 (Q\_s = 1) with even integer c ≥ 22 is the one written in full. **The general-q version above is the same computation with the Gaussian e^{−(π/q)(…)}, but it is written here only in this compact form (open item O4).** Its implementation is `dh_rect.gauss_fixed` and `dh_rect.lemI_tail_fixed` (fixed c), and `cert_extra.gauss_uniform` and `cert_extra.lemI_tail` (uniform in c ≥ c\_s).

A quadrature check (mpmath, q = 5, c ∈ {10.5, 15.5, 30.5}, δ ∈ {−½, 0, ½}, σ ∈ [−6, 6]) found J·c^σ below the bracket at every point tested.

**Step 5 (heat flow).**

- R\_{t,N}(σ + iT) := ∫ R\_{0,N}(σ + iT + √t v) e^{−v²}/√π dv.
- The contour shift v ↦ v + iπ√t/8 is valid: the integrand is holomorphic in the strip T ≤ Im s ≤ T + πt/8, since the poles of the Γ-factor lie on the real axis.
- It gives R\_{t,N}(σ + iT) = e^{π²t/64} ∫ e^{−iπ√t v/4} R\_{0,N}(σ + iT′ + √t v) e^{−v²}/√π dv, with T′ := T + πt/8. This uses Im α(iT′) = π/4 + O(1/T′).

On s\_v = σ\_v + iT′, with σ\_v := σ + √t v:

- |G(s\_v)| = |M\_0(iT′)| · w\_R^{σ\_v} · exp(O\_≤(1/(6(T′ − 0.66)) + |σ\_v|/T′² + σ\_v²/(4(T′ − 1))));
- Lemma I bounds J;
- (w\_R/c)^{σ\_v} = (1 + δ/c)^{σ\_v}.

Summing the two remainder integrals:

e\_C ≤ Σ\_{σ∈{(1−y)/2, (1+y)/2}} e^{π²t/64} (|M\_0(iT′)|/|B\_t|) K\_h ∫ (1 + δ/c)^{σ\_v} [G₀ + E + Tail](σ_v) e^{1/(6(T′−0.66)) + |σ\_v|/T′² + σ\_v²/(4(T′−1))} e^{−v²}/√π dv.

- The code evaluates the v-integral on 48 cells of |v| ≤ 12 with erf masses, plus a closed-form tail for |v| > 12 (`cert_barrier.tail_v`).
- Lemma T covers |σ\_v| > T′/4, where Lemma I is not used.
- The ratio |M\_0(iT′)|/|B\_t| is bounded by Lemma E(a).

## 3. Lemmas used by the certificates

**Lemma R′ (verified vector arithmetic; `rig.py`).**

- The model is IEEE-754 binary64 with round-to-nearest for +, −, ×, ÷ and √, plus exact rint and ldexp. No libm function is trusted.
- e^z is computed by Cody–Waite reduction with a 29-bit split of ln 2, then Horner on the degree-22 Taylor polynomial, then ldexp. The relative error is ≤ 110u + 1.02(e\_z + e\_r).
- e^{−iθ} is computed by Cody–Waite reduction with three 26-bit pieces of 2π, Taylor series on r/4, and two double-angle steps. The absolute error is ≤ e\_θ + e\_r + 200u.
- The constants come from Arb, with |split − true| ≤ 2^{−80}, checked at run time. log n comes from Arb.
- Every float summation of k terms carries Higham's γ\_k bound [Hi02, §4.2].

**Lemma H (the kernel bound K\_h; `cert_barrier.kh_rig`).** By the partial fractions of proof-rs-general.md (G1),

K(w)/e^{iπw²/q} = (1/2q) Σ\_{n mod 2q} a(n) e^{−iπn²/q}/(e^{iπ(w−n)/q} − 1).

For θ = a + ib, |e^{iθ} − 1|² = (1 − e^{−b})² + 4e^{−b} sin²(a/2). So on the line c + ue^{iπ/4},

|𝒦(w)|/|e^{iπw²/q}| ≤ |c₀|/2 + (1/2q) Σ\_{a(n)≠0} |a(n)| F\_n(u)^{−1/2},

with F\_n = (1 − e^{−b})² + 4e^{−b} sin²(a\_n/2), a\_n = π(c − n + u/√2)/q and b = πu/(√2q).

- The code uses the weights 1 ≥ |a(n)|.
- It takes the supremum over the 2q half-integer classes of c and over |u| ≤ 6. Each cell of width 10⁻³ is bounded by interval lower bounds of F; this is not a sampled grid.
- For |u| > 6 it uses #{a(n) ≠ 0}/(2q)/(1 − e^{−6π/(√2q)}).
- For q = 5 this gives K\_h ≤ 1.1794.

**Lemma H′ (analytic kernel bound, any q; `cert_barrier.kh_analytic`).** With b = πu/(√2q), β = |1 − e^{−b}| and θ\* = dist(½ + u/√2, ℤ),

(1/2q) Σ\_{n mod 2q} F\_n^{−1/2} ≤ min( e^{b/2}(H\_{q−1} + 1)/2 + min(1/(qβ), e^{b/2}/(4θ\*)), 1/β ).

**Proof.**

- Drop the condition a(n) ≠ 0.
- As n runs mod 2q, a\_n/2 mod π runs over π(m + θ)/(2q), with θ = frac(c + u/√2).
- sin x ≥ (2/π)min(x, π − x) gives sin(a\_n/2) ≥ d\_m/q, where {d\_m} = {j + θ, j + 1 − θ : 0 ≤ j < q}.
- G(d) = (β² + 4e^{−b}d²/q²)^{−1/2} is decreasing, so the sum is at most G(θ) + G(1 − θ) + q e^{b/2}H\_{q−1}.
- One of θ and 1 − θ is at least ½.
- F ≥ β² gives the bound 1/β.

The supremum over u is taken on cells of width 0.01 with outward interval bounds. For q = 5 this gives 2.9294, which the final region uses.

**Lemma M (monotonicity).** Let s₊ = (1 + y − ix)/2 and s₋ = (1 − y + ix)/2.

- **(a)** σ := Re s\* increases in x and in y, for x ≥ 10³:
  - ∂\_y σ = ½ + (t/4)Re α′(s₊) ≥ ½ − (t/4)(1/x + 1/x²) > 0;
  - ∂\_x σ = (t/4)Im α′(s₊), and Im(1/(2(s₊ + 𝔞))) = x/((1 + y + 2𝔞)² + x²) > 0, with the quadratic term O(1/x²).
- **(b)** ∂\_x log|γ| = −½[Im α\_t(s₋) + Im α\_t(s₊)], where α\_t = α + (t/2)αα′. To leading order this is −y/(2x) + πty/(8x²) + O(1/x²), which is negative for y ≥ 0.2 and x ≥ 300. **The constant in the O(1/x²) term is asserted (≤ 2), not derived, and the computation is written for y ≤ 1 (open item O3).**
- **(c)** ∂\_y log(|γ| m^y) = log(m/w\_R) + O(1/x) is negative for m ≤ N − 1. Since |∂²\_y log g\_m| ≤ 1/x, the certified blocks take every g\_m = |γ|m^y at both x-corners and both y-corners of a cell, and inflate the extremes by e^{Δy²/(4x\_lo)}. That covers any sign of ∂\_y.

**Lemma T (Gaussian tail in e\_C).** For |σ\_v| > T′/4 use crude global bounds instead of Lemma I:

- |G(s)| ≤ exp(|σ| log(|s| + q) + 1);
- |∫ w^{−s}𝒦 dw| ≤ K\_h·(c/√2)^{−|σ|}(1 + c)^{|σ|} e^{π|T|} (q/π)^{1/2} e^{πc²/(8q)·…}.

Against the Gaussian weight e^{−v²}, the contribution is ≤ exp(−T′²/(32t)) **for T′ ≥ 10⁶**. The certificates use T′ down to about 200, where the code bounds this region by 10⁻³⁰⁰. The expected size there is far smaller, but this has not been proved (open item O2).

**Lemma E (final-region error terms, uniform in x ≥ X; `dh_upper.final_region`).**

- **(a)** Write s₊ = −iT′ + Δ with Δ = (1 + y)/2 + iπt/8, and integrate α along the segment, using Re α ≥ log w\_R − 1/(2T′) and |Im α − π/4| ≤ 1/T′. This gives

  log(|M\_0(iT′)|/|M\_t(s₊)|) ≤ −((1 + y)/2)log w\_R(T′) + (1 + y)/(4T′) + (π/4 + 1/T′)πt/8 − (t/4)[(½log(qx/4π) − 1/x)² − (π/4 + 1/x)²].

  The total decreases in x and in y, since d/dy = −½log w\_R + 1/(4T′) < 0. So the corner (X, y₀) is the worst case.
- **(b)** In the variables 1/c ∈ [0, 1/c\_s], δ ∈ [−½, ½] and 1/T′ ∈ [0, 1/T′\_s], the Gaussian average of Lemma I's bracket is one Arb evaluation, uniform over all x ≥ X.
- **(c)** e\_A + e\_B ≤ ε(X)·F, with F ≤ (4π/t)^{1/2}(e^{((1−y₀)/2)²/t} + |γ(X, y₀)| e^{((1+Ytop)/2)²/t}) + 2. This is the sum-to-integral bound for Σ n^{−(1±y)/2}e^{−(t/4)log²n}, plus the maximal term.

**Lemma F (the tail N > N₀ of the final region, any y; `dh_upper.tail_any_y`).** For every N > N₀ = 10⁶ + 1 and y ∈ [y₁, y₂], the part of the plain bound coming from n ∈ (N₀, N] is at most

e^{−tL₀²/4 + c\_A L₀} + |γ(x(N₀), y₁)| e^{−tL₀²/4 + (c\_A + y₂)L₀},

where L₀ = log N₀, c\_A = (1 − y₁)/2 + tε/2 and ε = qt/(16N₀²) + 1/x(N₀), provided (t/2)L₀ > max(1, c\_A + y₂).

**Proof.**

- Re s\* ≥ (1 + y)/2 + (t/2)(log N − ε).
- Each sum over (N₀, N] is at most N·max of its summand at n = N₀ and at n = N, since the exponent is convex in log n.
- Both endpoint forms decrease in log N for log N ≥ L₀ under the stated condition, so the maximum is at N = N₀.
- The head (n ≤ N₀) is bounded by the last plain block, because Re s\* increases and |γ| decreases in x (Lemma M).

**Chebyshev certificate on an edge (`dh_rect.cheb_lb`).** Let f be analytic on the Bernstein ellipse E\_ρ, with |f| ≤ M there. Let p be its interpolant at the J + 1 Chebyshev points of the second kind, built from node values with errors ≤ e\_node. Then on [−1, 1]

|f − p| ≤ 4Mρ^{−J}/(ρ − 1) + Λ\_J e\_node + (float term),

by [Tr13, Thm 8.2], with Λ\_J ≤ (2/π)log(J + 1) + 1. The distance dist(p, (−∞, 0]) is bounded below on a grid in the angle, with step error ≤ (Σ\_k k|C\_k|)·Δθ/2.

**Direct evaluation of H\_t (`dh_rect.H_trap`).**

- The trapezoid rule on the Gaussian ξ-average with step h has error ≤ 2M/(e^{2πa/h} − 1), where M bounds ∫|F(w + ib)|dw for |b| < a [TW14, Thm 5.1].
- Here M ≤ e^{a²/(4t)}·sup over a strip of ∫(4πt)^{−1/2}e^{−w²/(4t)}|ξ| dw.
- The tail |w| > W uses `cert_dh.DH.tail_bound`.

## 4. Open items

- **O1.** Prop 6.1^χ (§2.3) is not written out; only the three changes to Polymath's proof are recorded.
- **O2.** Lemma T is proved only for T′ ≥ 10⁶; the certificates use T′ ≈ 200 and up.
- **O3.** Lemma M(b): the O(1/x²) constant is asserted, and the computation is written for y ≤ 1. The final region of papers/03-davenport-heilbronn-lambda uses y up to 1.485.
- **O4.** Lemma I is written in full only for q = 4 with even c ≥ 22. The general-q, half-integer-c form of §2.5 is the one the code uses.
- **O5.** For papers/02-lambda-unbounded-functional-equation-family (Lemma 1′), the qualitative version of §2 is needed uniformly for y up to about 2.36.

## References

- [P19] D.H.J. Polymath, "Effective approximation of heat flow evolution of the Riemann ξ function, and a new upper bound for the de Bruijn–Newman constant", arXiv:1904.12438.
- [dB50] N. G. de Bruijn, "The roots of trigonometric integrals", Duke Math. J. 17 (1950), 197–226.
- [Hi02] N. J. Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., SIAM, 2002.
- [Tr13] L. N. Trefethen, Approximation Theory and Approximation Practice, SIAM, 2013.
- [TW14] L. N. Trefethen and J. A. C. Weideman, "The exponentially convergent trapezoidal rule", SIAM Rev. 56 (2014), 385–458.
