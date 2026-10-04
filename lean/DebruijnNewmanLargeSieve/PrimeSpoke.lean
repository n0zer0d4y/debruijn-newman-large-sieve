import Mathlib

/-!
# The prime-denominator spoke (papers/01-prime-moduli-large-sieve)

This file formalises the finite core of the lower bound for the prime-moduli large sieve:

* `spoke_point`: if `k ∣ 1 + b * p` and `a = (1 + b * p) / k`, then `a / p = b / k + 1 / (k * p)`
  in `ℚ`. So the Farey point `a / p` sits on the spoke of `b / k` at distance `1 / (k * p)`.
* `spoke_bounds`: with `b < k` and `2 ≤ p`, the numerator satisfies `1 ≤ a < p`, so `a / p` is one
  of the points `a / p`, `1 ≤ a ≤ p - 1`, of the large sieve with modulus `p`.
* `exists_rich_class`: pigeonhole. Some reduced residue class mod `k` contains at least
  `#s / φ k` elements of any finite set `s` of integers coprime to `k`.

The analytic part (Jordan's inequality, the prime number theorem, Mertens) is in the written proof.
-/

namespace DebruijnNewmanLargeSieve.PrimeSpoke

open Finset

/-- The spoke identity: `a / p = b / k + 1 / (k p)` for `a = (1 + b p) / k`. -/
theorem spoke_point (k p b : ℕ) (hk : 0 < k) (hp : 0 < p) (hdiv : k ∣ 1 + b * p) :
    (((1 + b * p) / k : ℕ) : ℚ) / p = (b : ℚ) / k + 1 / ((k : ℚ) * p) := by
  obtain ⟨a, ha⟩ := hdiv
  have hk' : (k : ℚ) ≠ 0 := by exact_mod_cast hk.ne'
  have hp' : (p : ℚ) ≠ 0 := by exact_mod_cast hp.ne'
  have hq : (1 + b * p) / k = a := by
    rw [ha]; exact Nat.mul_div_cancel_left a hk
  rw [hq]
  have hcast : ((k : ℚ) * a) = 1 + (b : ℚ) * p := by exact_mod_cast ha.symm
  field_simp
  linear_combination hcast

/-- The spoke numerator is a genuine Farey numerator: `1 ≤ a < p`. -/
theorem spoke_bounds (k p b : ℕ) (hk : 0 < k) (hp : 2 ≤ p) (hb : b < k)
    (hdiv : k ∣ 1 + b * p) :
    1 ≤ (1 + b * p) / k ∧ (1 + b * p) / k < p := by
  obtain ⟨a, ha⟩ := hdiv
  have hq : (1 + b * p) / k = a := by
    rw [ha]; exact Nat.mul_div_cancel_left a hk
  rw [hq]
  constructor
  · rcases Nat.eq_zero_or_pos a with h0 | h0
    · rw [h0] at ha; omega
    · exact h0
  · -- k * a = 1 + b p ≤ 1 + (k - 1) p < k p
    by_contra hge
    push_neg at hge
    have h1 : k * p ≤ k * a := Nat.mul_le_mul_left k hge
    have h2 : b * p ≤ (k - 1) * p := Nat.mul_le_mul_right p (by omega)
    have h3 : (k - 1) * p + p = k * p := by
      rcases k with _ | k
      · omega
      · simp [Nat.succ_mul]
    omega

/-- Pigeonhole over reduced residue classes: for any finite set `s` of naturals coprime to `k`,
some reduced class `c < k` contains at least `#s / φ k` elements of `s`. -/
theorem exists_rich_class (k : ℕ) (hk : 1 < k) (s : Finset ℕ)
    (hs : ∀ p ∈ s, Nat.Coprime k p) :
    ∃ c ∈ (range k).filter (fun a => Nat.Coprime k a),
      #s / Nat.totient k ≤ #(s.filter (fun p => p % k = c)) := by
  set t := (range k).filter (fun a => Nat.Coprime k a) with ht
  have hmaps : ∀ p ∈ s, p % k ∈ t := by
    intro p hp
    rw [ht, mem_filter, mem_range]
    refine ⟨Nat.mod_lt p (by omega), ?_⟩
    have h := hs p hp
    unfold Nat.Coprime at h ⊢
    rw [Nat.gcd_comm, ← Nat.gcd_rec]
    exact h
  have hne : t.Nonempty := ⟨1, by rw [ht, mem_filter, mem_range]; exact ⟨hk, Nat.coprime_one_right k⟩⟩
  have hcard : #t = Nat.totient k := (Nat.totient_eq_card_coprime k).symm
  have hn : #t * (#s / Nat.totient k) ≤ #s := by
    rw [hcard]; exact Nat.mul_div_le (#s) (Nat.totient k)
  exact exists_le_card_fiber_of_mul_le_card_of_maps_to hmaps hne hn

/-- Witness: the three hypotheses of `spoke_point` and `spoke_bounds` hold together
(`k = 6`, `p = 5`, `b = 1`, since `6 ∣ 1 + 1·5`), with numerator `a = 1`. -/
theorem spoke_hyps_satisfiable :
    (0 < 6 ∧ 2 ≤ 5 ∧ 1 < 6 ∧ 6 ∣ 1 + 1 * 5) ∧ (1 + 1 * 5) / 6 = 1 := by
  decide

end DebruijnNewmanLargeSieve.PrimeSpoke
