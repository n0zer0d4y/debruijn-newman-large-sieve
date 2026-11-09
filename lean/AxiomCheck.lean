/-
  Axiom check: every theorem of the library (main theorems and witnesses) should depend only on `propext`, `Classical.choice` and `Quot.sound`
  (no `sorryAx`, no `Lean.ofReduceBool`). Run from `lean/`:  lake env lean AxiomCheck.lean
-/
import DebruijnNewmanLargeSieve

open DebruijnNewmanLargeSieve.PrimeSpoke

#print axioms spoke_point
#print axioms spoke_bounds
#print axioms exists_rich_class
#print axioms spoke_hyps_satisfiable
#print axioms exists_rich_class_rat
#print axioms spoke_example
#print axioms rich_class_example
