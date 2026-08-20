import Lean.Elab.Tactic.Omega

namespace U4Formal

/-!
Kernel-checked explicit infinite model for Ulrich's u4.
The syntax is the free algebra of finite implicational formulas.
-/

inductive Sym where
  | x | y | z | u | a | b | c
  | r : Nat → Sym
  deriving DecidableEq, Repr

inductive Fml where
  | v : Sym → Fml
  | arr : Fml → Fml → Fml
  deriving DecidableEq, Repr

namespace Fml

@[simp] def size : Fml → Nat
  | .v _ => 1
  | .arr p q => size p + size q + 1

@[simp] theorem size_pos (f : Fml) : 0 < size f := by
  induction f with
  | v s => simp
  | arr p q ihp ihq => simp [size]

end Fml

open Fml

abbrev I (p q : Fml) : Fml := .arr p q

def C (p q r : Fml) : Fml := I (I p (I q r)) (I p r)

def Phi (x y z u : Fml) : Fml := I (I (I x y) z) (C y z u)

theorem arr_inj {a b c d : Fml} (h : I a b = I c d) : a = c ∧ b = d := by
  cases h
  exact ⟨rfl, rfl⟩

theorem no_self_left (p q : Fml) : p ≠ I p q := by
  intro h
  have hs := congrArg Fml.size h
  have hq := Fml.size_pos q
  simp [Fml.size] at hs
  omega

theorem no_self_right (p q : Fml) : p ≠ I q p := by
  intro h
  have hs := congrArg Fml.size h
  have hq := Fml.size_pos q
  simp [Fml.size] at hs
  omega

/- Pair a b says that (a,b) is one of the recursive interlocking pairs,
   with arbitrary finite formulas substituted for the scheme variables. -/
inductive Pair : Fml → Fml → Prop where
  | base (a b c : Fml) : Pair a (C b a c)
  | step {a b : Fml} : Pair a b → (r : Fml) → Pair (I b r) (I a r)

namespace Pair

theorem interlocking {a b : Fml} (h : Pair a b) :
    (∀ w, a ≠ I b w) ∧ (∀ w, b ≠ I a w) := by
  induction h with
  | base a b c =>
      constructor
      · intro w eq
        have hs := congrArg Fml.size eq
        have ha := Fml.size_pos a
        have hb := Fml.size_pos b
        have hc := Fml.size_pos c
        have hw := Fml.size_pos w
        simp [C, Fml.size] at hs
        omega
      · intro w eq
        have leftEq := (arr_inj eq).1
        have hs := congrArg Fml.size leftEq
        have ha := Fml.size_pos a
        have hb := Fml.size_pos b
        have hc := Fml.size_pos c
        simp [Fml.size] at hs
        omega
  | step h r ih =>
      constructor
      · intro w eq
        have leftEq := (arr_inj eq).1
        exact ih.2 r leftEq
      · intro w eq
        have leftEq := (arr_inj eq).1
        exact ih.1 r leftEq

theorem no_D_eq_C {a b c d : Fml} (hab : Pair a b) (_hcd : Pair c d)
    (r s : Fml) : I a (I b r) ≠ C c d s := by
  intro eq
  have outer := arr_inj eq
  have inner := arr_inj outer.2
  have ha := outer.1
  have hbc := inner.1
  rw [← hbc] at ha
  exact (interlocking hab).1 (I d s) ha

theorem no_D_eq_Phi {a b : Fml} (hab : Pair a b) (r x y z u : Fml) :
    I a (I b r) ≠ Phi x y z u := by
  induction hab with
  | base a b c =>
      intro eq
      have outer := arr_inj eq
      have rest := arr_inj outer.2
      have cparts := arr_inj rest.1
      have bczu := arr_inj cparts.2
      have hsa := congrArg Fml.size outer.1
      have hsy := congrArg Fml.size cparts.1
      have hsz := congrArg Fml.size bczu.1
      have hx := Fml.size_pos x
      have ha := Fml.size_pos a
      have hb := Fml.size_pos b
      have hc := Fml.size_pos c
      have hu := Fml.size_pos u
      simp [Fml.size] at hsa hsy hsz
      omega
  | step hab s ih =>
      intro eq
      have outer := arr_inj eq
      have ls := (arr_inj outer.1).2
      have rs := (arr_inj (arr_inj outer.2).1).2
      have bad : z = I z u := ls.symm.trans rs
      exact no_self_left z u bad

end Pair

theorem Phi_no_refl (x y z u q : Fml) : Phi x y z u ≠ I q q := by
  intro eq
  have outer := arr_inj eq
  have both : I (I x y) z = C y z u := outer.1.trans outer.2.symm
  have parts := arr_inj both
  have hy := (arr_inj parts.1).2
  have hz := parts.2
  have hsy := congrArg Fml.size hy
  have hsz := congrArg Fml.size hz
  have hu := Fml.size_pos u
  simp [Fml.size] at hsy hsz
  omega

theorem C_no_refl (a b r q : Fml) : C a b r ≠ I q q := by
  intro eq
  have outer := arr_inj eq
  have both : I a (I b r) = I a r := outer.1.trans outer.2.symm
  have bad : I b r = r := (arr_inj both).2
  exact no_self_right r b bad.symm

/- The theorem predicate of the explicit infinite model. -/
def P (f : Fml) : Prop :=
  (∃ x y z u, Phi x y z u = f) ∨
  (∃ a b r, Pair a b ∧ C a b r = f)

theorem P_phi (x y z u : Fml) : P (Phi x y z u) := by
  exact Or.inl ⟨x, y, z, u, rfl⟩

theorem P_C {a b : Fml} (h : Pair a b) (r : Fml) : P (C a b r) := by
  exact Or.inr ⟨a, b, r, h, rfl⟩

theorem P_mp {X Y : Fml} (hx : P X) (hxy : P (I X Y)) : P Y := by
  rcases hxy with hmaj | hmaj
  · rcases hmaj with ⟨x, y, z, u, hmaj⟩
    have majorParts := arr_inj hmaj
    rcases hx with hmin | hmin
    · rcases hmin with ⟨x', y', z', u', hmin⟩
      have common : Phi x' y' z' u' = I (I x y) z := hmin.trans majorParts.1.symm
      have commonParts := arr_inj common
      have anteParts := arr_inj commonParts.1
      have hp : P (C z' (C y' z' u') u) := P_C (Pair.base z' y' u') u
      rw [commonParts.2, anteParts.2, majorParts.2] at hp
      exact hp
    · rcases hmin with ⟨a, b, r, hab, hmin⟩
      have common : C a b r = I (I x y) z := hmin.trans majorParts.1.symm
      have commonParts := arr_inj common
      have dParts := arr_inj commonParts.1
      have hp : P (C (I b r) (I a r) u) := P_C (Pair.step hab r) u
      rw [dParts.2, commonParts.2, majorParts.2] at hp
      exact hp
  · rcases hmaj with ⟨a, b, r, hab, hmaj⟩
    have majorParts := arr_inj hmaj
    rcases hx with hmin | hmin
    · rcases hmin with ⟨x, y, z, u, hmin⟩
      have bad : I a (I b r) = Phi x y z u := majorParts.1.trans hmin.symm
      exact False.elim (Pair.no_D_eq_Phi hab r x y z u bad)
    · rcases hmin with ⟨c, d, s, hcd, hmin⟩
      have bad : I a (I b r) = C c d s := majorParts.1.trans hmin.symm
      exact False.elim (Pair.no_D_eq_C hab hcd r s bad)

theorem P_no_refl (q : Fml) : ¬ P (I q q) := by
  intro h
  rcases h with h | h
  · rcases h with ⟨x, y, z, u, eq⟩
    exact Phi_no_refl x y z u q eq
  · rcases h with ⟨a, b, r, hab, eq⟩
    exact C_no_refl a b r q eq

structure U4Countermodel where
  pred : Fml → Prop
  u4_true : ∀ x y z u, pred (Phi x y z u)
  mp_closed : ∀ {X Y}, pred X → pred (I X Y) → pred Y
  refl_false : ∀ q, ¬ pred (I q q)

def explicitU4Countermodel : U4Countermodel where
  pred := P
  u4_true := P_phi
  mp_closed := P_mp
  refl_false := P_no_refl

theorem u4_countermodel_exists :
    ∃ pred : Fml → Prop,
      (∀ x y z u, pred (Phi x y z u)) ∧
      (∀ X Y, pred X → pred (I X Y) → pred Y) ∧
      (∀ q, ¬ pred (I q q)) := by
  exact ⟨P, P_phi, fun X Y => P_mp, P_no_refl⟩

/-!
## Direct correspondence with the original Hilbert-style problem

The countermodel theorem above is already enough semantically, but the following
definitions close the presentation gap that a reviewer could reasonably question:

* substitutions are explicit endomorphisms of the free formula algebra;
* the designated predicate is proved closed under every substitution;
* derivability contains exactly arbitrary `u4` instances, modus ponens, and an
  explicit substitution rule; and
* every derivable formula is designated in the countermodel.

Consequently no reflexive formula can be derived in the original system.
-/

abbrev Substitution := Sym → Fml

def applySubst (s : Substitution) : Fml → Fml
  | .v sym => s sym
  | .arr p q => I (applySubst s p) (applySubst s q)

@[simp] theorem applySubst_arr (s : Substitution) (p q : Fml) :
    applySubst s (I p q) = I (applySubst s p) (applySubst s q) := rfl

@[simp] theorem applySubst_C (s : Substitution) (p q r : Fml) :
    applySubst s (C p q r) =
      C (applySubst s p) (applySubst s q) (applySubst s r) := rfl

@[simp] theorem applySubst_Phi (s : Substitution) (x y z u : Fml) :
    applySubst s (Phi x y z u) =
      Phi (applySubst s x) (applySubst s y) (applySubst s z) (applySubst s u) := rfl

namespace Pair

theorem subst_closed {a b : Fml} (h : Pair a b) (s : Substitution) :
    Pair (applySubst s a) (applySubst s b) := by
  induction h with
  | base a b c =>
      simpa using Pair.base (applySubst s a) (applySubst s b) (applySubst s c)
  | step h r ih =>
      simpa using Pair.step ih (applySubst s r)

end Pair

theorem P_subst_closed {f : Fml} (hf : P f) (s : Substitution) :
    P (applySubst s f) := by
  rcases hf with hphi | hc
  · rcases hphi with ⟨x, y, z, u, rfl⟩
    simpa using P_phi (applySubst s x) (applySubst s y)
      (applySubst s z) (applySubst s u)
  · rcases hc with ⟨a, b, r, hab, rfl⟩
    simpa using P_C (Pair.subst_closed hab s) (applySubst s r)

inductive U4Derivable : Fml → Prop where
  | u4 (x y z u : Fml) : U4Derivable (Phi x y z u)
  | mp {X Y : Fml} : U4Derivable X → U4Derivable (I X Y) → U4Derivable Y
  | subst {f : Fml} : U4Derivable f → (s : Substitution) →
      U4Derivable (applySubst s f)

theorem U4Derivable.sound {f : Fml} (h : U4Derivable f) : P f := by
  induction h with
  | u4 x y z u => exact P_phi x y z u
  | mp hX hXY ihX ihXY => exact P_mp ihX ihXY
  | subst hf s ih => exact P_subst_closed ih s

theorem u4_not_derives_reflexivity (q : Fml) : ¬ U4Derivable (I q q) := by
  intro h
  exact P_no_refl q h.sound

theorem u4_fails_single_axiom_test : ∀ q : Fml, ¬ U4Derivable (I q q) := by
  exact u4_not_derives_reflexivity

/-!
## Separation from positive implicational logic

The countermodel already refutes the candidate, because reflexivity is a
theorem of positive implicational logic.  For a fully internalized comparison,
we also define the standard `K`/`S` Hilbert basis and derive reflexivity from it.
Thus the final separation theorem has no informal premise about the target
logic left outside the checked file.
-/

def K (p q : Fml) : Fml := I p (I q p)

def S (p q r : Fml) : Fml :=
  I (I p (I q r)) (I (I p q) (I p r))

inductive PositiveDerivable : Fml → Prop where
  | k (p q : Fml) : PositiveDerivable (K p q)
  | s (p q r : Fml) : PositiveDerivable (S p q r)
  | mp {X Y : Fml} :
      PositiveDerivable X → PositiveDerivable (I X Y) → PositiveDerivable Y
  | subst {f : Fml} : PositiveDerivable f → (s : Substitution) →
      PositiveDerivable (applySubst s f)

theorem positive_derives_reflexivity (q : Fml) :
    PositiveDerivable (I q q) := by
  have hs : PositiveDerivable (S q (I q q) q) :=
    PositiveDerivable.s q (I q q) q
  have hk₁ : PositiveDerivable (K q (I q q)) :=
    PositiveDerivable.k q (I q q)
  have hk₂ : PositiveDerivable (K q q) :=
    PositiveDerivable.k q q
  have hstep : PositiveDerivable (I (K q q) (I q q)) :=
    PositiveDerivable.mp hk₁ hs
  exact PositiveDerivable.mp hk₂ hstep

theorem positive_reflexivity_not_u4 (q : Fml) :
    PositiveDerivable (I q q) ∧ ¬ U4Derivable (I q q) := by
  exact ⟨positive_derives_reflexivity q, u4_not_derives_reflexivity q⟩

theorem u4_not_axiomatizes_positive_implication :
    ∃ f : Fml, PositiveDerivable f ∧ ¬ U4Derivable f := by
  exact ⟨I (Fml.v Sym.x) (Fml.v Sym.x),
    positive_reflexivity_not_u4 (Fml.v Sym.x)⟩

#check explicitU4Countermodel
#check u4_countermodel_exists
#check P_subst_closed
#check U4Derivable.sound
#check u4_not_derives_reflexivity
#check u4_fails_single_axiom_test
#check positive_derives_reflexivity
#check positive_reflexivity_not_u4
#check u4_not_axiomatizes_positive_implication

#print axioms u4_countermodel_exists
#print axioms U4Derivable.sound
#print axioms u4_not_axiomatizes_positive_implication

end U4Formal
