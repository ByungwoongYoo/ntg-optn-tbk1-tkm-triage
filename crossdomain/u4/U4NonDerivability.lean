import Std.Tactic.Omega

/-!
A machine-checked syntactic refutation of Ulrich's u4 as a single axiom
for positive implicational logic.

The Hilbert calculus below has exactly the substitution instances of u4 as
axioms and modus ponens as its only inference rule.  We classify every
formula derivable in this calculus and prove that no reflexivity instance
P → P is derivable.
-/

inductive Ty where
  | var : Nat → Ty
  | arr : Ty → Ty → Ty
  deriving DecidableEq, Repr

namespace Ty

def size : Ty → Nat
  | .var _ => 1
  | .arr a b => a.size + b.size + 1

end Ty

infixr:60 " ⟶ " => Ty.arr

/-- Ulrich's u4 type/formula. -/
def U4 (a b c d : Ty) : Ty :=
  (((a ⟶ b) ⟶ c) ⟶ ((b ⟶ (c ⟶ d)) ⟶ (b ⟶ d)))

/-- The three-parameter shape occurring in every nontrivial right-comb iterate. -/
def KShape (x y z : Ty) : Ty :=
  ((x ⟶ (y ⟶ z)) ⟶ (x ⟶ z))

/--
`Pair n x y` is the recursively generated pair of side-types in the
principal type of the right comb with `n` copies of u4.
-/
inductive Pair : Nat → Ty → Ty → Prop where
  | base (x b c : Ty) : Pair 2 x (KShape b x c)
  | step {n : Nat} {x y : Ty} (h : Pair n x y) (r : Ty) :
      Pair (n + 1) (y ⟶ r) (x ⟶ r)

/-- Neither member of a generated pair can be the other followed by an arrow. -/
theorem pair_no_cross_arrow {n : Nat} {x y : Ty} (h : Pair n x y) :
    (∀ w : Ty, x ≠ (y ⟶ w)) ∧ (∀ w : Ty, y ≠ (x ⟶ w)) := by
  induction h with
  | base x b c =>
      constructor
      · intro w hEq
        have hs := congrArg Ty.size hEq
        simp [Ty.size, KShape] at hs
        omega
      · intro w hEq
        injection hEq with hDom _
        have hs := congrArg Ty.size hDom
        simp [Ty.size] at hs
        omega
  | @step n x y h r ih =>
      constructor
      · intro w hEq
        injection hEq with hLeft _
        exact (ih.2 r) hLeft
      · intro w hEq
        injection hEq with hLeft _
        exact (ih.1 r) hLeft

/--
A generated pair cannot simultaneously have the two shapes forced by
matching a nontrivial right-comb function domain with a u4 instance.
-/
theorem pair_no_u4_domain {n : Nat} {x y : Ty} (h : Pair n x y) :
    ∀ a b c d : Ty,
      ¬ (x = ((a ⟶ b) ⟶ c) ∧ y = (b ⟶ (c ⟶ d))) := by
  induction h with
  | base x p q =>
      intro a b c d hBoth
      rcases hBoth with ⟨hx, hy⟩
      injection hy with h1 h2
      injection h2 with hp hq
      have sx := congrArg Ty.size hx
      have sb := congrArg Ty.size h1
      have sp := congrArg Ty.size hp
      have sq := congrArg Ty.size hq
      simp [Ty.size] at sx sb sp sq
      omega
  | @step n x y h r ih =>
      intro a b c d hBoth
      rcases hBoth with ⟨hx, hy⟩
      injection hx with _ hr1
      injection hy with _ hr2
      have s1 := congrArg Ty.size hr1
      have s2 := congrArg Ty.size hr2
      simp [Ty.size] at s1 s2
      omega

/-- No substitution instance of u4 is an identity type. -/
theorem u4_not_identity (p a b c d : Ty) :
    (p ⟶ p) ≠ U4 a b c d := by
  intro hEq
  injection hEq with hDom hCod
  have hCore : ((a ⟶ b) ⟶ c) = KShape b c d := hDom.symm.trans hCod
  injection hCore with h1 h2
  injection h1 with _ hb
  have sb := congrArg Ty.size hb
  have sc := congrArg Ty.size h2
  simp [Ty.size] at sb sc
  omega

/--
`RType n t` describes all substitution instances of the principal type of
the unique typable right-comb derivation with `n` axiom occurrences.
-/
def RType (n : Nat) (t : Ty) : Prop :=
  (n = 1 ∧ ∃ a b c d : Ty, t = U4 a b c d) ∨
  (∃ x y z : Ty, Pair n x y ∧ t = KShape x y z)

/--
If one u4 instance is used as a function on an `RType n` argument, the
result is an `RType (n+1)` formula.
-/
theorem u4_apply_closure {n : Nat} {A B : Ty}
    (hU : ∃ a b c d : Ty, (A ⟶ B) = U4 a b c d)
    (hA : RType n A) : RType (n + 1) B := by
  rcases hU with ⟨a, b, c, d, hFun⟩
  rcases hA with hOne | hMany
  · rcases hOne with ⟨hn, p, q, r, s, hArg⟩
    subst n
    injection hFun with hDomain hResult
    have hMatch : U4 p q r s = ((a ⟶ b) ⟶ c) := hArg.symm.trans hDomain
    injection hMatch with hAB hC
    injection hAB with hAeq hBeq
    subst a
    subst b
    subst c
    refine Or.inr ?_
    refine ⟨r, KShape q r s, d, Pair.base r q s, ?_⟩
    exact hResult
  · rcases hMany with ⟨x, y, z, hPair, hArg⟩
    injection hFun with hDomain hResult
    have hMatch : KShape x y z = ((a ⟶ b) ⟶ c) := hArg.symm.trans hDomain
    injection hMatch with hAB hC
    injection hAB with hAeq hBeq
    subst a
    subst b
    subst c
    refine Or.inr ?_
    refine ⟨(y ⟶ z), (x ⟶ z), d, Pair.step hPair z, ?_⟩
    exact hResult

/-- A nontrivial right-comb type can never serve as a function on any RType argument. -/
theorem pair_domain_not_rtype {m n : Nat} {x y z : Ty}
    (hPair : Pair m x y) : ¬ RType n (x ⟶ (y ⟶ z)) := by
  intro hArg
  rcases hArg with hOne | hMany
  · rcases hOne with ⟨_, a, b, c, d, hEq⟩
    injection hEq with hx hyz
    injection hyz with hy _
    exact (pair_no_u4_domain hPair a b c d) ⟨hx, hy⟩
  · rcases hMany with ⟨x₂, y₂, z₂, hPair₂, hEq⟩
    injection hEq with hx hyz
    injection hyz with hy _
    have hBad : x = (y ⟶ (y₂ ⟶ z₂)) := by
      calc
        x = (x₂ ⟶ (y₂ ⟶ z₂)) := hx
        _ = (y ⟶ (y₂ ⟶ z₂)) := by rw [hy]
    exact (pair_no_cross_arrow hPair).1 (y₂ ⟶ z₂) hBad

/-- Hilbert derivability from substitution instances of u4 using only modus ponens. -/
inductive Derives : Ty → Prop where
  | ax (a b c d : Ty) : Derives (U4 a b c d)
  | mp {A B : Ty} : Derives (A ⟶ B) → Derives A → Derives B

/-- Every derivable formula belongs to one of the recursively classified RType families. -/
theorem derives_rtype {t : Ty} (h : Derives t) : ∃ n : Nat, RType n t := by
  induction h with
  | ax a b c d =>
      exact ⟨1, Or.inl ⟨rfl, a, b, c, d, rfl⟩⟩
  | @mp A B hFun hArg ihFun ihArg =>
      rcases ihFun with ⟨m, hm⟩
      rcases ihArg with ⟨n, hn⟩
      rcases hm with hOne | hMany
      · rcases hOne with ⟨_, a, b, c, d, hEq⟩
        exact ⟨n + 1, u4_apply_closure ⟨a, b, c, d, hEq⟩ hn⟩
      · rcases hMany with ⟨x, y, z, hPair, hEq⟩
        injection hEq with hDomain _
        have hnDomain : RType n (x ⟶ (y ⟶ z)) := by
          rw [← hDomain]
          exact hn
        exact False.elim ((pair_domain_not_rtype hPair) hnDomain)

/-- No classified type is a reflexivity instance. -/
theorem rtype_not_refl {n : Nat} {p : Ty} : ¬ RType n (p ⟶ p) := by
  intro h
  rcases h with hOne | hMany
  · rcases hOne with ⟨_, a, b, c, d, hEq⟩
    exact u4_not_identity p a b c d hEq
  · rcases hMany with ⟨x, y, z, hPair, hEq⟩
    injection hEq with hDomain hCodomain
    have hCore : (x ⟶ (y ⟶ z)) = (x ⟶ z) := hDomain.symm.trans hCodomain
    injection hCore with _ hYZ
    have hs := congrArg Ty.size hYZ
    simp [Ty.size] at hs
    omega

/-- Main result: u4 does not derive even P → P. -/
theorem u4_does_not_derive_reflexivity (p : Ty) : ¬ Derives (p ⟶ p) := by
  intro h
  rcases derives_rtype h with ⟨n, hn⟩
  exact rtype_not_refl hn

#print axioms u4_does_not_derive_reflexivity
#check u4_does_not_derive_reflexivity
