import Std.Tactic.Omega

/- First kernel-checked fragment of the explicit infinite-model proof. -/

inductive Form where
  | v : Nat → Form
  | arr : Form → Form → Form
  deriving DecidableEq, Repr

namespace Form

@[simp] def subst (σ : Nat → Form) : Form → Form
  | .v n => σ n
  | .arr a b => .arr (subst σ a) (subst σ b)

@[simp] def size : Form → Nat
  | .v _ => 1
  | .arr a b => size a + size b + 1

@[simp] theorem size_pos (f : Form) : 0 < size f := by
  induction f with
  | v n => simp
  | arr a b iha ihb => simp [size]; omega

end Form

open Form

abbrev V (n : Nat) : Form := .v n
abbrev I (a b : Form) : Form := .arr a b

def C (a b r : Form) : Form := I (I a (I b r)) (I a r)

def phi : Form := I (I (I (V 0) (V 1)) (V 2)) (C (V 1) (V 2) (V 3))

def rv (n : Nat) : Nat := 100 + n

def R (n : Nat) : Form := V (rv n)

mutual
  def A : Nat → Form
    | 0 => V 4
    | n+1 => I (B n) (R n)
  def B : Nat → Form
    | 0 => C (V 5) (V 4) (V 6)
    | n+1 => I (A n) (R n)
end

def D (n : Nat) : Form := I (A n) (I (B n) (R n))
def Td (n : Nat) : Form := C (A n) (B n) (R n)
def T : Nat → Form
  | 0 => phi
  | n+1 => Td n

def E (P Q : Form) : Prop := ∀ (σ : Nat → Form) (W : Form), subst σ P ≠ I (subst σ Q) W

lemma E_base_AB : E (A 0) (B 0) := by
  intro σ W h
  have hs := congrArg Form.size h
  simp [A, B, C, Form.subst, Form.size] at hs
  have ha := Form.size_pos (σ 4)
  have hb := Form.size_pos (σ 5)
  have hc := Form.size_pos (σ 6)
  have hw := Form.size_pos W
  omega

lemma E_base_BA : E (B 0) (A 0) := by
  intro σ W h
  have hs := congrArg Form.size h
  simp [A, B, C, Form.subst, Form.size] at hs
  have ha := Form.size_pos (σ 4)
  have hb := Form.size_pos (σ 5)
  have hc := Form.size_pos (σ 6)
  have hw := Form.size_pos W
  omega

lemma interlocking : ∀ n, E (A n) (B n) ∧ E (B n) (A n)
  | 0 => ⟨E_base_AB, E_base_BA⟩
  | n+1 => by
      have ih := interlocking n
      constructor
      · intro σ W h
        simp [A, B, Form.subst] at h
        injection h with hleft hright
        exact ih.2 σ (σ (rv n)) hleft
      · intro σ W h
        simp [A, B, Form.subst] at h
        injection h with hleft hright
        exact ih.1 σ (σ (rv n)) hleft

example : ∀ n, E (A n) (B n) := fun n => (interlocking n).1
