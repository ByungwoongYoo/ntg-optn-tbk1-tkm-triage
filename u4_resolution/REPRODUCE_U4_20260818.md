# Ulrich u4 reproduction instructions

## Pinned formal-verification environment

- Operating system used in CI: Ubuntu 24.04
- Lean toolchain: `leanprover/lean4:v4.30.0`
- No Mathlib dependency
- Core arithmetic tactic module: `Lean.Elab.Tactic.Omega`

## 1. Install Lean

```bash
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh \
  | sh -s -- -y --default-toolchain leanprover/lean4:v4.30.0
export PATH="$HOME/.elan/bin:$PATH"
lean --version
```

## 2. Kernel-check the countermodel

From the repository root:

```bash
cd u4_resolution/lean
lean U4Core.lean 2>&1 | tee ../lean_check/U4Core_OUTPUT.txt
lean U4Formal.lean 2>&1 | tee ../lean_check/U4Formal_OUTPUT.txt
```

Both commands must return exit code `0`. `U4Formal.lean` ends with:

```lean
#check explicitU4Countermodel
#check u4_countermodel_exists
```

The second declaration must have the type:

```lean
∃ pred : Fml → Prop,
  (∀ x y z u, pred (Phi x y z u)) ∧
  (∀ X Y, pred X → pred (I X Y) → pred Y) ∧
  (∀ q, ¬ pred (I q q))
```

## 3. Run the independent symbolic checker

Python 3.12 was used, with no third-party packages.

```bash
cd u4_resolution
python3 u4_symbolic_checker.py \
  --family 250 \
  --enumerate 12 \
  --output U4_SYMBOLIC_CHECK.json \
  | tee U4_SYMBOLIC_CHECK.log
```

The process must return exit code `0` and the JSON field `all_checks_pass` must be `true`.

This checker is supporting evidence only. The all-formula proof object is the Lean file.

## 4. Audit the fixed TPTP input

```bash
cat u4.p
```

The u4 axiom must be:

```tptp
fof(u4,axiom,
    ! [X,Y,Z,U] :
      p(i(i(i(X,Y),Z), i(i(Y,i(Z,U)),i(Y,U)))) ).
```

This is exactly

```text
((X -> Y) -> Z) -> ((Y -> (Z -> U)) -> (Y -> U)).
```

The ATP timeout records are historical diagnostics and are not needed to validate the final semantic countermodel.

## 5. Verify hashes

```bash
find u4_resolution -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum > U4_SHA256SUMS.txt
sha256sum -c U4_SHA256SUMS.txt
```

The release artifact additionally contains the exact compiler outputs, exit codes, Lean version, symbolic-check JSON, and a manifest generated after all checks finish.
