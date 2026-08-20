# Reproducing the strengthened Ulrich u4 verification

## Pinned environment

- Lean toolchain: `leanprover/lean4:v4.33.0`
- Lake package: `u4_resolution/lean/lakefile.lean`
- Mathlib: not used
- CI operating system: Ubuntu 24.04

The dated Lean 4.30 logs in `lean_check/` are retained as historical evidence
for the earlier proof. They do not certify the strengthened source. The current
source is certified only by a green run of
`.github/workflows/u4-strengthened-independent-check-20260820.yml` at the same
Git commit.

## 1. Build and check with Lean

Install `elan`, then from the repository root run:

```bash
cd u4_resolution/lean
elan toolchain install leanprover/lean4:v4.33.0
lake build --wfail
lake env lean U4Formal.lean
```

The final command prints the types and axiom dependencies of the principal
theorems. It must exit with status 0.

## 2. Run the independent symbolic checker

```bash
cd u4_resolution
python3 u4_symbolic_checker.py \
  --family 100 \
  --enumerate 12 \
  --output U4_SYMBOLIC_CHECK.json
```

The process must exit with status 0 and the JSON field `all_checks_pass` must be
`true`. This finite, independently implemented check supports but does not
replace the all-formula Lean proof.

## 3. Run the independent kernel checks

The repository workflow invokes `leanprover/lean-action@v1` with all of the
following features explicitly enabled:

```yaml
leanchecker: true
nanoda: true
nanoda-allow-sorry: false
axiom-audit: true
axiom-audit-allow: 'propext,Classical.choice,Quot.sound'
axiom-audit-root: U4Formal
```

`leanchecker` rechecks the compiled environment. `nanoda` is a separate Lean 4
type checker written in Rust. The axiom audit inspects compiled declarations,
not merely source text. The job fails unless the build, direct compilation,
both independent checks, the axiom audit, and the symbolic checker all pass.

## 4. Bind results to a commit

For a successful workflow run, download the artifact named
`u4-strengthened-independent-verification-20260820`. Confirm that:

1. `VERIFICATION_STATUS.txt` names the intended Git commit;
2. every recorded status is successful;
3. `LEAN_DIRECT_EXIT.txt` is `0`;
4. the symbolic-check JSON reports `all_checks_pass: true`; and
5. the included source files match `SHA256SUMS.txt`.

Do not create or update a release, DOI, or preprint package from a commit that
fails any gate or from an artifact produced for a different commit.
