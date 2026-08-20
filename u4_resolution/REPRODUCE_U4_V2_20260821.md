# Reproducing the corrected v2 Ulrich u4 verification

## Pinned environment

- Lean toolchain: `leanprover/lean4:v4.33.0`
- Lake package: `u4_resolution/lean/lakefile.lean`
- Mathlib: not used
- CI operating system: Ubuntu 24.04

The Lean 4.30 logs dated 2026-08-18 in `lean_check/` are retained as historical
evidence for the earlier proof. They do not certify the corrected v2 source.
The current source and release files are certified only by an all-green run of
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

The repository workflow pins the official Lean action commit
`6835ef47e8423fa9e9eacfcbb8fdf83fae42c820` and enables the following checks:

```yaml
leanchecker: true
leanchecker-args: U4Formal
axiom-audit: true
axiom-audit-allow: 'propext,Classical.choice,Quot.sound'
axiom-audit-root: U4Formal
```

`leanchecker` rechecks the compiled environment. The axiom audit inspects
compiled declarations, not merely source text.

The workflow separately pins the Lean 4 exporter and the Rust implementation
of `nanoda` to the following exact commits:

- `lean4export`: `15f6055e299ad5b89345e533cc2192f4cc00f659`
  (the Lean 4.33.0 tag);
- `nanoda_lib`: `6ae1f0cd962f081f6c423454c5da729d841236a7`.

The `nanoda` allowlist excludes `sorryAx`. Unpermitted axioms are ignored only
when unused and cause a hard error if a checked declaration depends on them.
The compiled axiom audit independently enforces the stated allowlist. The job
fails unless the build, direct compilation, `leanchecker`, `nanoda`, the axiom
audit, and the symbolic checker all pass.

The baseline all-green strengthened run is
<https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32384585278>
at proof-source commit `fec164ed763794330d57f76ddd1b0c390f22db09`
(2026-08-20 15:10:53 UTC; 2026-08-21 00:10:53 KST). Its `nanoda` log reports
96,602 declarations with no typechecker errors; the single reported diagnostic
concerns pretty-printing axioms, not type checking. A corrected-v2 release must
also have an all-green run at its own exact release commit.

## 4. Bind results to a commit

For a successful workflow run, download the artifact named
`u4-strengthened-v2-20260821`. Confirm that:

1. `VERIFICATION_STATUS.txt` names the intended Git commit;
2. every recorded status is successful;
3. `LEAN_DIRECT_EXIT.txt` is `0`;
4. the symbolic-check JSON reports `all_checks_pass: true`; and
5. from the extracted artifact directory, `sha256sum -c SHA256SUMS.txt`
   succeeds for every listed file.

Do not create or update a release, DOI, or preprint package from a commit that
fails any gate or from an artifact produced for a different commit.

## 5. Public identifiers

- Corrected v2 DOI: <https://doi.org/10.5281/zenodo.22031656>
- Concept DOI: <https://doi.org/10.5281/zenodo.21987698>
- Historical v1 DOI: <https://doi.org/10.5281/zenodo.21987699>

The version DOI identifies corrected v2; the concept DOI identifies its version family.
