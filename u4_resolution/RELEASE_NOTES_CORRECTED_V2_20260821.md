# Corrected v2 release notes — 21 August 2026 KST

This is a new version under the existing Zenodo concept DOI. It does not delete
or rewrite the public v1 record dated 18 August 2026.

## Corrections and strengthening relative to v1

- Rebuilt the preprint so the complete `Recursive schemes` section and the
  definitions of `C`, `Phi`, and `T_1` are present in the PDF.
- Upgraded the formal verification target from the historical Lean 4.30.0
  snapshot to Lean 4.33.0.
- Added explicit uniform substitution, an explicit `U4Derivable` system,
  soundness, and a checked separation from the stated `K`/`S` basis.
- Added a source scan, direct compilation, `leanchecker`, the separately
  implemented Rust checker `nanoda`, and a compiled axiom-dependency audit.
- Expanded the separate Python checks to 100 recursive schemes and exhaustive
  enumeration of all 58,786 binary application-tree shapes at 12 leaves.
- Replaced the non-portable historical checksum layout with archive-relative
  SHA-256 entries that exclude the checksum file itself.
- Added explicit licensing, citation metadata, AI-use disclosure, release
  notes, and a machine-readable manifest.

## Time record

- Historical v1 publication: 18 August 2026.
- Baseline proof-source commit: 20 August 2026 15:10:53 UTC, equivalent to
  21 August 2026 00:10:53 KST.
- Corrected v2 publication date: 21 August 2026 KST.

## Claim boundary

The result is machine checked and supported by multiple project-run checkers.
Independent third-party specialist review, journal peer review, and comparative
worldwide-priority review remain pending.
