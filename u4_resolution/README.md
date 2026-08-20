# Ulrich u4 corrected v2 release

This directory contains a formal and executable certificate for the following
claim. Let

```text
u4 = ((x -> y) -> z) -> ((y -> (z -> u)) -> (y -> u)).
```

The displayed formula is not a single axiom for positive implicational logic
under uniform substitution and modus ponens.

## What is proved

The construction uses the free algebra of finite implicational formulas and an
explicit substitution-closed predicate `P`. The Lean development proves:

1. every substitution instance of `u4` is in `P`;
2. `P` is closed under modus ponens;
3. no reflexive formula `q -> q` is in `P`;
4. every derivation from `u4` using substitution and modus ponens is sound in
   `P`; and
5. reflexivity is derivable from the explicit standard `K`/`S` basis but not
   from `u4`.

The primary certificate is `lean/U4Formal.lean`. The human-readable proof is
in `preprint/U4_PREPRINT_CORRECTED_V2_20260821.pdf`.

## Verification

The read-only GitHub workflow checks the exact release commit with:

- Lean 4.33.0 build with warnings treated as errors;
- direct Lean compilation;
- `leanchecker` environment rechecking;
- the separately implemented Rust checker `nanoda`;
- a compiled axiom-dependency audit; and
- a separate Python symbolic checker through 100 recursive schemes and every
  binary application-tree shape through 12 leaves.

The historical proof-source run is [GitHub Actions run
32384585278](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32384585278)
at commit `fec164ed763794330d57f76ddd1b0c390f22db09`. The corrected v2 package must
also be paired with an all-green artifact from its own exact release commit.

## Public identifiers

- Corrected v2 DOI:
  [10.5281/zenodo.22031656](https://doi.org/10.5281/zenodo.22031656)
- Concept DOI:
  [10.5281/zenodo.21987698](https://doi.org/10.5281/zenodo.21987698)
- Historical v1 DOI:
  [10.5281/zenodo.21987699](https://doi.org/10.5281/zenodo.21987699)

## Claim boundary

These are machine-verification results produced within this project. No
independent third-party specialist review or journal peer review has yet
occurred. The release does not claim absolute worldwide priority over every
unpublished or non-indexed result.

## Licensing

Code is released under the MIT License. The preprint, documentation, logs,
tables, and metadata are released under CC BY 4.0. See `LICENSE_SCOPE.md`.
