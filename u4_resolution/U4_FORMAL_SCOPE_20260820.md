# Ulrich u4: formal scope and verification map

## Result represented by the formalization

Let

```text
Phi(x,y,z,u) = ((x -> y) -> z) -> ((y -> (z -> u)) -> (y -> u)).
```

The formalization constructs an explicit predicate `P` on the free algebra of
finite implicational formulas and proves all three clauses below.

| Mathematical obligation | Lean declaration |
|---|---|
| Every substitution instance of `u4` is designated | `U4Formal.P_phi` |
| Designated formulas are closed under modus ponens | `U4Formal.P_mp` |
| No reflexive formula `q -> q` is designated | `U4Formal.P_no_refl` |

Together these are packaged by `U4Formal.u4_countermodel_exists`. The
formalized formula agrees with the `u4` formula and MP clause in the official
TPTP encoding associated with Fitelson and Peltier's positive-implication
problem.

## Strengthened derivability statement

The strengthened file makes uniform substitution explicit:

1. `applySubst` is the homomorphic extension of a variable assignment.
2. `Pair.subst_closed` and `P_subst_closed` prove that the countermodel is
   closed under every such substitution.
3. `U4Derivable` contains arbitrary `u4` instances, modus ponens, and uniform
   substitution.
4. `U4Derivable.sound` maps every derivation into `P`.
5. `u4_not_derives_reflexivity` therefore rules out every `q -> q`, even in
   this explicitly strengthened derivability relation.

This also covers condensed-detachment derivations: each such step is a
uniform substitution followed by modus ponens.

## Internal comparison with positive implicational logic

The same file defines the standard positive-implication `K` and `S` schemes.
`positive_derives_reflexivity` gives a checked `S K K` derivation of every
`q -> q`. The final theorem

```text
U4Formal.u4_not_axiomatizes_positive_implication
```

exhibits a formula derivable in that positive-implication basis but not from
`u4`. This makes the separation statement self-contained within the checked
Lean environment.

## Verification gates

The workflow `u4-strengthened-independent-check-20260820.yml` accepts a commit
only if every gate succeeds:

- Lean 4.33.0 build with warnings treated as errors;
- direct compilation of `U4Formal.lean`;
- Lean's `leanchecker` environment check;
- the independent Rust checker `nanoda`, with `sorryAx` forbidden;
- compiled-environment axiom audit, allowing only the standard
  `propext`, `Classical.choice`, and `Quot.sound` axioms;
- the separately implemented Python symbolic checker; and
- a source scan rejecting `sorry`, `admit`, custom `axiom`, and `unsafe`.

The workflow records the exact commit, tool outcomes, compiler output, source
files, and SHA-256 hashes in a CI artifact. A release should be created only
from a commit for which all gates are green.

## Claim boundary

The formal result supports the following claim:

> The displayed `u4` formula is not a single axiom for positive implicational
> logic under uniform substitution and modus ponens.

It does not by itself establish journal acceptance, consensus in the research
community, or absolute priority over every unpublished result. Public
timestamps, literature comparison, and peer review are separate questions.

## References

- Fitelson and Peltier, *Automated Reasoning for Finding Short Single Axioms for
  Positive Implication*, Journal of Automated Reasoning (2026):
  <https://doi.org/10.1007/s10817-026-09752-1>
- Lean proof validation guidance:
  <https://lean-lang.org/doc/reference/latest/ValidatingProofs/>
- `lean-action` independent-check documentation:
  <https://github.com/leanprover/lean-action>
