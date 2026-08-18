# B2[2]/Z100 global normalization and branch-coverage proof

This note separates the mathematical reduction from computational leaf outcomes. It proves that every hypothetical valid 14-element set has an image in exactly one of the canonical branch families used by the proof-producing workflows. It does not by itself certify any terminal UNSAT result.

## Problem

For a subset `A` of `Z/100Z`, define

```text
r_A(d) = |{(a,b) in A x A : a != b and a-b = d mod 100}|,
```

for every nonzero residue `d`. A `B2[2]` set satisfies `r_A(d) <= 2` for every nonzero `d`. We seek the maximum cardinality `M`.

## Lemma 1 — every valid 14-set has a unit difference

A unit modulo 100 is an odd residue not divisible by 5. Suppose a valid 14-set had no ordered difference coprime to 100.

### Same parity

All `14*13 = 182` ordered differences would be nonzero even residues. There are 49 such residues and each has capacity 2, so total capacity is only `49*2 = 98`, contradiction.

### Both parities

Every cross-parity difference is odd. If no difference is a unit, every cross-parity difference is divisible by 5 and therefore congruent to 5 modulo 10. Fix an occupied residue class of one parity modulo 10. Every occupied class of the opposite parity must be exactly 5 away. Thus at most one class of each parity is occupied and all differences are multiples of 5. There are 19 nonzero multiples of 5 modulo 100, total capacity `19*2 = 38`, again contradicting 182.

Hence some pair difference is a unit.

## Lemma 2 — affine normalization to `{0,1}`

Let `a-b` be a unit. Translation by `-b`, followed by multiplication by `(a-b)^(-1)`, is an automorphism of `Z/100Z` and maps the pair to `{0,1}`. It bijectively permutes nonzero differences and preserves every multiplicity bound. Therefore every hypothetical valid 14-set has an affine image containing 0 and 1.

Uniqueness of the normalized image is not required; existence suffices for a nonexistence proof.

## Lemma 3 — reflection reduction for the third selected point

Write a normalized set in increasing representatives as

```text
0 < 1 < t < x4 < ... < x14 = L < 100.
```

The involution `x -> 1-x (mod 100)` preserves `{0,1}` and exchanges the gap `t-1` with the wrap gap `100-L`. Choose the orientation satisfying

```text
t-1 <= 100-L.
```

There are 11 positive gaps from `t` to `L`, hence `L >= t+11`. The reflection condition gives `L <= 101-t`. Therefore

```text
t+11 <= 101-t,
```

so

```text
2 <= t <= 45.
```

Thus the 44 branches `t=2,...,45` cover every normalized candidate. A reflection tie creates duplication, not a coverage hole.

## Lemma 4 — complete next-selected ranges

For fixed third point `t`, the fourth point `u` has ten points following it. Hence

```text
t+1 <= u <= 91-t.
```

For fixed `t,u`, the fifth point `v` has nine points following it. Hence

```text
u+1 <= v <= 92-t.
```

More generally, let a canonical prefix have length `m`:

```text
P = [0,1,t,p4,...,pm],   3 <= m < 14.
```

With reflection bound `max(A) <= 101-t`, there remain `14-m` selected points. The next selected representative `q` therefore ranges over the complete interval

```text
pm + 1 <= q <= (101-t) - ((14-m)-1).
```

The next selected representative is unique. Consequently these child branches are disjoint and their union is exactly the parent prefix branch.

## Consequence

After affine normalization the set contains consecutive residues 0 and 1, so the minimum cyclic gap is 1. A checked proof that all 44 canonical `t` branches are UNSAT excludes every 14-element `B2[2]` subset of `Z/100Z`.

Together with an independently checked 13-element witness, this yields `M=13` computationally once the full terminal search evidence is audited without gaps.

## Evidence update — complete frozen-source replay

GitHub Actions run `32100985982` freshly compiled the frozen exact-search sources and completed all 255 `(t,u)` branches for `t=2,3,4` and all 41 `t` branches for `t=5,...,45`. The aggregate records 13,668,473,356 search nodes, zero timeouts, and zero 14-set witnesses. This replay was performed by the same research project and therefore is not independent external reproduction. Its independent file-, coverage-, and result-level audit is performed by `b2_z100_verify_fresh_replay_v6.py`.

The permitted successful internal status is therefore `COMPUTATIONAL PROOF CANDIDATE FOR M=13 — PENDING INDEPENDENT EXTERNAL REVIEW`, not `Lean proved`, `peer-reviewed solved`, or `independently externally verified`.
