# B2[2] in Z/100Z: corrected computational resolution record

**Status:** candidate computational resolution with distinct third-party support; official P2624 status Open pending packet review  
**Evidence grade:** executable

## Result represented by the evidence

The preserved computation records no 14-element `B_2[2]` subset of
`Z/100Z`. The set

`{0,5,7,31,58,61,62,63,72,80,84,91,97}`

has 13 distinct elements, and every nonzero ordered difference occurs at most
twice. Therefore the recorded exact maximum is 13.

This is a computational result with complete raw evidence and replay
instructions. This exact canonical raw corpus has not yet received an
independently administered full replay, specialist peer review, or complete
formal verification.

TheoremDB R6088–R6090 separately record Jordan Boisclair's independent audit
supporting `M = 13` from `P2624_certified_resolution.zip`, SHA-256
`facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35`.
That distinct package has no public source URL or locator, shows
`Needs packet proposal`, and has not been mapped to this corpus. Thus the
Jordan records are genuine external support but not proof that this raw corpus
was replayed.

## Lower bound

For the displayed 13-set, direct enumeration of all ordered pairs `a != b`
gives maximum multiplicity 2 among the 99 nonzero residues modulo 100. Hence a
valid set of size 13 exists.

## Unit-difference reduction

Suppose a hypothetical 14-set had no difference coprime to 100. A machine-
audited classification of its occupied residue classes modulo 10 gives two
possibilities:

1. all occupied classes have one parity; or
2. the mixed-parity support is a pair separated by 5 modulo 10.

In the first case all nonzero differences are even. There are 49 nonzero even
residues modulo 100, with total multiplicity capacity `49*2=98`, less than the
required `14*13=182` ordered differences. In the second case all differences
are multiples of 5. The 19 nonzero multiples of 5 have capacity `19*2=38`,
again less than 182. Thus every hypothetical 14-set contains a unit difference.

Translation and multiplication by its inverse normalize that difference to a
contained pair `{0,1}` without changing difference multiplicities. The
documented reflection reduction leaves the 44 cases `t=2,...,45`.

## Complete final g=1 search

The 44 normalized third-element branches are covered by the following final
dependencies:

- `t=2,3,4`: 255 disjoint v10 `(t,u)` subbranches from run `32040699080`;
- `t=5,...,10`: 6 v8 branches from run `32038183046`;
- `t=11,...,45`: 35 v8 branches from run `32038657803`.

Each final raw result is recorded as completed exhaustively, not timed out, and
without a 14-element witness. No timeout or missing output is interpreted as a
proof of nonexistence.

The old run `32004571194` is not the final tail dependency. Historical v8
files for `t=2,3,4`, including an incomplete v8 `t=3` run, are superseded by
the final v10 evidence.

## Redundant cross-checks

- A separate v9 exact solver completed 67 normalized `g=2,3` branches in run
  `32040180627` without a witness.
- Actual CNF and DRAT proof objects for `g=4,5,6,7` are preserved.
- Their original and later fresh `drat-trim` logs report `s VERIFIED`.
- The fresh replay records `drat-trim` commit
  `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`.

These checks are redundant consistency evidence after the global unit-
difference normalization; they do not replace the complete `g=1` branch
partition.

## Canonical evidence identity

- release ZIP bytes: `219202341`
- release ZIP SHA-256:
  `254031fc2fab17027e389900ca63e704d170eba6ad6861e235db3fe9be46727a`
- release-manifest SHA-256:
  `dade214232ffdc2ce3386dcb68ac0e8325004b5345c256e853301e3e76ccdca5`
- manifest files: `4155`
- manifest bytes: `518219226`

## Claim boundary

The evidence supports an executable computational resolution. Distinct
third-party support exists in R6088–R6090, but the exact audited package is not
publicly locatable and no byte-level relationship to this corpus is
established. The evidence does not establish worldwide novelty, journal
acceptance, full Lean formalization, or independent replay of this exact raw
corpus. P2624 remains Open with official bound `13 <= M <= 14`.
