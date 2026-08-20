# Current public status — live audit on 2026-08-20

## Official result

For P2624, TheoremDB still displays **Open** and the established interval

`13 <= M <= 14`.

The official status therefore does not yet decide whether a valid 14-set
exists. The live Work view contains 17 items: 5 packet items and 12 recent-work
items.

Problem page: https://www.theoremdb.org/statements/b2-two-set-z100/

## Project records and record-ID correction

The B2 records associated with this contributor are:

| Records | Role | Displayed state |
|---|---|---|
| R5659–R5661 | first attempt, claim, and artifact | Partial / Supported / Available |
| R5662–R5664 | executable v2 attempt, claim, and artifact | Partial / Supported / Available |

The v2 series points to Zenodo DOI `10.5281/zenodo.21988177` and the corrected
raw-evidence identity. Records `R5747`, `R5748`, and `R5749` are unrelated to
P2624 and must not be used as B2 provenance or submission targets.

## Distinct third-party audit

TheoremDB records R6088–R6090 attribute to Jordan Boisclair an independent
audit of a different certified artifact:

- `P2624_certified_resolution.zip`;
- SHA-256 `facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35`;
- 7 quotient occupancy types, 204,360 raw vectors, 1,341 affine orbits, and
  299,903,736 proof-tree nodes;
- recompilation of the checker and rejection of a deliberately corrupted
  proof byte.

The three records support `M = 13`, but each displays `Needs packet proposal`.
R6090's structured packet has `source.url = null` and
`source.locator = null`, and the exact artifact was not publicly locatable in
this audit.

## Evidence boundary

It is inaccurate to say that no third-party validation exists: R6088–R6090 are
documented independent support. It is also inaccurate to say that this
project's canonical raw corpus was independently replayed by Jordan. The
Jordan audit concerns a distinct, currently unavailable package, and no
byte-for-byte mapping to this raw corpus is established.

Use this wording:

> We provide a candidate computational resolution supported internally by a
> complete recorded search corpus. A distinct third-party certified-package
> audit in TheoremDB also supports `M = 13`, but its exact artifact lacks a
> public locator and has not been mapped to this corpus. P2624 remains Open and
> the official established bound remains `13 <= M <= 14`.

Do not state that P2624 is closed, that this raw corpus is externally
reproduced, or that the result is peer reviewed or formally verified.
