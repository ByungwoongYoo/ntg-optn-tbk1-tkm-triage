# TheoremDB P2624 live status — 2026-08-20

The official P2624 status remains **Open**, and the official established bound
remains `13 <= M <= 14`. The live Work view contains 17 items: 5 packet items
and 12 recent-work items.

The contributor records associated with this project's B2 work are
`R5659`–`R5664`. Records `R5747`–`R5749` are unrelated to P2624 and must not be
treated as this project's provenance or submission targets.

## Distinct third-party audit lineage

TheoremDB records `R6088`, `R6089`, and `R6090` attribute to Jordan Boisclair
an independent audit of a **different** artifact:

- file name: `P2624_certified_resolution.zip`;
- SHA-256: `facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35`;
- 7 quotient occupancy types;
- 204,360 admissible raw occupancy vectors;
- 1,341 affine orbits;
- 299,903,736 proof-tree nodes; and
- a negative control in which a corrupted proof byte was rejected.

Those records support `M = 13`, but all three retain the badge
`Needs packet proposal`. In R6090's structured packet, `source.url` and
`source.locator` are `null`; the exact ZIP was not available from a public
locator during this audit.

## Required evidence boundary

The Jordan audit is genuine third-party support recorded by TheoremDB, so it is
incorrect to say that *no* third-party validation exists. It is equally
incorrect to claim that Jordan independently replayed this package's canonical
raw corpus. This package's raw corpus was not independently replayed in this
workspace, and the unavailable `facb…a35` artifact cannot currently be mapped
byte-for-byte to it. The two evidence lineages must remain distinct until the
certified artifact is publicly locatable and compared directly.

Accordingly, do not claim that TheoremDB has closed P2624 or that this raw
corpus is externally reproduced. The accurate statement is: a distinct
third-party certified-package audit supports `M = 13`, while the official
problem status remains Open and the current official bound remains
`13 <= M <= 14`.
