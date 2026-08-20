# Corrections and superseded records

This file distinguishes final evidence from historical material preserved in
the canonical raw archive. Historical records remain valuable provenance, but
they must not be silently promoted to final proof dependencies.

## 1. Final tail run correction

Historical files:

- `repository_snapshot/crossdomain/B2_Z100_FINAL_RESOLUTION_20260817.md`
- `repository_snapshot/crossdomain/B2_Z100_FINAL_RESULTS_20260817.json`

refer to run `32004571194` for a `t=10,...,45` tail. That is a stale v7-era
record.

The corrected final dependency is:

- run `32038657803`;
- solver generation v8;
- final range `t=11,...,45`;
- 35 final branches.

The final `t=10` branch is in the v8 head run `32038183046`, and its recorded
aggregate node count is `537251446`.

## 2. Historical v8 hard branches

The v8 head archive retains `t=2,3,4` for provenance. They are not final upper-
bound evidence:

- historical v8 `t=2`: completed, superseded;
- historical v8 `t=3`: incomplete timeout, superseded and never proof of UNSAT;
- historical v8 `t=4`: completed, superseded.

All final evidence for `t=2,3,4` is supplied by 255 completed v10 `(t,u)`
subbranches in run `32040699080`. The separate v10 `(7,8)` artifact is a
regression control, not part of those 255 branches.

## 3. Status-language correction

The phrase “reproducible computational mathematical resolution” in a
historical narrative can be misread as a claim of completed independent
reproduction. The unambiguous current wording is:

> candidate computational resolution with distinct third-party support;
> official P2624 status Open pending packet review

The project-managed fresh DRAT run is a valuable clean-environment replay, but
it is not a review by an independent third party. TheoremDB R6088–R6090 do
record Jordan Boisclair's independent support from a distinct
`P2624_certified_resolution.zip` (SHA-256 `facbdfec…8a35`), but that package is
not publicly locatable and has not been mapped to this raw corpus. Do not mark
this corpus `reproduced`, `formally_verified`, or `Lean verified` on either
basis.

## 4. Manifest corrections

Two historical manifest issues are preserved and disclosed:

1. 370 artifact-local `SHA256SUMS.txt` files include a self-referential line
   created while each hash file was being written. Those self-lines are not
   authoritative.
2. The historical pre-upload manifest listed 4,132 files, while 25 hidden
   `.github/workflows` paths were omitted by `actions/upload-artifact`.

The authoritative integrity record is:

- `manifest/RELEASE_MANIFEST.json`;
- manifest SHA-256
  `dade214232ffdc2ce3386dcb68ac0e8325004b5345c256e853301e3e76ccdca5`;
- 4,155 governed files;
- 518,219,226 governed uncompressed bytes.

## 5. Outer ZIP recovery correction

The surviving 333,280,842-byte object is not the intended release byte stream.
It contains a complete 219,202,341-byte canonical ZIP followed by a duplicate
114,078,501-byte tail. The canonical ZIP is recovered by taking the exact
prefix through the first valid EOCD.

- residual SHA-256: `cc3104b4273ff71456f5dd4f77c14c82b317ae3e2c4b13e72166a998a4f74094`
- canonical SHA-256: `254031fc2fab17027e389900ca63e704d170eba6ad6861e235db3fe9be46727a`

Do not rezip or edit the canonical archive when producing the recovered full
raw-evidence file.

## 6. Public-release blockers

- Four raw `RUN.json` files contain a personal email address.
- Draft Zenodo/TheoremDB/preprint materials contain unfilled DOI, GitHub, and
  release-SHA placeholders.
- MIT for author-generated proof code/scripts/workflows and CC BY 4.0 for
  author-generated docs/data/logs/tables/manifests are declared in the package
  scope map; third-party materials retain their original licenses.
- A primary-literature and priority audit is not complete.
- Corpus-specific independent mathematical review and a separately
  administered replay or direct comparison with the unavailable Jordan
  package remain pending.

Redaction of the canonical raw release would change file identities and
invalidate its manifest. If redaction is approved, create a separately named
derivative with a new manifest and retain a provenance link to the canonical
archive.
