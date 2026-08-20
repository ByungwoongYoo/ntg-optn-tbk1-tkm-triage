# Final proof audit — B2[2] in Z/100Z

**Audit date:** 2026-08-20  
**Claim status:** candidate computational resolution with distinct third-party support; official P2624 status Open pending packet review  
**Evidence grade:** executable

## Executive conclusion

The canonical raw release is internally intact and supports a complete
computational nonexistence record for 14-element `B_2[2]` subsets of
`Z/100Z`, together with an independently checkable 13-element witness. The
archive-integrity and internal-consistency checks pass subject to the disclosed
historical warnings below.

This audit is not an external peer review. It did not perform a new full exact
replay, did not independently rederive every pruning lemma, and did not
establish novelty or priority. Consequently, the result must remain described
as a **candidate computational resolution with distinct third-party support;
official P2624 status Open pending packet review**.

After this package audit, the live TheoremDB review identified R6088–R6090,
which record Jordan Boisclair's independent support for `M = 13` using a
different unavailable `P2624_certified_resolution.zip` (SHA-256
`facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35`).
Those records do not establish a replay of this canonical corpus: their source
URL and locator are null, all show `Needs packet proposal`, and no direct
artifact comparison is available. P2624 remains Open with official bound
`13 <= M <= 14`.

## 1. Canonical release integrity

The surviving file is 333,280,842 bytes with SHA-256
`cc3104b4273ff71456f5dd4f77c14c82b317ae3e2c4b13e72166a998a4f74094`.
Binary ZIP-structure inspection shows:

- first complete ZIP ends at byte `219202341`;
- its SHA-256 is
  `254031fc2fab17027e389900ca63e704d170eba6ad6861e235db3fe9be46727a`;
- the remaining `114078501` bytes are a duplicate of the original tail,
  starting from original offset `105123840` through the original EOCD;
- the canonical prefix contains 4,157 ZIP entries;
- all entries pass CRC testing;
- all 4,155 rows of `manifest/RELEASE_MANIFEST.json` match both recorded byte
  length and SHA-256;
- verified manifest-governed bytes total `518219226`;
- `RELEASE_MANIFEST.json` itself hashes to
  `dade214232ffdc2ce3386dcb68ac0e8325004b5345c256e853301e3e76ccdca5`.

Therefore the canonical release is recoverable deterministically by truncating
at the first EOCD. Recompression is neither necessary nor desirable.

## 2. Mathematical reduction represented by the archive

The explicit 13-element witness is

`{0,5,7,31,58,61,62,63,72,80,84,91,97}`.

Direct ordered-difference counting gives maximum nonzero multiplicity 2.

For the upper bound, the recorded elementary lemma is that every hypothetical
14-set contains a difference coprime to 100. If no unit difference existed,
the admissible residue support modulo 10 would be either all one parity or a
mixed pair separated by 5. In the first case 182 ordered differences would
have to fit into 49 nonzero even residues with capacity 98; in the second they
would have to fit into 19 nonzero multiples of 5 with capacity 38. Both are
impossible. Translation and multiplication by the inverse of a unit difference
normalize a pair to `{0,1}`. The documented reflection reduction then produces
44 normalized third-element branches, `t=2,...,45`.

The archive contains a machine audit of the modulo-10 support classification
and witness. An independent reviewer should still audit the reflection step,
branch disjointness, and every safe-pruning argument at source level.

## 3. Final exact-search coverage

| Evidence | Final dependency | Count | Role |
|---|---:|---:|---|
| v10 `g=1`, `t=2,3,4` | run `32040699080` | 255 final `(t,u)` subbranches | Primary hard-branch evidence |
| v8 `g=1`, `t=5,...,10` | run `32038183046` | 6 branches | Primary evidence |
| v8 `g=1`, `t=11,...,45` | run `32038657803` | 35 branches | Primary evidence |
| v9 `g=2,3` | run `32040180627` | 67 branches | Redundant exact cross-check |
| CNF/DRAT `g=4,...,7` | recorded source runs | 4 cases | Redundant proof-producing SAT cross-check |

The v10 run also contains one separate `(t,u)=(7,8)` regression-control
artifact, so its source-artifact category count is 256 rather than 255.

For all final exact branches, the raw verifier requires:

- `completed_exhaustively == true`;
- `timed_out == false`;
- `witness_found == false`;
- exact contiguous child-branch ranges;
- matching aggregate node formulas;
- source snapshot hashes matching the repository snapshots.

The stored raw-verification report records all these checks as PASS. Final
incomplete branches: 0. Final timeouts treated as UNSAT: 0. Final 14-element
witnesses: 0.

## 4. Stored CNF/DRAT evidence

The raw release contains the actual CNF and DRAT objects for `g=4,5,6,7`,
CaDiCaL exit code 20 and logs, original `drat-trim` logs, and later fresh replay
logs. Both generations of checker logs report `s VERIFIED` for each case.

The fresh replay records `drat-trim` commit
`2e3b2dc0ecf938addbd779d42877b6ed69d9a985`. This audit verified the stored
file identities and archived verification records; it did not itself rerun the
502 MB `g=4` certificate. Exact paths, sizes, hashes, and clause counts are in
`DRAT_PROOF_INDEX.json`.

## 5. Artifact provenance and counts

The authoritative artifact index contains 373 GitHub source artifacts:

- 256 v10 artifacts;
- 44 v8 artifacts, including historical `t=2,3,4`;
- 67 v9 artifacts;
- 2 collector summaries;
- 4 DRAT artifacts.

Recorded GitHub API digests match the downloaded artifact archive hashes. Four
principal exact-search run snapshots and seven selected workflow snapshots are
preserved. The underlying GitHub raw collection and fresh DRAT replay are
identified in `provenance/RAW_RELEASE_IDENTITY.json`.

## 6. Disclosed corrections and warnings

1. The historical resolution Markdown and result JSON name run `32004571194`
   as the tail source. That v7 record is stale. The final v8 tail run is
   `32038657803`.
2. Historical v8 `t=2,3,4` files are superseded by v10. The v8 `t=3` timeout is
   not final evidence and is not a nonexistence result.
3. 370 artifact-local `SHA256SUMS.txt` files contain a self-referential line
   generated while the file was being written. Those self-lines are ignored;
   the corrected release manifest is authoritative.
4. The historical pre-upload manifest listed 4,132 files, but 25 hidden
   `.github/workflows` paths were omitted by `actions/upload-artifact`.
   Preserved final-run workflow snapshots and the release manifest supersede
   that historical list.
5. Four raw `RUN.json` files contain a personal email address. This does not
   affect the computation but blocks automatic public release pending an
   explicit preservation/redaction decision.
6. The archived quick-replay log contains unrelated spreadsheet-runtime warmup
   diagnostics before a PASS result. A clean rerun should replace it in any
   new public-facing package.

## 7. What was and was not established

Established internally:

- canonical archive identity and full manifest integrity;
- explicit valid 13-set;
- complete recorded final branch coverage;
- no final timeout or incomplete branch used as proof;
- stored original and fresh DRAT-verification records;
- exact source, workflow, artifact, and run provenance.

Not established for this canonical raw corpus by this audit:

- independently administered full replay;
- independent implementation agreement;
- complete human proof of every source-level pruning rule;
- worldwide novelty or priority;
- journal peer review or acceptance;
- full Lean/Coq/Isabelle formalization;
- cross-artifact equivalence with the distinct Jordan certified package.

The package now declares MIT for author-generated proof code, scripts, and
workflows; CC BY 4.0 for author-generated documents, data, logs, tables, and
manifests; and original-license precedence for third-party materials.

## 8. Required next validation

1. Reextract the canonical ZIP into an immutable, unique directory and rerun
   `verification/verify_raw_evidence.py` and
   `verification/clean_room_quick_verify.sh` after confirming the 502,355,512
   byte `g=4` DRAT file is complete.
2. Build `drat-trim` at the pinned commit and replay all four stored proofs.
3. Have an additive-combinatorics reviewer audit the unit-difference lemma,
   reflection reduction, branch partition, and every pruning rule.
4. Prefer a separately administered full replay or independent solver/encoding
   for this corpus. If the exact `facb…8a35` package becomes available, compare
   both evidence lineages directly rather than assuming equivalence.
5. Complete primary-literature, PII, secret, placeholder, and third-party
   redistribution audits before any public action; preserve the applied
   file-scope license map.
