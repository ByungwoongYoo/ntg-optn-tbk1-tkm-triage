# Separate-Implementation Verification

This clean branch publishes small, inspectable source, verification, metadata,
and audit files for the B₂[2] problem in Z/100Z and a **candidate computational
resolution** with maximum size `M = 13`. Binary ZIPs are intentionally not
committed to Git history.

The new Zenodo record [10.5281/zenodo.22026665](https://zenodo.org/records/22026665)
contains exactly 13 physical files: three PUBLIC_FIXED compact ZIPs, nine
FINAL_FIXED byte-preserving transport parts, and
`P2624_FINAL_FIXED_TRANSPORT_README.txt`. The logical 93,453,983-byte
reconstruction ZIP exists only after reassembly; it is not a fourteenth
physical record file.

The 219,202,341-byte canonical raw-evidence ZIP is not in the new record. It
remains in the separate earlier record
[10.5281/zenodo.21988177](https://zenodo.org/records/21988177) under the remote
filename `B2_Z100_RAW_EVIDENCE_RELEASE_20260818.zip`.

The branch name `b2-z100-certified-resolution-20260820` is an internal release
label. It does not mean that TheoremDB, a journal, or a proof assistant has
certified the result.

## Public status and claim boundary

As checked on 2026-08-20, TheoremDB P2624 remains **Open** and its official
established bound remains

`13 <= M <= 14`.

This project's records are R5659–R5664. R5747–R5749 are unrelated. R6088–R6090
record Jordan Boisclair's independent audit of a distinct, currently
unavailable `P2624_certified_resolution.zip`, SHA-256
`facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35`.
That package has no public source URL or locator and has not been mapped to this
corpus. Accordingly, do not describe it as an independent replay of these raw
files or say that P2624 is officially established or closed.

Separately,
`P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip` is a
separate implementation that regenerates all seven static lift-tree
certificates from the public specification. Its seven tree SHA-256 values and
node counts match the per-file values reported by R6090 exactly. This is a new
reconstruction, not a recovered copy or authentication of Jordan's unavailable
original ZIP.

The fresh exact-replay bundle remains a **same-project replay**. Neither its
fresh run nor the static reconstruction should be described as an independently
administered replay of this project's canonical raw corpus.

The permitted short description is:

> candidate computational resolution with distinct third-party support;
> official P2624 status remains Open pending packet review

See [`docs/PUBLIC_CLAIM_BOUNDARY.md`](docs/PUBLIC_CLAIM_BOUNDARY.md) and
[`docs/THEOREMDB_LIVE_STATUS_20260820_EN.md`](docs/THEOREMDB_LIVE_STATUS_20260820_EN.md).

## Evidence objects and record scope

| Evidence object | Distribution | Bytes | SHA-256 |
|---|---|---:|---|
| `B2_Z100_RAW_EVIDENCE_RELEASE_20260818.zip` | Physical file in old raw record `10.5281/zenodo.21988177`; not in the new record | 219,202,341 | `254031fc2fab17027e389900ca63e704d170eba6ad6861e235db3fe9be46727a` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip` | Logical artifact reconstructed from nine parts; not itself a physical record file | 93,453,983 | `83d80155ef4b2dba5a4def970d13317413b763b3447c00a931f4c624b9e45cb5` |
| `B2_Z100_FINAL_PROOF_AUDIT_PUBLIC_FIXED_20260820.zip` | Physical file in new record `10.5281/zenodo.22026665` | 115,896 | `ae781ccd6d191e1e1a4e66fbe383a864aa7b65a049bf93cf1a335db57136620f` |
| `B2_Z100_PUBLIC_AUDIT_AUXILIARIES_PUBLIC_FIXED_20260820.zip` | Physical file in new record `10.5281/zenodo.22026665` | 88,271 | `47f37758c1f8bd3e7b40efabf2d788113b96b2968ad469e87c3818b02502f60c` |
| `B2_Z100_FRESH_REPLAY_AUDIT_BUNDLE_PUBLIC_FIXED_20260820.zip` | Physical file in new record `10.5281/zenodo.22026665` | 1,087,657 | `052357754e9a04ddf6983f24b59f6575114f371e850cedb3c918e6d842bcf9fb` |

The local filename `B2_Z100_FULL_RAW_EVIDENCE_20260818_REISSUED.zip` is an
optional alias for the verified raw bytes above. The authoritative filename in
the old Zenodo record is `B2_Z100_RAW_EVIDENCE_RELEASE_20260818.zip`; both names
refer to the same 219,202,341-byte SHA-256 identity.

The authoritative machine-readable inventory is
[`release/ARTIFACTS.json`](release/ARTIFACTS.json). Always verify a downloaded
archive against [`release/ARTIFACT_SHA256SUMS.txt`](release/ARTIFACT_SHA256SUMS.txt).

The FINAL_FIXED 93 MB reconstruction is represented in the new record by nine
byte-preserving parts. The exact transport inventory and evidence boundary are
also preserved in
[`release/P2624_FINAL_FIXED_TRANSPORT_README.txt`](release/P2624_FINAL_FIXED_TRANSPORT_README.txt).

| Transport file | Bytes | SHA-256 |
|---|---:|---|
| `P2624_FINAL_FIXED_TRANSPORT_README.txt` | 4,291 | `78e449aea9e605d1afbbfa973cde19ef7f2ae370b19cd176b76792114afe8906` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.00.part` | 10,485,760 | `a7e774e1d9a496aa940bc1e5d85c10c9d89dea6f5b07c5245637392959cade57` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.01.part` | 10,485,760 | `0112e827f76389b9ab5b9ef2a567e91dff634fb2b5096683907a406d2d7aa500` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.02.part` | 10,485,760 | `6f11df8fe8e2599b41077acd07fbdb79d6065fffd8b42985a13e466175edbed4` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.03.part` | 10,485,760 | `65b96e6c0b3fd09087c9a915a9229d94eeeeef50f6c16f3e2e71e909bd851bfe` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.04.part` | 10,485,760 | `ead1cbd52ccb65f6e8ff73660a245a248efee8613094478ec99786a44b1e64a3` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.05.part` | 10,485,760 | `adf6033d99928fd831fb2b626798e558ee6e07e61032b33947c04b5c4f7e4959` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.06.part` | 10,485,760 | `d4870a40483010613d5ef91161bfbb3c8c52e5e50c83e65dda08b27eed4f41c2` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.07.part` | 10,485,760 | `be00304e1ed0d9beb9ed46b8e101b7b50f7a8049fcbcd6563478f78644ab1bcb` |
| `P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.08.part` | 9,567,903 | `b16114d5d2e974b3f2eeb49ab1a19006236e560a977219bc6f5801ba04a3ad88` |

Reassemble in numeric suffix order and verify the authoritative whole-ZIP
identity:

```bash
cat P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.??.part \
  > P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip
printf '%s  %s\n' \
  83d80155ef4b2dba5a4def970d13317413b763b3447c00a931f4c624b9e45cb5 \
  P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip | sha256sum -c -
unzip -t P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip
```

The `.00.part` through `.08.part` objects are transport pieces only; their
sizes sum to 93,453,983 bytes. The reassembled ZIP's size and SHA-256 above are
the scientific artifact identity. The final deterministic-build, ZIP,
manifest, component, witness, quotient, and negative-control checks are in
[`qa/P2624_FINAL_FIXED_QA_REPORT_20260820.json`](qa/P2624_FINAL_FIXED_QA_REPORT_20260820.json).

## Quick checks

The following checks need only this branch:

```bash
python3 verification/verify_corrected_result.py
python3 verification/verify_witness_and_reductions.py
```

They validate the explicit 13-set, elementary reductions, and consistency of
the recorded final run map. They do not rerun the exhaustive searches.

For a raw-corpus audit, download
`B2_Z100_RAW_EVIDENCE_RELEASE_20260818.zip` from the separate old record
`10.5281/zenodo.21988177`, extract it, then run:

```bash
python3 verification/verify_raw_evidence.py /path/to/extracted/raw-evidence
bash verification/clean_room_quick_verify.sh /path/to/extracted/raw-evidence
```

The full exact-search replay is expensive:

```bash
bash verification/reproduce_all_exact_searches.sh /path/to/output
```

Stored DRAT objects can be checked after extracting the raw archive:

```bash
bash verification/replay_drat_all.sh /path/to/extracted/raw-evidence
```

## Repository contents

- `claim/`: corrected claim synopsis and evidence coverage;
- `docs/`: public-status, proof-audit, privacy, and claim-boundary documents;
- `metadata/`: sanitized provenance and DRAT/raw-release identities;
- `source/`: exact-search and SAT source code;
- `verification/`: executable checks and recorded verification summaries;
- `controls/`, `summaries/`, and `replay/`: compact audit outputs;
- `workflows/`: archived workflow definitions, deliberately outside
  `.github/workflows/` so a publication push cannot execute them; and
- `qa/`: final three-package QA plus the FINAL_FIXED reconstruction QA report.

Historical superseded claim files, draft submission text, API responses with
personal email, and unrelated research files are excluded.

## Licensing

- Author-generated proof code, scripts, and workflows: MIT.
- Author-generated documents, data, logs, tables, and manifests: CC BY 4.0.
- Third-party materials remain under their original licenses and are not
  relicensed here.

See `LICENSE_SCOPE.md` and `THIRD_PARTY_NOTICES.md` for file-scope precedence.

## Citation

Please cite [10.5281/zenodo.22026665](https://zenodo.org/records/22026665).
Machine-readable citation metadata is provided in `CITATION.cff`.
