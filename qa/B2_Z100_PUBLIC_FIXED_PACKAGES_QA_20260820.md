# B2/Z100 PUBLIC_FIXED package QA — 2026-08-20

The three publication-facing compact archives passed deterministic double
build, ZIP CRC and safe-path checks, fixed-timestamp checks, internal-manifest
verification, syntax checks, and privacy/secret/placeholder scans.

| Public archive | Bytes | SHA-256 | Internal manifest SHA-256 | Manifest files |
|---|---:|---|---|---:|
| `B2_Z100_FINAL_PROOF_AUDIT_PUBLIC_FIXED_20260820.zip` | 115,896 | `ae781ccd6d191e1e1a4e66fbe383a864aa7b65a049bf93cf1a335db57136620f` | `8a138bced85debac924b5cdf8db4520a27de54a8897e768d068cb392be6e8511` | 66 |
| `B2_Z100_PUBLIC_AUDIT_AUXILIARIES_PUBLIC_FIXED_20260820.zip` | 88,271 | `47f37758c1f8bd3e7b40efabf2d788113b96b2968ad469e87c3818b02502f60c` | `4d8b180231f570b120094de00b46d33dd517ac82455335020658b03b0f16fcc3` | 47 |
| `B2_Z100_FRESH_REPLAY_AUDIT_BUNDLE_PUBLIC_FIXED_20260820.zip` | 1,087,657 | `052357754e9a04ddf6983f24b59f6575114f371e850cedb3c918e6d842bcf9fb` | `95d4129fd1851f34b8d139460b957ed586df7f64e8b30b94cf095251c3c14857` | 30 |

All three manifests report zero email, secret, and placeholder matches. The
proof and auxiliary corrected-result checks pass for the 13-point witness,
maximum ordered-difference multiplicity 2, 67 mod-10 supports, 44 normalized
third branches, and the three distinct external-audit records.

The fresh-replay archive explicitly retains the same-project/not-external
boundary. It is not an independently administered replay of the canonical raw
corpus.

The superseded compact ZIPs dated only `20260818` or lacking `PUBLIC_FIXED` are
preserved byte-for-byte but are excluded from the public Zenodo and GitHub
artifact indexes.

## Static-certificate reconstruction transport

`P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip` is
93,453,983 bytes
with SHA-256
`83d80155ef4b2dba5a4def970d13317413b763b3447c00a931f4c624b9e45cb5`.
Its ZIP integrity passes. The verified Zenodo transport uses nine
byte-preserving files, `.00.part` through `.08.part`; their sizes sum to the
whole archive's 93,453,983 bytes. Exact per-part sizes and SHA-256 values are
recorded in `release/P2624_FINAL_FIXED_TRANSPORT_README.txt` and
`release/ARTIFACTS.json`. Scientific identity is determined only after numeric
suffix-order concatenation and verification of the whole-ZIP hash.

That package independently reconstructs all seven certificate files and
matches the per-file R6090 tree hashes. It is not a recovered copy or
authentication of the unavailable original `P2624_certified_resolution.zip`.

The complete machine-readable PUBLIC_FIXED QA is
`B2_Z100_PUBLIC_FIXED_PACKAGES_QA_20260820.json`.
The separate FINAL_FIXED reconstruction QA is
`P2624_FINAL_FIXED_QA_REPORT_20260820.json`.
