P2624 B2[2]/Z100 certified reconstruction -- FINAL_FIXED transport inventory
Date: 2026-08-20
Reserved DOI: 10.5281/zenodo.22026665
Record URL: https://zenodo.org/records/22026665

PURPOSE
=======
The logical reconstruction ZIP is 93,453,983 bytes.  It is represented here
as nine byte-preserving 10 MiB transport parts so each file can be uploaded
reliably.  The parts are not separate archives.  Concatenate parts 00 through
08, in numeric order, to reconstruct the logical ZIP exactly.

LOGICAL ZIP (result after concatenation)
========================================
Filename: P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip
Bytes:    93453983
SHA-256:  83d80155ef4b2dba5a4def970d13317413b763b3447c00a931f4c624b9e45cb5

NUMBERED TRANSPORT PARTS
========================
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.00.part
  Bytes: 10485760
  SHA-256: a7e774e1d9a496aa940bc1e5d85c10c9d89dea6f5b07c5245637392959cade57
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.01.part
  Bytes: 10485760
  SHA-256: 0112e827f76389b9ab5b9ef2a567e91dff634fb2b5096683907a406d2d7aa500
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.02.part
  Bytes: 10485760
  SHA-256: 6f11df8fe8e2599b41077acd07fbdb79d6065fffd8b42985a13e466175edbed4
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.03.part
  Bytes: 10485760
  SHA-256: 65b96e6c0b3fd09087c9a915a9229d94eeeeef50f6c16f3e2e71e909bd851bfe
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.04.part
  Bytes: 10485760
  SHA-256: ead1cbd52ccb65f6e8ff73660a245a248efee8613094478ec99786a44b1e64a3
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.05.part
  Bytes: 10485760
  SHA-256: adf6033d99928fd831fb2b626798e558ee6e07e61032b33947c04b5c4f7e4959
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.06.part
  Bytes: 10485760
  SHA-256: d4870a40483010613d5ef91161bfbb3c8c52e5e50c83e65dda08b27eed4f41c2
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.07.part
  Bytes: 10485760
  SHA-256: be00304e1ed0d9beb9ed46b8e101b7b50f7a8049fcbcd6563478f78644ab1bcb
P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.08.part
  Bytes: 9567903
  SHA-256: b16114d5d2e974b3f2eeb49ab1a19006236e560a977219bc6f5801ba04a3ad88

COMPACT AUDIT PACKAGES IN THE SAME DEPOSIT
==========================================
B2_Z100_FINAL_PROOF_AUDIT_PUBLIC_FIXED_20260820.zip
  Bytes: 115896
  SHA-256: ae781ccd6d191e1e1a4e66fbe383a864aa7b65a049bf93cf1a335db57136620f
B2_Z100_PUBLIC_AUDIT_AUXILIARIES_PUBLIC_FIXED_20260820.zip
  Bytes: 88271
  SHA-256: 47f37758c1f8bd3e7b40efabf2d788113b96b2968ad469e87c3818b02502f60c
B2_Z100_FRESH_REPLAY_AUDIT_BUNDLE_PUBLIC_FIXED_20260820.zip
  Bytes: 1087657
  SHA-256: 052357754e9a04ddf6983f24b59f6575114f371e850cedb3c918e6d842bcf9fb

REASSEMBLY
==========
POSIX shell (Linux/macOS):

  cat P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip.??.part > P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip
  sha256sum P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip

Cross-platform Python 3 (streams in numeric order):

  python3 -c "from pathlib import Path; n='P2624_independent_certified_reconstruction_FINAL_FIXED_20260820.zip'; o=open(n,'wb'); [o.write(Path(f'{n}.{i:02d}.part').read_bytes()) for i in range(9)]; o.close()"

Then verify that the reconstructed ZIP has the exact byte count and SHA-256
shown under LOGICAL ZIP.  Test the archive with `unzip -t` (or an equivalent
ZIP integrity checker), extract it, and run:

  python3 P2624_independent_certified_reconstruction/src/verify_manifest.py P2624_independent_certified_reconstruction/SHA256SUMS

EVIDENCE BOUNDARY
=================
This is a newly assembled independent reconstruction package with its own
archive hash.  It is not, and does not claim byte identity with, the unavailable
original P2624_certified_resolution.zip container reported under SHA-256
facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35.
The seven certificate tree components do match the seven SHA-256 values
published in TheoremDB record R6090.  The fresh-replay bundle is same-project
replay evidence, not an external third-party reproduction or certification.
