# Privacy exclusion record

## Excluded source container

The original file
`github_artifact_9330620059_internal_audit_v6.zip` is not distributed in this
bundle.

- bytes: `99081`
- SHA-256:
  `79a787967a6c421e86f51248816705de7fc4abdf15598cd050c700cea9d441a0`

## Excluded member

- path: `REPLAY_RUN_API.json`
- bytes: `14193`
- SHA-256:
  `0a6a0f4f48870d7b3d9f9f6dec07f777e35656f16bc7ed71516ce7c1845a0753`
- reason: contains a personal email address in copied API metadata

No redacted approximation of that JSON is supplied. All other extracted
members are copied unchanged into `sanitized_internal_audit/`. Original
internal file/hash lists are preserved as provenance records and may name the
excluded member; they are not the distribution manifest. The authoritative
distributed inventory is `manifest/MANIFEST.json`.

The original full-replay artifact included in `artifacts/` was separately
scanned and contained no email-address match.
