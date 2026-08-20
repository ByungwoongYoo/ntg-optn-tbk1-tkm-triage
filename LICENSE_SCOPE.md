# File-scope license map

The license follows the file's role, not merely the outer ZIP.

## MIT (`LICENSE-CODE-MIT`)

The following author-generated software is licensed under MIT:

- `source/*.py` and `source/*.cpp`;
- `verification/*.py` and `verification/*.sh`;
- `verification/Dockerfile`;
- `tools/*.py`; and
- `workflows/*.yml` and `workflows/*.yaml`.

Historical project code remains MIT when it falls under those paths; its
historical status is independent of its license.

## CC BY 4.0 (`LICENSE-DOCS-DATA-CC-BY-4.0`)

The following author-generated material is licensed under CC BY 4.0:

- author-generated Markdown, TeX, JSON, CSV, TSV, TXT, log, and SHA sidecar files
  that are not one of the explicitly MIT-mapped software files above;
- `audit/`, `claim/`, `controls/`, `historical/`, `metadata/`, `provenance/`,
  `release_materials/`, and `summaries/`;
- author-generated CNF/DRAT objects and computational outputs, if present; and
- `manifest/MANIFEST.json` and `manifest/MANIFEST.sha256`.

## Exclusions and precedence

- These grants apply only to material for which Byungwoong Yoo owns or may
  license the relevant rights.
- CaDiCaL, `drat-trim`, compiler/toolchain components, libraries, cited works,
  and any other third-party source or binary remain under their original
  licenses. They are not relicensed by this package.
- Original third-party copyright and license notices take precedence for the
  corresponding material.
- The two license texts and this scope map are notices of their own terms and
  are not reclassified by the path rules above.
