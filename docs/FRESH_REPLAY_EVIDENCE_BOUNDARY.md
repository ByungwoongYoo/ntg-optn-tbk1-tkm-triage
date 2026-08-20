# Evidence boundary

This bundle proves only the following provenance and verification facts:

1. the included full-replay ZIP is byte-identical to recovered GitHub artifact
   9312169839 at SHA-256 `8baab190…79bf`;
2. its recorded summary covers 255 v10 and 41 v8 descriptors, totaling
   13,668,473,356 nodes with zero timeouts and zero 14-set witnesses;
3. the sanitized internal-audit tree and local recheck agree on the aggregate
   counts and valid 13-set; and
4. the distributed bundle excludes the only file known here to contain a
   personal email address.

This does **not** establish that:

- the replay was administered by an independent third party;
- every mathematical reduction and pruning rule has survived specialist peer
  review;
- this replay is the same computation or artifact as Jordan Boisclair's
  unavailable `facbdfec…8a35` certified package;
- P2624 is officially closed; or
- the result is formally verified or journal accepted.

The accurate description is a same-project fresh full exact replay and
sanitized internal audit supporting a candidate `M = 13`, alongside distinct
third-party support in R6088–R6090. Official P2624 status remains Open with
bound `13 <= M <= 14`.
