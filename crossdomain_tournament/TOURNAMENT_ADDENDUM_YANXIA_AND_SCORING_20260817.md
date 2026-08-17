# Cross-domain tournament addendum — Yanxia Shengxiaofang + frozen scoring

Date frozen: 2026-08-17 KST, before the unified all-track run.

## Scope

The original eight tracks in `TOURNAMENT_PROTOCOL_20260817.md` remain active and must all be run. No track is removed because another track looks promising.

Within T2 (lost East Asian medical texts), add an explicit subtrack:

### T2b — `煙霞聖效方` (Yanxia Shengxiaofang / 연하성효방)

Goal: recover all explicit `煙霞聖效方` source-attributed blocks accessible in the scanned Uibangyuchwi corpus, segment candidate formula names, and compare them against the rest of the scanned corpus. A candidate is not a discovery unless source boundary, formula boundary, page/text witness, parallel-source comparison, and prior-art audit all survive.

The pre-existing `簡奇方` extraction remains the positive control for T2. T2b may advance only if the control extraction works.

## Frozen tournament scoring

Ranking is by evidence stage first. Ties use control status, then medical/biomedical relevance as a small bonus. Candidate count never overrides a failed control.

Evidence stage:
- 0 = inaccessible / execution failure / legal or authentication gate with no analyzable data
- 1 = public data or target inventory successfully accessed
- 2 = nontrivial known-truth/control reconstruction or blind benchmark completed
- 3 = control passed and an unverified novel/discovery candidate was generated under frozen rules
- 4 = independently verified witness/certificate or official blind result

Numerical score used only for ordering ties:
`10 * evidence_stage + 3 * control_pass + 2 * candidate_ready + medical_bonus`

`medical_bonus` is fixed before execution:
- T1 Nineveh medicine = 2
- T2 Uibangyuchwi/Jianqifang/Yanxia = 2
- T3 historical medical/alchemical ciphers = 1
- T4 u4 logic = 0
- T5 RePAIR = 0
- T6 Vesuvius = 0
- T7 CASP17 = 2
- T8 B2[2] = 0

`candidate_ready=1` means a concrete candidate/witness target is saved for falsification, not that it is novel or correct.

## Hard claim boundary

No track is called solved at stage 3. Stage 4 requires the track-specific certificate from the original protocol. A timeout is no result. Text similarity alone is not a fragment join. A source-name hit alone is not reconstruction of a lost book. A CASP candidate not officially timestamped/submitted before truth release is not blind validation. A B2 near-miss is not a witness.
