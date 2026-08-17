# Master 16-track exact-answer/open-problem tournament

Frozen: 2026-08-17 KST before reading the fresh biomedical smoke-probe outputs from the new crossdomain-branch run.

This master tournament combines **all eight previously defined cross-domain tracks** with **all eight biomedical tracks** that were previously smoke-ranked. `煙霞聖效方` is retained as an explicit T2 subtrack and cannot replace or remove another track.

## Cross-domain tracks

- C1: Nineveh/eBL medical tablet fragment join
- C2: Uibangyuchwi lost medical texts; Jianqifang positive control + Yanxia Shengxiaofang subtrack
- C3: DECRYPT historical medical/alchemical cipher
- C4: Positive implicational logic u4
- C5: Pompeii RePAIR fresco fragment join
- C6: Vesuvius unopened-scroll text recovery
- C7: CASP17 blind structure prediction
- C8: B2[2] subset of Z_100, 13-vs-14 exact finite problem

## Biomedical tracks

- B1: ProteinGym human missense interpretation
- B2: CAMI III metagenome reconstruction/profiling
- B3: Virtual Cell perturbation prediction
- B4: CAGI splicing
- B5: CAGI lentiMPRA/regulatory variant prediction
- B6: CAFA protein-function prediction
- B7: CAMDA AMR genome-to-phenotype
- B8: rare-disease causal-variant/genome interpretation

## Frozen master ranking rules

The master rank is **not** a publication-impact beauty contest. It prioritizes proximity to a decisive, independently checkable answer.

### Evidence stage
- 0: inaccessible, failed execution, or gated with no analyzable dataset
- 1: public/authorized data and objective target/ground truth structure are demonstrably accessible
- 2: a nontrivial known-truth control, retrospective benchmark, or reproducible reconstruction has actually been completed
- 3: an outcome-blind concrete candidate/prediction/witness has been frozen and is awaiting independent truth or expert verification
- 4: a new result has an independently checkable certificate/witness or official blind outcome

### Failure/closure rule
A track whose core proposed method already failed its predeclared harder validation is marked `CLOSED_CURRENT_METHOD`. It may retain its evidence stage but receives a master-score penalty and does not advance unless a genuinely new protocol is created before outcomes.

ProteinGym B1 is **pre-marked CLOSED_CURRENT_METHOD** because the five-model ensemble's prespecified ClinVar historical temporal gate failed and strict novel-human DMS validation was negative. The surviving retrospective family-held-out result is preserved; it is not erased.

### Exactness bonus, frozen before fresh probe results
- +4: formal proof/countermodel, SAT/UNSAT certificate, finite witness with deterministic verifier
- +3: physical fragment join, exact historical cipher mapping, source/page text witness with independent expert/material verification
- +2: blinded experimental/competition truth with numeric scoring
- +1: future annotation/reclassification truth but weaker temporal/circularity control

### Medical/Korean-medicine bonus
- +2: direct Korean/East-Asian medicine or ancient medicine
- +1: biomedical/medical/drug/health relevance
- +0: otherwise

### Master score
`20*evidence_stage + exactness_bonus + medical_bonus - 20*CLOSED_CURRENT_METHOD`

Evidence stage dominates. The score is only a triage rank. A stage-3 candidate is never called solved. `SOLVED`/`DISCOVERY` requires stage 4 and the track-specific certificate from the original protocol.

## Advancement

After all 16 are evaluated, advance every non-gated track tied for the highest evidence stage, plus any lower-stage track with a currently active blind deadline where waiting would irreversibly lose the opportunity (e.g. CASP17). Do not discard the remaining tracks; record their hard stop or next executable gate.
