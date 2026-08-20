# 연구 완료 — 사전 정의된 합리적 계산 경로를 소진했으나 엄격한 성공 기준 전체를 독립 holdout에서 통과하지 못했습니다.

- 최종 기계 판정: **FINAL_NEGATIVE_RESULT**
- 엄격한 holdout 생존자: **0개**
- 종료 시 활성 핵심 workflow: **0개**
- 성공 숫자 기준·분할·threshold의 사후 완화: **없음**

## 가장 강한 확인 후보

- `source`: `final_work/evidence/9408181719_LANTERN_CAMI3_V7_UNTOUCHED_HOLDOUT_18_19_FINAL_V2_20260820/decision/HOLDOUT_DECISION.json`
- `original_status`: `PARTIAL`
- `is_holdout`: `True`
- `genome_fraction_gain_pp`: `0.1971779907390392`
- `mean_recovery_gain_pp`: `0.8235867074069952`
- `low_abundance_gain_pp`: `0.343797130120449`
- `relative_chimera_change`: `-0.44187475457456615`
- `longitudinal_drop_pp`: `0.8175554801072593`
- `paired_bootstrap_ci_low_pp`: `0.3989347065603817`
- `paired_bootstrap_ci_high_pp`: `1.3009854950872048`
- `pseudo_all_positive`: `True`
- `metrics_complete`: `True`
- `strict_holdout_success`: `False`

## 주장 경계

SUCCESS requires every original numerical gate on an untouched holdout, no truth leakage, and actual restricted CAMI application plus submission readiness. Toy-only success is labeled separately.

기술 실패는 failure audit에 남겼으며, 성공으로 바꾸기 위해 수치 문턱을 낮추지 않았습니다.
