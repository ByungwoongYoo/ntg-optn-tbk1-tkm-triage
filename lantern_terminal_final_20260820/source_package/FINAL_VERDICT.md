# 연구 종료 판정 보류 — 활성 계산 또는 기술적·경로완료 관문이 남아 있어 성공이나 최종 음성으로 판정할 수 없습니다.

- 기계 판정: **INCOMPLETE_ACTIVE_RUNS**
- 엄격한 holdout 성공 후보: **0개**
- 분석된 decision/result JSON: **254개**
- 활성 workflow: **25개**
- 최신 기술 실패 workflow: **35개**
- 모든 사전 정의 경로 완료: **예**
- CAMI submission-ready 확인: **아니오**
- restricted blind data 적용 확인: **아니오**

## 변경하지 않은 성공 기준

- `minimum_genome_fraction_gain_pp`: `0.5`
- `minimum_mean_genome_recovery_gain_pp`: `1.0`
- `minimum_low_abundance_gain_pp`: `2.0`
- `maximum_relative_chimera_increase`: `0.1`
- `minimum_longitudinal_ablation_drop_pp`: `0.5`
- `minimum_paired_bootstrap_ci_low_pp`: `0.0`

## 가장 강한 확인 후보

- `source`: `final_work/evidence/9408181719_LANTERN_CAMI3_V7_UNTOUCHED_HOLDOUT_18_19_FINAL_V2_20260820/decision/HOLDOUT_DECISION.json`
- `source_sha256`: `2c56f5461bb4d4600302b213f72e014a38b179120cb1c837cce1c7a8bcf7f42a`
- `original_status`: `PARTIAL`
- `is_holdout`: `True`
- `is_development`: `False`
- `leakage_ok`: `True`
- `genome_fraction_gain_pp`: `0.1971779907390392`
- `mean_recovery_gain_pp`: `0.8235867074069952`
- `low_abundance_gain_pp`: `0.343797130120449`
- `relative_chimera_change`: `-0.44187475457456615`
- `longitudinal_drop_pp`: `0.8175554801072593`
- `paired_bootstrap_ci_low_pp`: `0.3989347065603817`
- `paired_bootstrap_ci_high_pp`: `1.3009854950872048`
- `pseudo_all_positive`: `True`
- `metrics_complete`: `True`
- `gate_genome_fraction`: `False`
- `gate_mean_recovery`: `False`
- `gate_low_abundance`: `False`
- `gate_chimera`: `True`
- `gate_longitudinal`: `True`
- `gate_bootstrap`: `True`
- `gate_pseudo`: `True`
- `strict_holdout_success`: `False`

## 경로 완료 감사

- `read_only_pairing`: **완료**
- `additive_holdout`: **완료**
- `hybrid_holdout`: **완료**
- `precision_holdout`: **완료**
- `extension_holdout`: **완료**
- `high_depth_long`: **완료**
- `domain_specific`: **완료**

## 주장 경계

SUCCESS requires every original numerical gate on an untouched holdout, no truth leakage, and actual restricted CAMI application plus submission readiness. Toy-only success is labeled separately.
