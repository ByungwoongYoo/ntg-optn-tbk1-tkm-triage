#!/usr/bin/env python3
"""Aggregate immutable LANTERN evidence and apply the original strict gates.

This script never tunes an assembly, threshold, sample split, or model. It only reads
completed evidence artifacts, finds development and untouched-holdout decisions, and
applies the pre-specified numerical gates. Missing metrics fail a complete success gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

THRESHOLDS = {
    "minimum_genome_fraction_gain_pp": 0.5,
    "minimum_mean_genome_recovery_gain_pp": 1.0,
    "minimum_low_abundance_gain_pp": 2.0,
    "maximum_relative_chimera_increase": 0.10,
    "minimum_longitudinal_ablation_drop_pp": 0.5,
    "minimum_paired_bootstrap_ci_low_pp": 0.0,
}

ALIASES = {
    "gf_gain": [
        "genome_fraction_gain_percentage_points", "genome_fraction_gain_pp",
        "gf_gain_pp", "genome_fraction_gain",
    ],
    "mean_gain": [
        "mean_genome_recovery_gain_percentage_points", "mean_recovery_gain_pp",
        "mean_genome_recovery_gain_pp", "mean_gain_pp",
    ],
    "low_gain": [
        "low_abundance_gain_percentage_points", "low_abundance_gain_pp",
        "low_abundance_recovery_gain_pp", "low_abundance_mean_gain_pp",
    ],
    "chimera_change": [
        "relative_chimera_change", "relative_chimera_increase",
        "relative_misassembly_change", "chimera_relative_change",
    ],
    "longitudinal_drop": [
        "longitudinal_ablation_drop_percentage_points", "longitudinal_drop_pp",
        "longitudinal_ablation_drop_pp", "ablation_drop_pp",
    ],
    "ci_low": [
        "paired_bootstrap_ci_low", "paired_bootstrap_95_ci_low",
        "bootstrap_ci_low", "paired_ci_low",
    ],
    "ci_high": [
        "paired_bootstrap_ci_high", "paired_bootstrap_95_ci_high",
        "bootstrap_ci_high", "paired_ci_high",
    ],
}

DECISION_NAME = re.compile(r"(decision|verdict|result|summary|gate)", re.I)
HOLDOUT_WORD = re.compile(r"(untouched|holdout|validation|reserve|blind)", re.I)
DEVELOPMENT_WORD = re.compile(r"(development|train|training|model_freeze)", re.I)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--evidence-root", type=Path, required=True)
    p.add_argument("--runs-json", type=Path, required=True)
    p.add_argument("--artifacts-json", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def scalar(v: Any) -> float | None:
    if isinstance(v, bool) or v is None:
        return None
    if isinstance(v, (int, float)) and math.isfinite(float(v)):
        return float(v)
    if isinstance(v, str):
        s = v.strip().replace("%", "")
        try:
            x = float(s)
            return x if math.isfinite(x) else None
        except ValueError:
            return None
    return None


def flatten(obj: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out[key.lower()] = v
            out.update(flatten(v, key))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            key = f"{prefix}.{i}" if prefix else str(i)
            out[key.lower()] = v
            out.update(flatten(v, key))
    return out


def lookup(flat: dict[str, Any], aliases: Iterable[str]) -> float | None:
    aliases_l = [a.lower() for a in aliases]
    exact: list[tuple[int, float]] = []
    for key, value in flat.items():
        tail = key.rsplit(".", 1)[-1]
        if tail in aliases_l:
            x = scalar(value)
            if x is not None:
                exact.append((len(key), x))
    if exact:
        exact.sort()
        return exact[0][1]
    for key, value in flat.items():
        if any(a in key for a in aliases_l):
            x = scalar(value)
            if x is not None:
                return x
    return None


def bool_lookup(flat: dict[str, Any], fragments: Iterable[str]) -> bool | None:
    fragments_l = [x.lower() for x in fragments]
    for key, value in flat.items():
        if any(x in key for x in fragments_l):
            if isinstance(value, bool):
                return value
            if isinstance(value, str) and value.strip().lower() in {"true", "false", "pass", "fail", "yes", "no"}:
                return value.strip().lower() in {"true", "pass", "yes"}
    return None


def ci_from_obj(obj: Any, flat: dict[str, Any]) -> tuple[float | None, float | None]:
    lo = lookup(flat, ALIASES["ci_low"])
    hi = lookup(flat, ALIASES["ci_high"])
    if lo is not None or hi is not None:
        return lo, hi
    for key, value in flat.items():
        if "paired_bootstrap_95_ci" in key or "paired_bootstrap_ci" in key:
            if isinstance(value, list) and len(value) >= 2:
                return scalar(value[0]), scalar(value[1])
    return None, None


def pseudo_positive(obj: Any, flat: dict[str, Any]) -> bool | None:
    explicit = bool_lookup(flat, ["pseudo_novel_all_positive", "pseudo_all_positive", "pseudo_novel_or_isolation_gain"])
    if explicit is not None:
        return explicit
    tiers = None
    if isinstance(obj, dict):
        for key in ("pseudo_novel_or_isolation_tiers", "pseudo_novel_tiers", "pseudo_tiers"):
            if isinstance(obj.get(key), dict):
                tiers = obj[key]
                break
    if isinstance(tiers, dict) and tiers:
        vals: list[float] = []
        for value in tiers.values():
            if isinstance(value, dict):
                x = None
                for k in ("mean_recovery_gain_percentage_points", "mean_gain_pp", "gain_pp"):
                    if k in value:
                        x = scalar(value[k])
                        break
                if x is not None:
                    vals.append(x)
        if vals:
            return all(x > 0 for x in vals)
    gates = None
    if isinstance(obj, dict) and isinstance(obj.get("gates"), dict):
        gates = obj["gates"]
    if gates:
        for k, v in gates.items():
            if "pseudo" in k.lower():
                return bool(v)
    return None


def status_from_obj(obj: Any, flat: dict[str, Any]) -> str:
    if isinstance(obj, dict):
        for key in ("overall_status", "status", "verdict", "decision"):
            if isinstance(obj.get(key), str):
                return obj[key]
    for key, value in flat.items():
        if key.endswith(".status") and isinstance(value, str):
            return value
    return "UNKNOWN"


def leakage_ok(path: Path, obj: Any, flat: dict[str, Any]) -> bool:
    bad = bool_lookup(flat, ["truth_leakage", "leakage_detected", "gold_used_for_construction"])
    if bad is True:
        return False
    no_leak = bool_lookup(flat, ["no_truth_leakage", "leakage_free", "truth_blind"])
    if no_leak is False:
        return False
    text = path.as_posix().lower() + " " + json.dumps(obj, ensure_ascii=False).lower()[:100000]
    if "truth_accessed\": true" in text and "after" not in text and "evaluation" not in text:
        return False
    return True


@dataclass
class Candidate:
    source: str
    source_sha256: str
    original_status: str
    is_holdout: bool
    is_development: bool
    leakage_ok: bool
    genome_fraction_gain_pp: float | None
    mean_recovery_gain_pp: float | None
    low_abundance_gain_pp: float | None
    relative_chimera_change: float | None
    longitudinal_drop_pp: float | None
    paired_bootstrap_ci_low_pp: float | None
    paired_bootstrap_ci_high_pp: float | None
    pseudo_all_positive: bool | None
    metrics_complete: bool
    gate_genome_fraction: bool
    gate_mean_recovery: bool
    gate_low_abundance: bool
    gate_chimera: bool
    gate_longitudinal: bool
    gate_bootstrap: bool
    gate_pseudo: bool
    strict_holdout_success: bool


def build_candidate(path: Path, obj: Any) -> Candidate:
    flat = flatten(obj)
    gf = lookup(flat, ALIASES["gf_gain"])
    mean = lookup(flat, ALIASES["mean_gain"])
    low = lookup(flat, ALIASES["low_gain"])
    chim = lookup(flat, ALIASES["chimera_change"])
    longdrop = lookup(flat, ALIASES["longitudinal_drop"])
    ci_lo, ci_hi = ci_from_obj(obj, flat)
    pseudo = pseudo_positive(obj, flat)
    text = path.as_posix()
    is_development = bool(DEVELOPMENT_WORD.search(text)) and not bool(HOLDOUT_WORD.search(text))
    is_holdout = bool(HOLDOUT_WORD.search(text)) and not is_development
    complete = all(v is not None for v in (gf, mean, low, chim, longdrop, ci_lo)) and pseudo is not None
    g_gf = gf is not None and gf >= THRESHOLDS["minimum_genome_fraction_gain_pp"]
    g_mean = mean is not None and mean >= THRESHOLDS["minimum_mean_genome_recovery_gain_pp"]
    g_low = low is not None and low >= THRESHOLDS["minimum_low_abundance_gain_pp"]
    g_chim = chim is not None and chim <= THRESHOLDS["maximum_relative_chimera_increase"]
    g_long = longdrop is not None and longdrop >= THRESHOLDS["minimum_longitudinal_ablation_drop_pp"]
    g_ci = ci_lo is not None and ci_lo > THRESHOLDS["minimum_paired_bootstrap_ci_low_pp"]
    g_pseudo = pseudo is True
    leak = leakage_ok(path, obj, flat)
    strict = bool(is_holdout and complete and leak and g_gf and g_mean and g_low and g_chim and g_long and g_ci and g_pseudo)
    return Candidate(
        source=path.as_posix(), source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        original_status=status_from_obj(obj, flat), is_holdout=is_holdout,
        is_development=is_development, leakage_ok=leak,
        genome_fraction_gain_pp=gf, mean_recovery_gain_pp=mean,
        low_abundance_gain_pp=low, relative_chimera_change=chim,
        longitudinal_drop_pp=longdrop, paired_bootstrap_ci_low_pp=ci_lo,
        paired_bootstrap_ci_high_pp=ci_hi, pseudo_all_positive=pseudo,
        metrics_complete=complete, gate_genome_fraction=g_gf,
        gate_mean_recovery=g_mean, gate_low_abundance=g_low,
        gate_chimera=g_chim, gate_longitudinal=g_long,
        gate_bootstrap=g_ci, gate_pseudo=g_pseudo,
        strict_holdout_success=strict,
    )


def load_runs(path: Path) -> list[dict[str, Any]]:
    obj = json.loads(path.read_text())
    if isinstance(obj, dict):
        for key in ("workflow_runs", "runs"):
            if isinstance(obj.get(key), list):
                return obj[key]
    return obj if isinstance(obj, list) else []


def critical_coverage(runs: list[dict[str, Any]]) -> tuple[dict[str, bool], list[dict[str, Any]]]:
    categories = {
        "read_only_pairing": re.compile(r"pairing|pair inference|mash.*pair", re.I),
        "additive_holdout": re.compile(r"additive.*holdout", re.I),
        "hybrid_holdout": re.compile(r"v4.*holdout|hybrid.*holdout", re.I),
        "precision_holdout": re.compile(r"v5.*holdout|precision.*holdout", re.I),
        "extension_holdout": re.compile(r"v6.*extension.*holdout|extension.*holdout", re.I),
        "high_depth_long": re.compile(r"high.depth.*(metaflye|metamdbg)|metamdbg.*holdout", re.I),
        "domain_specific": re.compile(r"domain.specific|viral|plasmid", re.I),
    }
    latest: dict[tuple[str, str], dict[str, Any]] = {}
    for run in runs:
        name = str(run.get("name", ""))
        branch = str(run.get("head_branch", ""))
        key = (name, branch)
        old = latest.get(key)
        stamp = str(run.get("created_at", ""))
        if old is None or stamp > str(old.get("created_at", "")):
            latest[key] = run
    latest_runs = list(latest.values())
    coverage: dict[str, bool] = {}
    for category, regex in categories.items():
        matching = [r for r in latest_runs if regex.search(str(r.get("name", "")))]
        coverage[category] = any(r.get("status") == "completed" for r in matching)
    return coverage, latest_runs


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    a = parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    candidates: list[Candidate] = []
    parse_errors: list[dict[str, str]] = []
    for path in sorted(a.evidence_root.rglob("*.json")):
        if not DECISION_NAME.search(path.name) and not DECISION_NAME.search(path.parent.name):
            continue
        if path.stat().st_size > 20_000_000:
            continue
        try:
            obj = json.loads(path.read_text(encoding="utf-8", errors="strict"))
            if isinstance(obj, (dict, list)):
                candidates.append(build_candidate(path, obj))
        except Exception as exc:
            parse_errors.append({"path": path.as_posix(), "error": repr(exc)})

    runs = load_runs(a.runs_json)
    coverage, latest_runs = critical_coverage(runs)
    active = [r for r in latest_runs if r.get("status") in {"queued", "in_progress", "waiting", "pending", "requested"}]
    technical_failures = [r for r in latest_runs if r.get("status") == "completed" and r.get("conclusion") in {"failure", "timed_out", "cancelled", "action_required", "startup_failure"}]
    successes = [c for c in candidates if c.strict_holdout_success]
    complete_paths = all(coverage.values())

    challenge_ready = False
    blind_applied = False
    for path in sorted(a.evidence_root.rglob("*.json")):
        if path.stat().st_size > 5_000_000:
            continue
        try:
            obj = json.loads(path.read_text(encoding="utf-8", errors="replace"))
            flat = flatten(obj)
            if bool_lookup(flat, ["submission_ready", "cami_submission_ready"]) is True:
                challenge_ready = True
            if bool_lookup(flat, ["blind_applied", "challenge_data_applied", "restricted_data_processed"]) is True:
                blind_applied = True
        except Exception:
            pass

    if successes and challenge_ready and blind_applied:
        overall = "SUCCESS"
    elif successes:
        overall = "TOY_STRICT_SUCCESS_AWAITING_RESTRICTED_CAMI_APPLICATION"
    elif active:
        overall = "INCOMPLETE_ACTIVE_RUNS"
    elif technical_failures or not complete_paths:
        overall = "INCOMPLETE_TECHNICAL_OR_PATH_COVERAGE"
    else:
        overall = "FINAL_NEGATIVE_RESULT"

    def rank(c: Candidate) -> tuple[Any, ...]:
        passed = sum([
            c.gate_genome_fraction, c.gate_mean_recovery, c.gate_low_abundance,
            c.gate_chimera, c.gate_longitudinal, c.gate_bootstrap, c.gate_pseudo,
        ])
        return (c.strict_holdout_success, c.is_holdout, c.metrics_complete, passed,
                c.genome_fraction_gain_pp or -1e9, c.mean_recovery_gain_pp or -1e9)

    candidates.sort(key=rank, reverse=True)
    best = candidates[0] if candidates else None
    machine = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "overall_status": overall,
        "thresholds": THRESHOLDS,
        "strict_success_count": len(successes),
        "best_candidate": asdict(best) if best else None,
        "critical_path_coverage": coverage,
        "all_critical_paths_completed": complete_paths,
        "active_run_count": len(active),
        "technical_failure_count": len(technical_failures),
        "challenge_submission_ready": challenge_ready,
        "restricted_blind_data_applied": blind_applied,
        "n_decision_candidates": len(candidates),
        "n_parse_errors": len(parse_errors),
        "claim_boundary": (
            "SUCCESS requires every original numerical gate on an untouched holdout, no truth leakage, "
            "and actual restricted CAMI application plus submission readiness. Toy-only success is labeled separately."
        ),
    }
    (a.out / "FINAL_MACHINE_SUMMARY.json").write_text(json.dumps(machine, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_csv(a.out / "DECISION_CANDIDATES.csv", [asdict(c) for c in candidates])
    write_csv(a.out / "LATEST_WORKFLOW_RUNS.csv", latest_runs)
    write_csv(a.out / "TECHNICAL_FAILURES.csv", technical_failures)
    write_csv(a.out / "ACTIVE_RUNS.csv", active)
    write_csv(a.out / "PARSE_ERRORS.csv", parse_errors)

    if overall == "SUCCESS":
        opening = "연구 완료 — LANTERN이 엄격한 독립 holdout 관문을 통과했고 실제 CAMI restricted data 적용 및 submission-ready 상태까지 확인되었습니다."
    elif overall == "TOY_STRICT_SUCCESS_AWAITING_RESTRICTED_CAMI_APPLICATION":
        opening = "연구 완료 — 공개 Toy의 엄격한 미사용 holdout 성공 기준은 통과했으며, 실제 CAMI restricted data 적용과 최종 제출만 사용자 승인·접근 단계로 남았습니다."
    elif overall == "FINAL_NEGATIVE_RESULT":
        opening = "연구 완료 — 사전 정의된 합리적 계산 경로를 소진했으나 엄격한 성공 기준 전체를 재현 가능하게 통과하지 못했습니다."
    else:
        opening = "연구 종료 판정 보류 — 활성 계산 또는 기술적·경로완료 관문이 남아 있어 성공이나 최종 음성으로 판정할 수 없습니다."

    lines = [
        f"# {opening}", "",
        f"- 기계 판정: **{overall}**",
        f"- 엄격한 holdout 성공 후보: **{len(successes)}개**",
        f"- 분석된 decision/result JSON: **{len(candidates)}개**",
        f"- 활성 workflow: **{len(active)}개**",
        f"- 최신 기술 실패 workflow: **{len(technical_failures)}개**",
        f"- 모든 사전 정의 경로 완료: **{'예' if complete_paths else '아니오'}**",
        f"- CAMI submission-ready 확인: **{'예' if challenge_ready else '아니오'}**",
        f"- restricted blind data 적용 확인: **{'예' if blind_applied else '아니오'}**",
        "", "## 변경하지 않은 성공 기준", "",
    ]
    for key, value in THRESHOLDS.items():
        lines.append(f"- `{key}`: `{value}`")
    if best:
        lines += ["", "## 가장 강한 확인 후보", ""]
        for key, value in asdict(best).items():
            lines.append(f"- `{key}`: `{value}`")
    lines += ["", "## 경로 완료 감사", ""]
    for key, value in coverage.items():
        lines.append(f"- `{key}`: **{'완료' if value else '미완료'}**")
    lines += ["", "## 주장 경계", "", machine["claim_boundary"], ""]
    (a.out / "FINAL_VERDICT.md").write_text("\n".join(lines), encoding="utf-8")

    manifest = {
        "runs_json_sha256": hashlib.sha256(a.runs_json.read_bytes()).hexdigest(),
        "artifacts_json_sha256": hashlib.sha256(a.artifacts_json.read_bytes()).hexdigest(),
        "evidence_root": a.evidence_root.as_posix(),
        "candidate_sources": [c.source for c in candidates],
    }
    (a.out / "RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(machine, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
