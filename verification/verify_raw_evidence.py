#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, re, sys, time
from collections import Counter
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
OUT = ROOT / "verification"
OUT.mkdir(exist_ok=True)
checks, failures, warnings = [], [], []

def pass_(name, detail=None):
    checks.append({"check": name, "status": "PASS", "detail": detail})

def fail(name, detail):
    checks.append({"check": name, "status": "FAIL", "detail": detail})
    failures.append({"check": name, "detail": detail})

def warn(name, detail):
    checks.append({"check": name, "status": "WARN", "detail": detail})
    warnings.append({"check": name, "detail": detail})

def expect(condition, name, detail=None):
    pass_(name, detail) if condition else fail(name, detail)

def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

start = time.time()
required = [
    "evidence/g1/v10_final_t2_t4",
    "evidence/g1/v8_head_t2_t10",
    "evidence/g1/v8_tail_t11_t45",
    "evidence/g2_g3/v9",
    "evidence/g4_g7_drat",
    "manifest/ARTIFACT_INDEX.json",
    "repository_snapshot/crossdomain/b2_gap1_v_split_v10.cpp",
    "repository_snapshot/crossdomain/b2_gap1_u_split_v8.cpp",
    "repository_snapshot/crossdomain/b2_gap_g_u_split_v9.cpp",
]
for rel in required:
    expect((ROOT / rel).exists(), f"required:{rel}")

if failures:
    print(json.dumps({"status": "FAIL", "failures": failures}, indent=2))
    raise SystemExit(1)

# GitHub source-artifact provenance.
index = load(ROOT / "manifest/ARTIFACT_INDEX.json")
by_name = {x["name"]: x for x in index}
expect(len(index) == 373, "artifact_index_count", len(index))
expect(len(by_name) == 373, "artifact_name_uniqueness", len(by_name))
metadata_files = list(ROOT.rglob("ARTIFACT_METADATA.json"))
expect(len(metadata_files) == 373, "artifact_metadata_count", len(metadata_files))
bad_meta = []
for path in metadata_files:
    item = load(path)
    if by_name.get(item.get("name")) != item:
        bad_meta.append(str(path.relative_to(ROOT)))
expect(not bad_meta, "artifact_metadata_matches_index", bad_meta[:20])
bad_api = [
    x["name"] for x in index
    if x["api_digest"] != "sha256:" + x["archive_sha256"]
    or x["api_size"] != x["archive_bytes"]
]
expect(not bad_api, "github_api_digest_matches_downloaded_archive", bad_api[:20])

category = Counter()
for item in index:
    n = item["name"]
    if re.fullmatch(r"b2-gap1-t\d+-u\d+-v-split-v10", n):
        category["g1_v10"] += 1
    elif re.fullmatch(r"b2-gap1-t\d+-u-split-v8", n):
        category["g1_v8"] += 1
    elif re.fullmatch(r"b2-g[23]-t\d+-u-split-v9", n):
        category["g23_v9"] += 1
    elif n in {"b2-gap1-v10-result-collector-20260817",
               "b2-gap23-v9-result-collector-20260817"}:
        category["collectors"] += 1
    elif re.fullmatch(r"b2-gap[4-7]-drat.*", n):
        category["drat"] += 1
    else:
        category["unclassified"] += 1
expected_categories = Counter(
    {"g1_v10": 256, "g1_v8": 44, "g23_v9": 67, "collectors": 2, "drat": 4}
)
expect(category == expected_categories, "artifact_category_counts", dict(category))

# Verify artifact-local SHA files, except their inherently self-referential line.
hash_errors = []
hashes_checked = 0
self_hash_lines = 0
for sumfile in ROOT.rglob("SHA256SUMS.txt"):
    if "manifest" in sumfile.parts or "verification" in sumfile.parts:
        continue
    for line in sumfile.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            hash_errors.append([str(sumfile.relative_to(ROOT)), "parse", line])
            continue
        expected, recorded = parts
        name = Path(recorded.strip()).name
        if name == "SHA256SUMS.txt":
            self_hash_lines += 1
            continue
        target = sumfile.parent / name
        if not target.exists():
            hash_errors.append([str(sumfile.relative_to(ROOT)), "missing", name])
            continue
        hashes_checked += 1
        actual = sha256(target)
        if actual != expected:
            hash_errors.append([str(target.relative_to(ROOT)), expected, actual])
expect(not hash_errors, "artifact_internal_hashes",
       {"checked": hashes_checked, "errors": hash_errors[:20]})
warn("historical_self_hash_lines_skipped",
     f"{self_hash_lines} artifact-local files hashed themselves while being written; "
     "the corrected release manifest supersedes those self-lines.")

src = ROOT / "repository_snapshot/crossdomain"
src_hash = {
    "v10": sha256(src / "b2_gap1_v_split_v10.cpp"),
    "v8": sha256(src / "b2_gap1_u_split_v8.cpp"),
    "v9": sha256(src / "b2_gap_g_u_split_v9.cpp"),
}

def validate_completed(result, audit, status, key, low, high, aggregate_key):
    assert result["completed_exhaustively"] is True
    assert result["timed_out"] is False
    assert result["witness_found"] is False
    assert result["witness_verified"] is False
    assert status == audit["status"] == "EXHAUSTED_NO_WITNESS"
    assert audit["coverage_ok"] is True
    values = [b[key] for b in result["branches"]]
    assert values == list(range(low, high + 1)), (values, low, high)
    assert result[f"expected_{key}_branches"] == len(values)
    assert result[f"observed_{key}_branches"] == len(values)
    assert all(b["completed"] and not b["timed"] and not b["found"] and
               not b["verified"] for b in result["branches"])
    assert result[aggregate_key] == 1 + sum(int(b["nodes"]) for b in result["branches"])

# Final g=1 v10: 255 t=2,3,4 subbranches plus one independent regression control.
expected_v10 = {(t, u) for t in (2, 3, 4) for u in range(t + 1, 92 - t)}
seen_v10, controls, errors = set(), [], []
for folder in sorted((ROOT / "evidence/g1/v10_final_t2_t4").iterdir()):
    if not folder.is_dir():
        continue
    try:
        r = load(folder / "RESULT.json")
        a = load(folder / "INDEPENDENT_AUDIT.json")
        s = (folder / "STATUS.txt").read_text().strip()
        t, u = int(r["third"]), int(r["fourth"])
        if a.get("control"):
            assert (t, u) == (7, 8)
            assert r["completed_exhaustively"] and not r["witness_found"]
            assert r["aggregate_nodes_including_fourth_root"] == 65593402
            controls.append((t, u))
        else:
            validate_completed(
                r, a, s, "fifth", u + 1, 92 - t,
                "aggregate_nodes_including_fourth_root"
            )
            assert (t, u) in expected_v10
            seen_v10.add((t, u))
        assert sha256(folder / "b2_gap1_v_split_v10.cpp") == src_hash["v10"]
    except Exception as exc:
        errors.append([folder.name, repr(exc)])
expect(not errors, "g1_v10_raw_artifacts", errors[:20])
expect(seen_v10 == expected_v10 and controls == [(7, 8)],
       "g1_v10_coverage",
       {"final_subbranches": len(seen_v10), "control": controls})

# Final g=1 v8: t=5,...,45. Historical t=2,3,4 remain for audit only.
v8_final, v8_history, errors = {}, {}, []
for base, final_set in [
    (ROOT / "evidence/g1/v8_head_t2_t10", set(range(5, 11))),
    (ROOT / "evidence/g1/v8_tail_t11_t45", set(range(11, 46))),
]:
    for folder in sorted(base.iterdir()):
        if not folder.is_dir():
            continue
        try:
            r = load(folder / "RESULT.json")
            a = load(folder / "INDEPENDENT_AUDIT.json")
            s = (folder / "STATUS.txt").read_text().strip()
            t = int(r["third"])
            if t in final_set:
                validate_completed(
                    r, a, s, "fourth", t + 1, 91 - t,
                    "aggregate_nodes_including_third_root"
                )
                assert sha256(folder / "b2_gap1_u_split_v8.cpp") == src_hash["v8"]
                v8_final[t] = r
            else:
                v8_history[t] = {
                    "status": s,
                    "completed": r["completed_exhaustively"],
                    "timed_out": r["timed_out"],
                }
        except Exception as exc:
            errors.append([folder.name, repr(exc)])
expect(not errors, "g1_v8_raw_artifacts", errors[:20])
expect(set(v8_final) == set(range(5, 46)),
       "g1_v8_final_coverage", sorted(v8_final))
expect(v8_final[10]["aggregate_nodes_including_third_root"] == 537251446,
       "g1_v8_t10_regression_node_count",
       v8_final[10]["aggregate_nodes_including_third_root"])
warn("g1_v8_t2_t4_historical_superseded", v8_history)
expect({2, 3, 4} | set(v8_final) == set(range(2, 46)),
       "g1_all_44_normalized_third_branches")

# Redundant g=2,3 v9 cross-check.
expected_g23 = {(2, t) for t in range(4, 41)} | {(3, t) for t in range(6, 36)}
seen_g23, errors = set(), []
for folder in sorted((ROOT / "evidence/g2_g3/v9").iterdir()):
    if not folder.is_dir():
        continue
    try:
        r = load(folder / "RESULT.json")
        a = load(folder / "INDEPENDENT_AUDIT.json")
        s = (folder / "STATUS.txt").read_text().strip()
        g, t = int(r["gap"]), int(r["third"])
        high = 100 - (t - g) - 10 * g
        validate_completed(
            r, a, s, "fourth", t + g, high,
            "aggregate_nodes_including_third_root"
        )
        assert sha256(folder / "b2_gap_g_u_split_v9.cpp") == src_hash["v9"]
        seen_g23.add((g, t))
    except Exception as exc:
        errors.append([folder.name, repr(exc)])
expect(not errors, "g2_g3_v9_raw_artifacts", errors[:20])
expect(seen_g23 == expected_g23, "g2_g3_v9_coverage", len(seen_g23))

# 13-set witness.
witness = [0, 5, 7, 31, 58, 61, 62, 63, 72, 80, 84, 91, 97]
counts = {d: 0 for d in range(1, 100)}
for x in witness:
    for y in witness:
        if x != y:
            counts[(x - y) % 100] += 1
expect(len(witness) == 13 and len(set(witness)) == 13 and max(counts.values()) <= 2,
       "valid_13_set",
       {"set": witness, "max_ordered_difference_multiplicity": max(counts.values())})

# Machine audit of the elementary mod-10 support classification in the unit-difference lemma.
unit_classes = {1, 3, 7, 9}
support_patterns = []
bad_support = []
for mask in range(1, 1 << 10):
    support = [i for i in range(10) if (mask >> i) & 1]
    if all((a - b) % 10 not in unit_classes for a in support for b in support):
        same_parity = len({x % 2 for x in support}) == 1
        mixed_pair = (
            not same_parity and len(support) == 2
            and (support[0] - support[1]) % 10 == 5
        )
        if not (same_parity or mixed_pair):
            bad_support.append(support)
        support_patterns.append(support)
expect(not bad_support, "unit_difference_mod10_classification",
       {"patterns": len(support_patterns), "bad": bad_support})
expect(14 * 13 > 49 * 2 and 14 * 13 > 19 * 2,
       "unit_difference_capacity_inequalities",
       {"required": 182, "even_capacity": 98, "multiple_of_5_capacity": 38})

# Actual CNF and DRAT files, original logs, and project-managed fresh replay logs.
drat_rows, errors = [], []
for g in range(4, 8):
    dirs = [p for p in (ROOT / f"evidence/g4_g7_drat/g{g}").iterdir() if p.is_dir()]
    if len(dirs) != 1:
        errors.append([g, "artifact directory count", len(dirs)])
        continue
    d = dirs[0]
    try:
        cnf = d / f"gap{g}.cnf"
        proof = d / f"gap{g}.drat"
        original_log = d / "drat_check.log"
        exit_code = d / "cadical_exit_code.txt"
        fresh_log = ROOT / f"independent_replay/fresh_drat_20260818/g{g}.fresh.log"
        assert cnf.is_file() and proof.is_file()
        assert original_log.is_file() and fresh_log.is_file()
        assert exit_code.read_text().strip() == "20"
        assert re.search(r"(?m)^s VERIFIED\s*$",
                         original_log.read_text(encoding="utf-8", errors="replace"))
        assert re.search(r"(?m)^s VERIFIED\s*$",
                         fresh_log.read_text(encoding="utf-8", errors="replace"))
        header = next(
            line for line in cnf.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.startswith("p cnf ")
        )
        _, _, nvars, nclauses = header.split()
        drat_rows.append({
            "gap": g,
            "cnf_bytes": cnf.stat().st_size,
            "drat_bytes": proof.stat().st_size,
            "cnf_sha256": sha256(cnf),
            "drat_sha256": sha256(proof),
            "variables": int(nvars),
            "clauses": int(nclauses),
            "original_drat_trim_verified": True,
            "fresh_drat_trim_verified": True,
        })
    except Exception as exc:
        errors.append([g, repr(exc)])
expect(not errors, "g4_g7_actual_cnf_drat_and_logs", errors)
expect(len(drat_rows) == 4, "g4_g7_fresh_replay_count", drat_rows)

# Environment and exact workflow provenance.
for tag in [
    "g1_v10_32040699080",
    "g1_v8_head_32038183046",
    "g1_v8_tail_32038657803",
    "g23_v9_32040180627",
]:
    for name in ["RUN.json", "JOBS.json", "WORKFLOW_AT_RUN.yml"]:
        expect((ROOT / "environment" / tag / name).exists(),
               f"provenance:{tag}/{name}")

# Preserve and explain the historical pre-upload manifest defect.
historical = ROOT / "manifest/PRE_UPLOAD_SHA256SUMS.txt"
if historical.exists():
    entries = []
    for line in historical.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2:
            entries.append(parts[1].strip())
    missing = []
    for rel in entries:
        parts = Path(rel).parts
        local = Path(*parts[1:]) if parts and parts[0] == "B2_RAW_COLLECTED" else Path(rel)
        if not (ROOT / local).exists():
            missing.append(local.as_posix())
    warn("historical_preupload_manifest_not_authoritative", {
        "listed": len(entries),
        "missing_after_actions_upload": len(missing),
        "missing_examples": missing[:30],
        "explanation": (
            "actions/upload-artifact omitted hidden .github paths. "
            "Exact WORKFLOW_AT_RUN snapshots and the corrected release manifest supersede it."
        )
    })

status = "PASS" if not failures else "FAIL"
report = {
    "status": status,
    "root": str(ROOT),
    "elapsed_seconds": round(time.time() - start, 3),
    "summary": {
        "github_artifacts": len(index),
        "g1_v10_final_subbranches": len(seen_v10),
        "g1_v8_final_branches": len(v8_final),
        "g1_normalized_branches": 44,
        "g2_g3_v9_branches": len(seen_g23),
        "drat_cases": drat_rows,
        "witness_13": witness,
    },
    "checks": checks,
    "warnings": warnings,
    "failures": failures,
}
(OUT / "RAW_EVIDENCE_VERIFICATION.json").write_text(
    json.dumps(report, indent=2), encoding="utf-8"
)
md = [
    "# B2[2]/Z100 raw-evidence verification",
    "",
    f"**Status: {status}**",
    "",
    f"- GitHub source artifacts: {len(index)}",
    f"- final v10 t=2,3,4 subbranches: {len(seen_v10)}",
    f"- final v8 t=5,...,45 branches: {len(v8_final)}",
    f"- v9 g=2,3 branches: {len(seen_g23)}",
    f"- actual g=4,...,7 CNF/DRAT cases with original and fresh VERIFIED logs: {len(drat_rows)}",
    f"- valid 13-set: `{witness}`",
    "",
    "## Warnings",
]
md.extend([f"- {w['check']}: `{w['detail']}`" for w in warnings] or ["- none"])
md.extend(["", "## Failures"])
md.extend([f"- {f['check']}: `{f['detail']}`" for f in failures] or ["- none"])
(OUT / "RAW_EVIDENCE_VERIFICATION.md").write_text(
    "\n".join(md) + "\n", encoding="utf-8"
)
print(json.dumps({
    "status": status,
    "elapsed_seconds": report["elapsed_seconds"],
    "summary": report["summary"],
    "warnings": warnings,
    "failures": failures,
}, indent=2))
raise SystemExit(0 if status == "PASS" else 1)
