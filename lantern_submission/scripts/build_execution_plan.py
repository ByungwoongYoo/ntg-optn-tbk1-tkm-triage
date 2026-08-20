#!/usr/bin/env python3
"""Create a deterministic active-CAMI execution plan from frozen metadata."""
from __future__ import annotations

import argparse
import hashlib
import json
import shlex
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--mapping-freeze", type=Path, required=True)
    p.add_argument("--input-freeze", type=Path, required=True)
    p.add_argument("--rule", type=Path, required=True)
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--memory-gb", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    if a.threads <= 0 or a.memory_gb <= 0:
        raise ValueError("threads and memory must be positive")
    mapping = json.loads(a.mapping_freeze.read_text(encoding="utf-8"))
    inputs = json.loads(a.input_freeze.read_text(encoding="utf-8"))
    rule = json.loads(a.rule.read_text(encoding="utf-8"))
    if mapping.get("mapping_method") != "explicit_metadata_only":
        raise ValueError("non-explicit mapping rejected")
    if inputs.get("status") != "ACTIVE_INPUT_MANIFEST_FROZEN":
        raise ValueError("input manifest is not frozen")
    if rule.get("post_result_tuning_allowed") is not False:
        raise ValueError("rule does not prohibit post-result tuning")
    if set(mapping["sample_ids"]) != set(inputs["sample_ids"]):
        raise ValueError("mapping and input sample sets differ")

    runner = Path("lantern_submission/scripts/run_individual_frozen_v7.py")
    commands = []
    for group in mapping["groups"]:
        individual = str(group["individual_id"])
        command = [
            "python", str(runner),
            "--individual-id", individual,
            "--mapping-freeze", str(a.mapping_freeze),
            "--input-freeze", str(a.input_freeze),
            "--rule", str(a.rule),
            "--out", str(a.work_dir),
            "--threads", str(a.threads),
            "--memory-gb", str(a.memory_gb),
        ]
        commands.append({
            "individual_id": individual,
            "samples": list(group["samples"]),
            "command": command,
            "shell": shlex.join(command),
        })

    a.out.mkdir(parents=True, exist_ok=True)
    plan = {
        "status": "PREREGISTERED_EXECUTION_PLAN",
        "candidate": rule["candidate"],
        "algorithm_source_commit": rule["algorithm_source_commit"],
        "threads_per_individual": a.threads,
        "memory_gb_per_individual": a.memory_gb,
        "work_dir": str(a.work_dir),
        "n_individuals": len(commands),
        "n_samples": len(mapping["sample_ids"]),
        "mapping_freeze_sha256": sha256_file(a.mapping_freeze),
        "input_freeze_sha256": sha256_file(a.input_freeze),
        "rule_sha256": sha256_file(a.rule),
        "commands": commands,
        "parallelism": "not prespecified; scheduler may serialize or parallelize without changing commands",
        "result_dependent_changes_allowed": False,
        "restricted_storage": "private local/HPC only",
    }
    json_path = a.out / "RUN_PLAN.json"
    json_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    shell_path = a.out / "RUN_PLAN.sh"
    shell_path.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n\n" +
        "\n".join(item["shell"] for item in commands) + "\n",
        encoding="utf-8",
    )
    shell_path.chmod(0o755)
    (a.out / "MANIFEST.sha256").write_text(
        f"{sha256_file(a.mapping_freeze)}  {a.mapping_freeze}\n"
        f"{sha256_file(a.input_freeze)}  {a.input_freeze}\n"
        f"{sha256_file(a.rule)}  {a.rule}\n"
        f"{sha256_file(json_path)}  {json_path.name}\n"
        f"{sha256_file(shell_path)}  {shell_path.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(plan, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
