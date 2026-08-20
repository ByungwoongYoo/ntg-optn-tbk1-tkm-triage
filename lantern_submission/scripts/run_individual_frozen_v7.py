#!/usr/bin/env python3
"""Execute the unchanged frozen LANTERN-v7 rule for one explicit individual group.

Restricted reads and all derived files must remain on private local/HPC storage.
This runner never uploads data and never changes the frozen rule in response to output.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--individual-id", required=True)
    p.add_argument("--mapping-freeze", type=Path, required=True)
    p.add_argument("--input-freeze", type=Path, required=True)
    p.add_argument("--rule", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--memory-gb", type=int, required=True)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"required executable not found: {name}")
    return path


def run(command: list[str], log: Path, cwd: Path | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\t" + shlex.join(command) + "\n")
        handle.flush()
        completed = subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, text=True)
        handle.write(f"\nEXIT_CODE\t{completed.returncode}\nWALL_SECONDS\t{time.time()-started:.3f}\n")
    if completed.returncode != 0:
        raise RuntimeError(f"command failed ({completed.returncode}); see {log}")


def merge_gzip(inputs: list[Path], output: Path, threads: int, log: Path) -> None:
    pigz = require_tool("pigz")
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with log.open("w", encoding="utf-8") as handle, output.open("wb") as out_handle:
        handle.write("INPUTS\t" + "\t".join(str(path) for path in inputs) + "\n")
        decompressor = subprocess.Popen([pigz, "-dc", *map(str, inputs)], stdout=subprocess.PIPE, stderr=handle)
        assert decompressor.stdout is not None
        compressor = subprocess.Popen([pigz, "-p", str(max(1, threads))], stdin=decompressor.stdout, stdout=out_handle, stderr=handle)
        decompressor.stdout.close()
        compressor_rc = compressor.wait()
        decompressor_rc = decompressor.wait()
        handle.write(
            f"DECOMPRESS_EXIT\t{decompressor_rc}\nCOMPRESS_EXIT\t{compressor_rc}\n"
            f"WALL_SECONDS\t{time.time()-started:.3f}\n"
        )
    if decompressor_rc != 0 or compressor_rc != 0 or output.stat().st_size == 0:
        raise RuntimeError(f"failed to materialize merged gzip: {output}")
    with gzip.open(output, "rb") as handle:
        while handle.read(8 * 1024 * 1024):
            pass


def choose_spades_executable() -> str:
    for name in ("metaspades.py", "spades.py"):
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError("neither metaspades.py nor spades.py is available")


def main() -> None:
    a = parse_args()
    if a.threads <= 0 or a.memory_gb <= 0:
        raise ValueError("threads and memory must be positive")
    mapping = json.loads(a.mapping_freeze.read_text(encoding="utf-8"))
    inputs = json.loads(a.input_freeze.read_text(encoding="utf-8"))
    rule = json.loads(a.rule.read_text(encoding="utf-8"))
    if mapping.get("mapping_method") != "explicit_metadata_only":
        raise ValueError("mapping is not explicit-metadata-only")
    if rule.get("candidate") != "LANTERN-v7-short-terminal-500-s1":
        raise ValueError("unexpected candidate rule")
    group = next((g for g in mapping["groups"] if str(g["individual_id"]) == str(a.individual_id)), None)
    if group is None:
        raise ValueError(f"individual absent from mapping: {a.individual_id}")
    input_by_sample = {str(row["sample_id"]): row for row in inputs["records"]}
    samples = [str(sample) for sample in group["samples"]]
    if any(sample not in input_by_sample for sample in samples):
        raise ValueError("mapped sample absent from input freeze")

    root = a.out / f"individual_{a.individual_id}"
    dirs = {name: root / name for name in ("reads", "assemblies", "candidates", "extension", "logs", "provenance", "submission")}
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    plan = {
        "status": "DRY_RUN" if a.dry_run else "EXECUTING",
        "individual_id": str(a.individual_id),
        "sample_ids": samples,
        "threads": a.threads,
        "memory_gb": a.memory_gb,
        "mapping_freeze_sha256": sha256_file(a.mapping_freeze),
        "input_freeze_sha256": sha256_file(a.input_freeze),
        "rule_sha256": sha256_file(a.rule),
        "algorithm_source_commit": rule["algorithm_source_commit"],
        "post_result_tuning_allowed": False,
    }
    (dirs["provenance"] / "INDIVIDUAL_PLAN.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    if a.dry_run:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return

    for tool in ("megahit", "flye", "minimap2", "pigz"):
        require_tool(tool)
    spades = choose_spades_executable()

    r1_files = [Path(input_by_sample[s]["files"]["short_r1"]["path"]) for s in samples]
    r2_files = [Path(input_by_sample[s]["files"]["short_r2"]["path"]) for s in samples]
    long_files = [Path(input_by_sample[s]["files"]["long_reads"]["path"]) for s in samples]
    r1 = dirs["reads"] / "combined_R1.fastq.gz"
    r2 = dirs["reads"] / "combined_R2.fastq.gz"
    long_reads = dirs["reads"] / "combined_long.fastq.gz"
    merge_gzip(r1_files, r1, a.threads, dirs["logs"] / "merge_R1.log")
    merge_gzip(r2_files, r2, a.threads, dirs["logs"] / "merge_R2.log")
    merge_gzip(long_files, long_reads, a.threads, dirs["logs"] / "merge_long.log")

    megahit_dir = dirs["assemblies"] / "megahit_short"
    run([
        "megahit", "-1", str(r1), "-2", str(r2), "-o", str(megahit_dir),
        "-t", str(a.threads), "--min-contig-len", "500",
    ], dirs["logs"] / "megahit.log")
    megahit_fasta = megahit_dir / "final.contigs.fa"

    spades_short_dir = dirs["assemblies"] / "metaspades_short"
    run([
        spades, "--meta", "-1", str(r1), "-2", str(r2), "-o", str(spades_short_dir),
        "-t", str(a.threads), "-m", str(a.memory_gb),
    ], dirs["logs"] / "metaspades_short.log")
    short_fasta = spades_short_dir / "contigs.fasta"

    spades_hybrid_dir = dirs["assemblies"] / "metaspades_hybrid"
    run([
        spades, "--meta", "-1", str(r1), "-2", str(r2), "--nanopore", str(long_reads),
        "-o", str(spades_hybrid_dir), "-t", str(a.threads), "-m", str(a.memory_gb),
    ], dirs["logs"] / "metaspades_hybrid.log")
    hybrid_fasta = spades_hybrid_dir / "contigs.fasta"

    flye_dir = dirs["assemblies"] / "flye_long"
    run([
        "flye", "--nano-raw", str(long_reads), "--meta", "--out-dir", str(flye_dir),
        "--threads", str(a.threads),
    ], dirs["logs"] / "flye.log")
    flye_fasta = flye_dir / "assembly.fasta"

    for path in (megahit_fasta, short_fasta, hybrid_fasta, flye_fasta):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"expected assembly missing: {path}")

    repo_root = Path(__file__).resolve().parents[2]
    prepare = repo_root / "lantern_v6" / "prepare_extension_candidates.py"
    extension = repo_root / "lantern_v6" / "extension_union.py"
    validator = repo_root / "lantern_cami3" / "scripts" / "verify_cami_assembly.py"
    run([
        sys.executable, str(prepare),
        "--source", f"metaspades_hybrid={hybrid_fasta}",
        "--source", f"megahit_short={megahit_fasta}",
        "--source", f"flye_long={flye_fasta}",
        "--min-length", "500", "--max-n-fraction", "0.05",
        "--out", str(dirs["candidates"]),
    ], dirs["logs"] / "prepare_candidates.log")

    paf = dirs["candidates"] / "candidate_to_backbone.paf"
    with paf.open("wb") as paf_handle, (dirs["logs"] / "candidate_to_backbone.log").open("wb") as err_handle:
        completed = subprocess.run([
            "minimap2", "-x", "asm5", "-c", "-N", "100", "-t", str(a.threads),
            str(short_fasta), str(dirs["candidates"] / "extension_candidates.fasta"),
        ], stdout=paf_handle, stderr=err_handle)
    if completed.returncode != 0:
        raise RuntimeError("candidate-to-backbone minimap2 failed")

    params = rule["extension_parameters"]
    command = [
        sys.executable, str(extension),
        "--backbone", str(short_fasta),
        "--candidates", str(dirs["candidates"] / "extension_candidates.fasta"),
        "--metadata", str(dirs["candidates"] / "extension_candidate_metadata.tsv"),
        "--candidate-to-backbone-paf", str(paf),
        "--mode", str(params["mode"]),
        "--min-candidate-length", str(params["min_candidate_length"]),
        "--min-terminal-bp", str(params["min_terminal_bp"]),
        "--alignment-identity", str(params["alignment_identity"]),
        "--alignment-min-bp", str(params["alignment_min_bp"]),
        "--whole-max-aligned-fraction", str(params["whole_max_aligned_fraction"]),
        "--containment-fraction", str(params["containment_fraction"]),
        "--min-source-count", str(params["min_source_count"]),
        "--out", str(dirs["extension"]),
    ]
    if params["single_target_terminal"]:
        command.append("--single-target-terminal")
    run(command, dirs["logs"] / "extension_union.log")

    output_fasta = dirs["extension"] / "LANTERN_V6_EXTENSION_UNION.fasta"
    validated = dirs["submission"] / f"LANTERN_V7_INDIVIDUAL_{a.individual_id}.fasta"
    shutil.copyfile(output_fasta, validated)
    run([
        sys.executable, str(validator), str(validated),
        "--out", str(dirs["submission"] / "VALIDATION.json"),
    ], dirs["logs"] / "validate_fasta.log")
    (dirs["submission"] / "SHA256.txt").write_text(
        f"{sha256_file(validated)}  {validated.name}\n", encoding="utf-8"
    )
    plan["status"] = "COMPLETED"
    plan["output_fasta"] = str(validated.resolve())
    plan["output_sha256"] = sha256_file(validated)
    (dirs["provenance"] / "INDIVIDUAL_PLAN.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps(plan, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
