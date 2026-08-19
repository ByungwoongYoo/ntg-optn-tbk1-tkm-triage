#!/usr/bin/env python3
"""Map statistically replicated unitigs to representative assemblies and extract context.

Exact sequence mapping is performed in both orientations. The output supports annotation and
known-mechanism review; it does not assign function or novelty by itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--strict", required=True); p.add_argument("--manifest", required=True); p.add_argument("--all-rtab", required=True)
    p.add_argument("--reference-gbff", required=True); p.add_argument("--out", default="unitig_contexts"); p.add_argument("--flank", type=int, default=10000); p.add_argument("--max-representatives", type=int, default=5)
    return p.parse_args()


def locate(seq: str, record: Any) -> tuple[int, int, str] | None:
    s = str(record.seq).upper(); q = seq.upper(); i = s.find(q)
    if i >= 0: return i, i + len(q), "+"
    rc = str(Seq(q).reverse_complement()); i = s.find(rc)
    if i >= 0: return i, i + len(q), "-"
    return None


def nearest_features(rec: Any, start: int, end: int, limit: int = 8) -> list[dict[str, Any]]:
    out = []
    for f in rec.features:
        if f.type not in {"CDS", "gene", "rRNA", "tRNA", "ncRNA"}: continue
        fs, fe = int(f.location.start), int(f.location.end)
        distance = 0 if not (fe <= start or fs >= end) else min(abs(fe - start), abs(fs - end))
        if distance > 20000: continue
        out.append({"type": f.type, "start": fs, "end": fe, "strand": int(f.location.strand or 0), "distance": distance, "gene": ";".join(f.qualifiers.get("gene", [])), "locus_tag": ";".join(f.qualifiers.get("locus_tag", [])), "product": ";".join(f.qualifiers.get("product", [])), "note": ";".join(f.qualifiers.get("note", []))[:500]})
    return sorted(out, key=lambda x: (x["distance"], x["type"] != "CDS", x["start"]))[:limit]


def main() -> None:
    a = args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True); contexts = out / "contexts"; contexts.mkdir(exist_ok=True)
    strict = pd.read_csv(a.strict, dtype={"candidate_id": str, "canonical_sequence": str, "variant": str}); seq_col = "canonical_sequence" if "canonical_sequence" in strict.columns else "variant"; strict[seq_col] = strict[seq_col].astype(str).str.upper()
    manifest = pd.read_csv(a.manifest, dtype={"assembly_ID": str}); rtab = pd.read_csv(a.all_rtab, sep="\t", index_col=0); rtab.index = rtab.index.astype(str).str.upper(); rtab.columns = rtab.columns.astype(str); rtab = rtab.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    refs = list(SeqIO.parse(a.reference_gbff, "genbank")); ref_rows = []; mappings = []; context_records = []
    for _, cand in strict.iterrows():
        cid = str(cand.candidate_id); seq = str(cand[seq_col]);
        for rec in refs:
            hit = locate(seq, rec)
            if hit:
                st, en, strand = hit; ref_rows.append({"candidate_id": cid, "record": rec.id, "start": st, "end": en, "strand": strand, "nearby_features": json.dumps(nearest_features(rec, st, en), ensure_ascii=False)})
        if seq not in rtab.index: continue
        present = [x for x in rtab.columns[rtab.loc[seq].astype(int) == 1].tolist() if x in set(manifest.assembly_ID)]
        sub = manifest[manifest.assembly_ID.isin(present)].copy(); sub["is_R"] = sub.phenotype.astype(str) == "R"; sort_cols = [c for c in ["is_R", "source_group", "Kleborate_ST", "assembly_ID"] if c in sub.columns]; sub = sub.sort_values(sort_cols, ascending=[False] + [True] * (len(sort_cols) - 1))
        chosen = []; seen_groups = set()
        for _, row in sub.iterrows():
            g = str(row.get("source_group", ""))
            if g and g in seen_groups and len(chosen) < max(2, a.max_representatives // 2): continue
            chosen.append(row); seen_groups.add(g)
            if len(chosen) >= a.max_representatives: break
        for j, row in enumerate(chosen, 1):
            p = Path(str(row.assembly_path)); found = False
            for rec in SeqIO.parse(str(p), "fasta"):
                hit = locate(seq, rec)
                if not hit: continue
                st, en, strand = hit; left = max(0, st - a.flank); right = min(len(rec.seq), en + a.flank); subrec = rec[left:right]; context_id = f"{cid}__{row.assembly_ID}__rep{j}"; subrec.id = context_id; subrec.name = context_id; subrec.description = f"candidate={cid} assembly={row.assembly_ID} phenotype={row.phenotype} source_group={row.get('source_group','')} contig={rec.id} interval={left+1}-{right} unitig={st+1}-{en} strand={strand}"; context_records.append(subrec)
                mappings.append({"candidate_id": cid, "sequence": seq, "assembly_ID": row.assembly_ID, "phenotype": row.phenotype, "source_group": row.get("source_group"), "Kleborate_ST": row.get("Kleborate_ST"), "assembly_path": str(p), "contig": rec.id, "start_1based": st + 1, "end_1based": en, "strand": strand, "context_id": context_id, "context_left_1based": left + 1, "context_right_1based": right}); found = True; break
            if not found: mappings.append({"candidate_id": cid, "sequence": seq, "assembly_ID": row.assembly_ID, "phenotype": row.phenotype, "source_group": row.get("source_group"), "Kleborate_ST": row.get("Kleborate_ST"), "assembly_path": str(p), "mapping_error": "exact_sequence_not_found"})
    pd.DataFrame(mappings).to_csv(out / "unitig_assembly_mappings.csv", index=False); pd.DataFrame(ref_rows).to_csv(out / "unitig_reference_mappings.csv", index=False)
    if context_records: SeqIO.write(context_records, out / "replicated_unitig_contexts.fna", "fasta")
    summary = {"n_strict_candidates": int(len(strict)), "n_reference_exact_mappings": int(len(ref_rows)), "n_assembly_context_mappings": int(sum("context_id" in x for x in mappings)), "n_context_records": int(len(context_records)), "flank": a.flank, "boundary": "Context extraction and automated annotation are localization aids. Functional mechanism and novelty require manual sequence review, database/literature comparison, and biological validation."}
    (out / "UNITIG_CONTEXT_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n"); hashes = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]; (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n"); print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__": main()
