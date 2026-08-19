#!/usr/bin/env python3
"""Build an auditable targeted colistin-pathway variant matrix from K. pneumoniae assemblies.

This script is deliberately discovery-oriented but conservative. It aligns a pinned panel of
reference proteins (plus selected promoter+gene nucleotide regions) to each assembly, calls
reference-relative sequence features, records assembly quality, and emits an Rtab matrix for
population-structure-aware association with pyseer.

It does not label any sequence change as causal or novel. Reference-relative differences can be
lineage markers, neutral polymorphisms, assembly artefacts, or established mechanisms.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

TARGETS: dict[str, tuple[str, ...]] = {
    "mgrB": ("mgrb", "negative regulator of pho q", "negative regulator of phoq"),
    "phoP": ("phop",), "phoQ": ("phoq",), "pmrA": ("pmra", "basr"),
    "pmrB": ("pmrb", "bass"), "pmrD": ("pmrd",), "crrA": ("crra",),
    "crrB": ("crrb",), "crrC": ("crrc",),
    "eptA": ("epta", "pmrc", "phosphoethanolamine transferase epta"),
    "ugd": ("ugd", "pmre"), "arnA": ("arna", "pmrh"),
    "arnB": ("arnb", "pmri"), "arnC": ("arnc", "pmrf"),
    "arnD": ("arnd",), "arnT": ("arnt", "pmrk"), "lpxA": ("lpxa",),
    "lpxC": ("lpxc",), "lpxD": ("lpxd",), "lpxM": ("lpxm", "msbb"),
    "pagP": ("pagp",), "ramA": ("rama",), "acrB": ("acrb",),
}
PROMOTER_UPSTREAM = {"mgrB": 250, "pmrD": 200, "eptA": 200}
ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")


@dataclass
class RefGene:
    name: str
    protein: str
    nucleotide: str
    promoter_gene: str | None
    upstream: int
    record_id: str
    locus_tag: str
    product: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--labels", required=True)
    p.add_argument("--assemblies-dir", required=True)
    p.add_argument("--reference-gbff", required=True)
    p.add_argument("--out", default="targeted_variant_scan")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-qcov", type=float, default=0.80)
    p.add_argument("--min-pident", type=float, default=70.0)
    p.add_argument("--min-n50", type=int, default=5000)
    p.add_argument("--max-contigs", type=int, default=1000)
    p.add_argument("--min-genome", type=int, default=4_500_000)
    p.add_argument("--max-genome", type=int, default=7_500_000)
    p.add_argument("--max-n-fraction", type=float, default=0.02)
    return p.parse_args()


def norm(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", x.lower())


def feature_strings(feat: Any) -> list[str]:
    vals: list[str] = []
    for key in ("gene", "gene_synonym", "locus_tag", "product", "note", "standard_name"):
        vals.extend(str(v) for v in feat.qualifiers.get(key, []))
    return vals


def match_target(strings: Iterable[str]) -> str | None:
    raw = [str(x).lower() for x in strings]
    compact = [norm(x) for x in raw]
    for target, aliases in TARGETS.items():
        for alias in (target, *aliases):
            a = norm(alias)
            if any(c == a for c in compact):
                return target
    joined = " | ".join(raw)
    for target, aliases in TARGETS.items():
        for alias in aliases:
            if len(alias) >= 8 and alias.lower() in joined:
                return target
    return None


def extract_reference(gbff: Path, out: Path) -> dict[str, RefGene]:
    found: dict[str, list[tuple[int, Any, Any]]] = defaultdict(list)
    records = list(SeqIO.parse(str(gbff), "genbank"))
    if not records:
        raise RuntimeError(f"No GenBank records parsed from {gbff}")
    for rec in records:
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            target = match_target(feature_strings(feat))
            if not target:
                continue
            trans = feat.qualifiers.get("translation", [])
            if not trans:
                continue
            score = 0
            genes = [norm(x) for x in feat.qualifiers.get("gene", [])]
            synonyms = [norm(x) for x in feat.qualifiers.get("gene_synonym", [])]
            if norm(target) in genes:
                score += 100
            if norm(target) in synonyms:
                score += 80
            if "plasmid" not in rec.description.lower():
                score += 10
            score += min(len(trans[0]), 1000) // 100
            found[target].append((score, rec, feat))

    refs: dict[str, RefGene] = {}
    for target in TARGETS:
        if target not in found:
            continue
        _, rec, feat = sorted(found[target], key=lambda x: x[0], reverse=True)[0]
        protein = str(feat.qualifiers["translation"][0]).replace("*", "")
        nt = str(feat.extract(rec.seq)).upper()
        upstream = PROMOTER_UPSTREAM.get(target, 0)
        pg: str | None = None
        if upstream:
            start = int(feat.location.start); end = int(feat.location.end)
            strand = int(feat.location.strand or 1)
            if strand == 1:
                promoter = rec.seq[max(0, start - upstream):start]
                if len(promoter) == upstream:
                    pg = str(promoter + feat.extract(rec.seq)).upper()
            else:
                promoter = rec.seq[end:min(len(rec.seq), end + upstream)].reverse_complement()
                if len(promoter) == upstream:
                    pg = str(promoter + feat.extract(rec.seq)).upper()
        refs[target] = RefGene(
            target, protein, nt, pg, upstream, rec.id,
            ";".join(feat.qualifiers.get("locus_tag", ["-"])),
            ";".join(feat.qualifiers.get("product", ["-"])),
        )

    out.mkdir(parents=True, exist_ok=True)
    SeqIO.write([SeqRecord(Seq(v.protein), id=k, description=f"{v.record_id} {v.locus_tag} {v.product}") for k, v in refs.items()], out / "reference_targets.faa", "fasta")
    SeqIO.write([SeqRecord(Seq(v.nucleotide), id=k, description=f"{v.record_id} {v.locus_tag} {v.product}") for k, v in refs.items()], out / "reference_targets.fna", "fasta")
    promoter_records = [SeqRecord(Seq(v.promoter_gene), id=f"{k}__PROMOTER_GENE", description=f"upstream={v.upstream}") for k, v in refs.items() if v.promoter_gene]
    if promoter_records:
        SeqIO.write(promoter_records, out / "reference_promoter_gene.fna", "fasta")
    pd.DataFrame([{"gene": k, "protein_length": len(v.protein), "nt_length": len(v.nucleotide), "promoter_upstream": v.upstream, "record": v.record_id, "locus_tag": v.locus_tag, "product": v.product} for k, v in refs.items()]).to_csv(out / "reference_target_manifest.csv", index=False)
    missing = sorted(set(TARGETS) - set(refs))
    (out / "REFERENCE_TARGETS_MISSING.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    return refs


def fasta_stats(path: Path) -> dict[str, Any]:
    lengths: list[int] = []; n_bases = 0; total = 0
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq).upper(); lengths.append(len(seq)); total += len(seq); n_bases += seq.count("N")
    lengths.sort(reverse=True); half = total / 2; cum = 0; n50 = 0
    for x in lengths:
        cum += x
        if cum >= half:
            n50 = x; break
    return {"genome_length": total, "contigs": len(lengths), "n50": n50, "n_fraction": (n_bases / total) if total else 1.0}


def find_assemblies(labels: pd.DataFrame, assembly_dir: Path) -> tuple[dict[str, Path], list[str]]:
    paths = list(assembly_dir.rglob("*.fna")); by_base: dict[str, list[Path]] = defaultdict(list)
    for p in paths:
        m = ACC_RE.search(str(p))
        if m:
            by_base[m.group(0).split(".")[0]].append(p)
    mapping: dict[str, Path] = {}; missing: list[str] = []
    for acc in labels["assembly_ID"].astype(str):
        exact = [p for p in paths if acc in str(p)]
        if exact:
            mapping[acc] = sorted(exact, key=lambda x: len(str(x)))[0]
        elif by_base.get(acc.split(".")[0]):
            mapping[acc] = sorted(by_base[acc.split(".")[0]], key=lambda x: len(str(x)))[0]
        else:
            missing.append(acc)
    return mapping, missing


FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qseq", "sseq"]


def run_one(item: tuple[str, Path], protein_query: Path, promoter_query: Path | None, raw_dir: Path) -> tuple[str, str | None]:
    acc, fasta = item
    try:
        subprocess.run(["tblastn", "-query", str(protein_query), "-subject", str(fasta), "-seg", "no", "-evalue", "1e-6", "-max_target_seqs", "10", "-max_hsps", "10", "-outfmt", "6 " + " ".join(FIELDS), "-out", str(raw_dir / f"{acc}.tblastn.tsv")], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if promoter_query and promoter_query.exists():
            subprocess.run(["blastn", "-task", "blastn", "-query", str(promoter_query), "-subject", str(fasta), "-dust", "no", "-evalue", "1e-20", "-max_target_seqs", "10", "-max_hsps", "10", "-outfmt", "6 " + " ".join(FIELDS), "-out", str(raw_dir / f"{acc}.blastn.tsv")], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return acc, None
    except subprocess.CalledProcessError as e:
        return acc, (e.stderr or str(e))[-2000:]


def read_blast(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=FIELDS)
    x = pd.read_csv(path, sep="\t", names=FIELDS, dtype={"qseqid": str, "sseqid": str, "qseq": str, "sseq": str})
    for c in ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    return x


def best_hsp(x: pd.DataFrame) -> pd.Series | None:
    return None if x.empty else x.sort_values(["bitscore", "pident", "length"], ascending=False).iloc[0]


def call_aa_features(gene: str, row: pd.Series) -> tuple[set[str], list[dict[str, Any]]]:
    features: set[str] = set(); details: list[dict[str, Any]] = []; qpos = int(row.qstart) - 1
    insertion: list[str] = []; insertion_anchor = qpos
    for q, s in zip(str(row.qseq), str(row.sseq)):
        if q == "-":
            if s != "-": insertion.append(s); insertion_anchor = qpos
            continue
        if insertion:
            seq = "".join(insertion); fid = f"{gene}:INS_AFTER_{insertion_anchor}:{seq}"; features.add(fid)
            details.append({"feature": fid, "gene": gene, "type": "aa_insertion", "position": insertion_anchor, "ref": "-", "alt": seq}); insertion = []
        qpos += 1
        if s == "-": fid, typ = f"{gene}:DEL_{qpos}{q}", "aa_deletion"
        elif s == "*": fid, typ = f"{gene}:STOP_{qpos}{q}", "premature_stop"
        elif q != s and q != "X" and s != "X": fid, typ = f"{gene}:{q}{qpos}{s}", "aa_substitution"
        else: continue
        features.add(fid); details.append({"feature": fid, "gene": gene, "type": typ, "position": qpos, "ref": q, "alt": s})
    if insertion:
        seq = "".join(insertion); fid = f"{gene}:INS_AFTER_{insertion_anchor}:{seq}"; features.add(fid)
        details.append({"feature": fid, "gene": gene, "type": "aa_insertion", "position": insertion_anchor, "ref": "-", "alt": seq})
    return features, details


def call_promoter_features(gene: str, upstream: int, row: pd.Series) -> tuple[set[str], list[dict[str, Any]]]:
    features: set[str] = set(); details: list[dict[str, Any]] = []; qpos = int(row.qstart) - 1
    for q, s in zip(str(row.qseq), str(row.sseq)):
        if q == "-": continue
        qpos += 1
        if qpos > upstream: continue
        rel = qpos - upstream - 1
        if s == "-": fid, typ = f"{gene}_PROM:{rel}{q}>DEL", "promoter_deletion"
        elif q != s and s != "N" and q != "N": fid, typ = f"{gene}_PROM:{rel}{q}>{s}", "promoter_substitution"
        else: continue
        features.add(fid); details.append({"feature": fid, "gene": gene, "type": typ, "position": rel, "ref": q, "alt": s})
    return features, details


def main() -> None:
    a = parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    raw_dir = out / "raw_alignments"; raw_dir.mkdir(exist_ok=True)
    labels = pd.read_csv(a.labels, dtype={"assembly_ID": str})
    if not {"assembly_ID", "phenotype"}.issubset(labels.columns): raise SystemExit("labels must contain assembly_ID and phenotype")
    labels = labels.drop_duplicates("assembly_ID", keep="first"); labels = labels[labels.phenotype.astype(str).isin(["R", "S"])]
    refs = extract_reference(Path(a.reference_gbff), out / "reference")
    if len(refs) < 10: raise RuntimeError(f"Only {len(refs)} target genes found in reference")
    mapping, missing = find_assemblies(labels, Path(a.assemblies_dir)); (out / "MISSING_ASSEMBLIES.txt").write_text("\n".join(missing) + ("\n" if missing else ""))

    qc_rows = []; valid: dict[str, Path] = {}
    for acc, pth in mapping.items():
        st = fasta_stats(pth)
        qc_pass = a.min_genome <= st["genome_length"] <= a.max_genome and st["contigs"] <= a.max_contigs and st["n50"] >= a.min_n50 and st["n_fraction"] <= a.max_n_fraction
        qc_rows.append({"assembly_ID": acc, "path": str(pth), **st, "qc_pass": qc_pass})
        if qc_pass: valid[acc] = pth
    qc = pd.DataFrame(qc_rows); qc.to_csv(out / "assembly_qc.csv", index=False)
    if len(valid) < 100: raise RuntimeError(f"Only {len(valid)} assemblies passed QC")

    protein_query = out / "reference" / "reference_targets.faa"; promoter_query = out / "reference" / "reference_promoter_gene.fna"
    errors = []
    with cf.ThreadPoolExecutor(max_workers=max(1, a.threads)) as ex:
        for acc, err in (f.result() for f in cf.as_completed([ex.submit(run_one, item, protein_query, promoter_query if promoter_query.exists() else None, raw_dir) for item in valid.items()])):
            if err: errors.append({"assembly_ID": acc, "error": err})
    pd.DataFrame(errors).to_csv(out / "alignment_errors.csv", index=False)

    sample_features: dict[str, set[str]] = defaultdict(set); call_rows = []; feature_meta: dict[str, dict[str, Any]] = {}
    for acc in sorted(valid):
        tb = read_blast(raw_dir / f"{acc}.tblastn.tsv")
        for gene, ref in refs.items():
            hit = best_hsp(tb[tb.qseqid == gene])
            if hit is None:
                fid = f"{gene}:ABSENT_OR_UNDETECTED"; sample_features[acc].add(fid)
                feature_meta.setdefault(fid, {"feature": fid, "gene": gene, "type": "gene_absent_or_undetected", "position": None, "ref": None, "alt": None})
                call_rows.append({"assembly_ID": acc, "gene": gene, "status": "absent", "qcov": 0.0}); continue
            qcov = (float(hit.qend) - float(hit.qstart) + 1.0) / float(hit.qlen)
            status = "full" if qcov >= a.min_qcov and float(hit.pident) >= a.min_pident else "partial_or_divergent"
            if status != "full":
                fid = f"{gene}:PARTIAL_OR_DIVERGENT"; sample_features[acc].add(fid)
                feature_meta.setdefault(fid, {"feature": fid, "gene": gene, "type": "partial_or_divergent", "position": None, "ref": None, "alt": None})
            feats, det = call_aa_features(gene, hit); sample_features[acc].update(feats)
            for d in det: feature_meta.setdefault(d["feature"], d)
            if feats:
                nf = f"{gene}:ANY_NONREFERENCE_AA"; sample_features[acc].add(nf)
                feature_meta.setdefault(nf, {"feature": nf, "gene": gene, "type": "gene_burden_any_nonreference", "position": None, "ref": None, "alt": None})
            lof = any(feature_meta[f]["type"] in {"premature_stop", "aa_deletion", "aa_insertion"} for f in feats)
            if lof or status != "full":
                lf = f"{gene}:POTENTIAL_LOF"; sample_features[acc].add(lf)
                feature_meta.setdefault(lf, {"feature": lf, "gene": gene, "type": "potential_loss_of_function", "position": None, "ref": None, "alt": None})
            call_rows.append({"assembly_ID": acc, "gene": gene, "status": status, "qcov": qcov, "pident": float(hit.pident), "bitscore": float(hit.bitscore), "subject": str(hit.sseqid), "sstart": int(hit.sstart), "send": int(hit.send), "n_features": len(feats)})

        bn = read_blast(raw_dir / f"{acc}.blastn.tsv")
        for gene, ref in refs.items():
            if not ref.promoter_gene: continue
            hit = best_hsp(bn[bn.qseqid == f"{gene}__PROMOTER_GENE"])
            if hit is None:
                fid = f"{gene}_PROM:UNDETECTED"; sample_features[acc].add(fid)
                feature_meta.setdefault(fid, {"feature": fid, "gene": gene, "type": "promoter_undetected", "position": None, "ref": None, "alt": None}); continue
            promoter_cov_end = min(int(hit.qend), ref.upstream); promoter_cov_start = min(max(int(hit.qstart), 1), ref.upstream)
            pcov = max(0, promoter_cov_end - promoter_cov_start + 1) / ref.upstream
            if pcov < 0.8:
                fid = f"{gene}_PROM:PARTIAL_OR_STRUCTURAL"; sample_features[acc].add(fid)
                feature_meta.setdefault(fid, {"feature": fid, "gene": gene, "type": "promoter_partial_or_structural", "position": None, "ref": None, "alt": None})
            feats, det = call_promoter_features(gene, ref.upstream, hit); sample_features[acc].update(feats)
            for d in det: feature_meta.setdefault(d["feature"], d)

    pd.DataFrame(call_rows).to_csv(out / "target_gene_calls.csv", index=False)
    samples = sorted(valid); features = sorted({f for fs in sample_features.values() for f in fs})
    with open(out / "targeted_variants.Rtab", "w") as fh:
        fh.write("Gene\t" + "\t".join(samples) + "\n")
        for f in features: fh.write(f + "\t" + "\t".join("1" if f in sample_features[s] else "0" for s in samples) + "\n")
    meta = pd.DataFrame([feature_meta[f] for f in features]); counts = {f: sum(f in sample_features[s] for s in samples) for f in features}
    meta["n_present"] = meta.feature.map(counts); meta["frequency"] = meta.n_present / len(samples)
    meta["known_mechanism_screen"] = meta.apply(lambda r: bool(r["gene"] in {"mgrB", "pmrB"} and r["type"] in {"potential_loss_of_function", "premature_stop", "aa_deletion", "aa_insertion", "gene_absent_or_undetected", "partial_or_divergent", "promoter_partial_or_structural"}), axis=1)
    meta.to_csv(out / "targeted_variant_metadata.csv", index=False)
    sample_meta = labels[labels.assembly_ID.isin(samples)].merge(qc, on="assembly_ID", how="left", validate="one_to_one"); sample_meta.to_csv(out / "targeted_sample_metadata.csv", index=False)
    with open(out / "phenotypes.tsv", "w") as fh:
        fh.write("samples\tphenotype\n"); ymap = sample_meta.set_index("assembly_ID").phenotype.astype(str).to_dict()
        for s in samples: fh.write(f"{s}\t{1 if ymap[s] == 'R' else 0}\n")
    summary = {"n_labels": int(len(labels)), "n_assemblies_found": int(len(mapping)), "n_missing": int(len(missing)), "n_qc_pass": int(len(samples)), "phenotypes_qc_pass": sample_meta.phenotype.value_counts().to_dict(), "n_reference_targets": int(len(refs)), "reference_targets": sorted(refs), "n_features": int(len(features)), "n_alignment_errors": int(len(errors)), "boundary": "All calls are reference-relative genomic features. They are not causal or novel determinants. Association, source-held-out replication, known-mechanism/literature review, and biological validation remain required."}
    (out / "TARGETED_VARIANT_MATRIX_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    hashes = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n"); print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
