#!/usr/bin/env python3
"""Build immutable A1/A2/B nucleotide, long-ORF, and RdRP query sets.

The source contigs and previously screened RdRP segments are already tracked in
the repository.  This script extracts exact records, reconstructs the specified
reverse-strand partial ORFs, and fails closed if any sequence has changed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


CODON_TABLE = {
    codon: aa
    for aa, codons in {
        "F": ("TTT", "TTC"), "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
        "I": ("ATT", "ATC", "ATA"), "M": ("ATG",), "V": ("GTT", "GTC", "GTA", "GTG"),
        "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"), "P": ("CCT", "CCC", "CCA", "CCG"),
        "T": ("ACT", "ACC", "ACA", "ACG"), "A": ("GCT", "GCC", "GCA", "GCG"),
        "Y": ("TAT", "TAC"), "*": ("TAA", "TAG", "TGA"), "H": ("CAT", "CAC"),
        "Q": ("CAA", "CAG"), "N": ("AAT", "AAC"), "K": ("AAA", "AAG"),
        "D": ("GAT", "GAC"), "E": ("GAA", "GAG"), "C": ("TGT", "TGC"),
        "W": ("TGG",), "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
        "G": ("GGT", "GGC", "GGA", "GGG"),
    }.items()
    for codon in codons
}

CANDIDATES = {
    "PNX_Picorna_A1": {
        "contig_id": "DRR853912_10707", "rdrp_id": "DRR853912_10707_frame=-1", "reverse_offset": 0,
        "nt_length": 1988, "orf_length": 606, "rdrp_raw_length": 408, "rdrp_length": 408,
        "rdrp_trim_suffix": "", "coordinates": "171-1988",
        "nt_sha256": "62c9bf75a32cbf28d74295350bda538ce702f5702722af0c8d56f0eada56d0b3",
        "orf_nt_sha256": "33c5c19cef2e558785d33ff8ee33d2149aab3e5491fd2069092930ed52e7a664",
        "orf_sha256": "d957c7a5b2c276125caf417a373719f66ec4f81da12a139415d24e9152885bbb",
        "rdrp_raw_sha256": "b8b183715ada4e9fb64186562e8c408a64146be3ac27a9fe0d8b2ef09beaad59",
        "rdrp_sha256": "b8b183715ada4e9fb64186562e8c408a64146be3ac27a9fe0d8b2ef09beaad59",
    },
    "PNX_Picorna_A2": {
        "contig_id": "DRR853912_79526", "rdrp_id": "DRR853912_79526_frame=-1", "reverse_offset": 0,
        "nt_length": 2468, "orf_length": 606, "rdrp_raw_length": 408, "rdrp_length": 408,
        "rdrp_trim_suffix": "", "coordinates": "651-2468",
        "nt_sha256": "9ac77c66fdf0e77577a089b363fdf36ef5bc7df88745197c899baf8dd0c0b4f9",
        "orf_nt_sha256": "21e2d30493473d921ab7b2cc12e716788d7397cd51ee80aa3c0df3c42ae2b8ff",
        "orf_sha256": "079fc741207f3fdd09532b47925771f317271c6e4973dd354e169a449d1a0052",
        "rdrp_raw_sha256": "0bec0c140edf1bd9f98c21a95186f4c78f6d6c5461a7b5cc96ddded01b9fdcc9",
        "rdrp_sha256": "0bec0c140edf1bd9f98c21a95186f4c78f6d6c5461a7b5cc96ddded01b9fdcc9",
    },
    "PNX_Picorna_B": {
        "contig_id": "DRR853910_21434", "rdrp_id": "DRR853910_21434_frame=-3", "reverse_offset": 2,
        "nt_length": 3458, "orf_length": 1135, "rdrp_raw_length": 438, "rdrp_length": 435,
        "rdrp_trim_suffix": "**L", "coordinates": "52-3456",
        "nt_sha256": "d3db01132aacf3e860ef247496969f7df0c39d178ceb6a77e9966e866dc82a29",
        "orf_nt_sha256": "48e121c6f861e459b79affe412db7d319d82f00deb5fdcfe69caca87c2f6677a",
        "orf_sha256": "3f15c030914b6e83962c925f8f6d21b745191f180027a51a3b192fe1866f8199",
        "rdrp_raw_sha256": "a4b4735d9fa7bacda45271a7765346c871431f95c88911ce4e8149f12b7c2742",
        "rdrp_sha256": "4b7e33aa3e4be2e42196b01224f9fba93baeb3d1d244ce8d3f0838b82ebe00b1",
    },
}

PIPELINE_CONTROLS = {
    "PNX_Duplo_A_control": {
        "source_id": "PNX_Duplo_A", "rdrp_id": "DRR853908_24399_frame=+3", "forward_offset": 2,
        "nt_length": 1489, "orf_length": 476, "coordinates": "60-1487",
        "nt_sha256": "52fc1f2657ace79429770f3ce12867b77d5dafc8b58e236fbc5643766c50bbb8",
        "orf_sha256": "279b18559b7b426cc68513b33b81b0ab2b720640c844ebbbb0187f2141d48d81",
        "rdrp_raw_length": 364, "rdrp_raw_sha256": "cc6f8bc4d580ab9b4fad9a5e1595ddd1580375048447284c2870bd982ee28552",
        "rdrp_supported_length": 363, "rdrp_supported_sha256": "93054b91cf405ca1932d550893806ef5e78e45afb4072962d92f0ef0a03a922e",
        "rdrp_trim_suffix": "E", "rdrp_cleaning_reason": "terminal residue arose from incomplete terminal-codon padding and is not contig-supported",
    },
    "PNX_Duplo_B_control": {
        "source_id": "PNX_Duplo_B", "rdrp_id": "DRR853908_33091_frame=+1", "forward_offset": 0,
        "nt_length": 1256, "orf_length": 340, "coordinates": "235-1254",
        "nt_sha256": "ed4f061115d1e4ca341048f73ec35d4b4af88ec8487964585185f6baed855b67",
        "orf_sha256": "1b9586afecce2eb34a5de42f7cfb0bcdc030e548d589406f9fb66a95578ff4b1",
        "rdrp_raw_length": 334, "rdrp_raw_sha256": "520024a0cdf37f1aa3b8aba4358165fd1cda298277db898a53c744f1f3beea3a",
        "rdrp_supported_length": 334, "rdrp_supported_sha256": "520024a0cdf37f1aa3b8aba4358165fd1cda298277db898a53c744f1f3beea3a",
        "rdrp_trim_suffix": "", "rdrp_cleaning_reason": "none",
    },
}


def read_fasta(paths: list[Path]) -> dict[str, str]:
    records: dict[str, str] = {}
    for path in paths:
        name: str | None = None
        for raw in path.read_text().splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise SystemExit(f"duplicate FASTA identifier: {name}")
                records[name] = ""
            elif name is None:
                raise SystemExit(f"sequence before FASTA header in {path}")
            else:
                records[name] += line.upper()
    return records


def sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def translate(sequence: str, offset: int) -> str:
    return "".join(CODON_TABLE.get(sequence[i:i + 3], "X") for i in range(offset, len(sequence) - 2, 3))


def write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    with path.open("w") as handle:
        for name, description, sequence in records:
            handle.write(f">{name} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--contig-part", action="append", type=Path, required=True)
    parser.add_argument("--rdrp", type=Path, required=True)
    parser.add_argument("--control-contigs", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path, required=True)
    parser.add_argument("--strict-members", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise SystemExit(f"output directory must be empty: {args.out}")

    expected_input_hashes = {
        args.candidate_manifest: "fdbdff3d43af176f6304666a115845ff25d0f5f9dc17562afc6c643e29fa1420",
        args.strict_members: "c81481b938ba72a658408ca844f452d86d3c507911c3dae9169a668fb9f7352f",
    }
    for path, expected_hash in expected_input_hashes.items():
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed != expected_hash:
            raise SystemExit(f"lineage-definition input changed: {path} {observed} != {expected_hash}")
    parent_manifest = json.loads(args.candidate_manifest.read_text())
    with args.strict_members.open(newline="") as handle:
        strict_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_parent_members = {
        "PNX_Picorna_A": {
            "DRR853910_26168_frame=-1", "DRR853912_10707_frame=-1", "DRR853912_79526_frame=-1",
        },
        "PNX_Picorna_B": {"DRR853910_21434_frame=-3", "DRR853912_43327_frame=-1"},
    }
    strict_by_query = {row["query"]: row for row in strict_rows}
    for lineage, expected_members in expected_parent_members.items():
        observed_manifest = set(parent_manifest.get(lineage, {}).get("members", []))
        observed_strict = {row["query"] for row in strict_rows if row["lineage"] == lineage}
        if observed_manifest != expected_members or observed_strict != expected_members:
            raise SystemExit(f"parent-lineage crosswalk changed for {lineage}")

    contigs = read_fasta(args.contig_part)
    rdrps = read_fasta([args.rdrp])
    control_contigs = read_fasta([args.control_contigs])
    nt_records: list[tuple[str, str, str]] = []
    orf_nt_records: list[tuple[str, str, str]] = []
    orf_records: list[tuple[str, str, str]] = []
    rdrp_records: list[tuple[str, str, str]] = []
    manifest: list[dict[str, object]] = []
    sequences: dict[str, dict[str, str]] = {}

    for candidate, expected in CANDIDATES.items():
        nt = contigs.get(str(expected["contig_id"]), "")
        rdrp_raw = rdrps.get(str(expected["rdrp_id"]), "")
        if not nt or not rdrp_raw:
            raise SystemExit(f"missing exact source sequence for {candidate}")
        if len(rdrp_raw) != expected["rdrp_raw_length"] or sha(rdrp_raw) != expected["rdrp_raw_sha256"]:
            raise SystemExit(f"{candidate} raw screened RdRP sequence changed")
        trim_suffix = str(expected["rdrp_trim_suffix"])
        if trim_suffix:
            if not rdrp_raw.endswith(trim_suffix):
                raise SystemExit(f"{candidate} expected terminal translation artifact is absent")
            rdrp = rdrp_raw[:-len(trim_suffix)]
        else:
            rdrp = rdrp_raw
        rc = reverse_complement(nt)
        offset = int(expected["reverse_offset"])
        translated = translate(rc, offset)
        orf = translated.split("*", 1)[0]
        orf_nt = rc[offset:offset + len(orf) * 3]
        stop_codon = rc[offset + len(orf) * 3:offset + len(orf) * 3 + 3]
        checks = {
            "nt_length": len(nt), "orf_length": len(orf), "rdrp_length": len(rdrp),
            "nt_sha256": sha(nt), "orf_nt_sha256": sha(orf_nt),
            "orf_sha256": sha(orf), "rdrp_sha256": sha(rdrp),
        }
        for key, observed in checks.items():
            if observed != expected[key]:
                raise SystemExit(f"{candidate} {key} mismatch: {observed} != {expected[key]}")
        if rdrp not in orf:
            raise SystemExit(f"{candidate} screened RdRP segment is not contained in the reconstructed ORF")
        if "GDD" not in rdrp or "GDY" not in rdrp:
            raise SystemExit(f"{candidate} lacks the expected screened RdRP motif pattern")
        if stop_codon not in {"TAA", "TAG", "TGA"}:
            raise SystemExit(f"{candidate} has no immediate downstream in-frame stop after the partial ORF")
        observed_start = len(nt) - (offset + len(orf_nt)) + 1
        observed_end = len(nt) - offset
        observed_coordinates = f"{observed_start}-{observed_end}"
        if observed_coordinates != expected["coordinates"]:
            raise SystemExit(f"{candidate} reconstructed reverse-strand coordinates changed: {observed_coordinates}")
        rdrp_aa_start = orf.index(rdrp) + 1
        rdrp_aa_end = rdrp_aa_start + len(rdrp) - 1
        nt_records.append((candidate, f"source={expected['contig_id']}; associated_root_RNAseq_contig", nt))
        orf_nt_records.append((candidate, f"source={expected['contig_id']}; partial_ORF_CDS; strand=-; coordinates={observed_coordinates}", orf_nt))
        orf_records.append((candidate, f"source={expected['contig_id']}; partial_ORF; strand=-; coordinates={observed_coordinates}", orf))
        rdrp_records.append((candidate, f"source={expected['rdrp_id']}; screened_RdRP_segment", rdrp))
        manifest.append({
            "candidate": candidate, "source_contig": expected["contig_id"], "source_rdrp": expected["rdrp_id"],
            "strand": "-", "frame": int(expected["reverse_offset"]) + 1,
            "original_coordinates_1based": observed_coordinates,
            "nt_length": len(nt), "orf_aa_length": len(orf), "rdrp_aa_length": len(rdrp),
            "rdrp_within_orf_aa_coordinates_1based": f"{rdrp_aa_start}-{rdrp_aa_end}",
            "orf_nt_length": len(orf_nt), "immediate_downstream_stop_codon": stop_codon,
            "source_rdrp_raw_aa_length": len(rdrp_raw),
            "removed_post_stop_noncontiguous_suffix": trim_suffix,
            "sequence_cleaning_reason": "screening-window suffix after the reconstructed ORF stop was excluded" if trim_suffix else "none",
            "nt_sha256": sha(nt), "orf_nt_sha256": sha(orf_nt),
            "orf_sha256": sha(orf), "source_rdrp_raw_sha256": sha(rdrp_raw), "rdrp_sha256": sha(rdrp),
            "boundary_interpretation": "N-terminal contig boundary open; downstream in-frame stop present",
            "claim_boundary": "partial Picornavirales-like RNA sequence candidate; not a complete genome or formal virus species",
        })
        sequences[candidate] = {"nt": nt, "orf": orf, "rdrp": rdrp}

    control_nt_records: list[tuple[str, str, str]] = []
    control_orf_records: list[tuple[str, str, str]] = []
    control_manifest: list[dict[str, object]] = []
    for control, expected in PIPELINE_CONTROLS.items():
        nt = control_contigs.get(str(expected["source_id"]), "")
        rdrp_raw = rdrps.get(str(expected["rdrp_id"]), "")
        if len(nt) != expected["nt_length"] or sha(nt) != expected["nt_sha256"]:
            raise SystemExit(f"{control} exact nucleotide control changed")
        if len(rdrp_raw) != expected["rdrp_raw_length"] or sha(rdrp_raw) != expected["rdrp_raw_sha256"]:
            raise SystemExit(f"{control} raw screened RdRP sequence changed")
        trim_suffix = str(expected["rdrp_trim_suffix"])
        if trim_suffix:
            if not rdrp_raw.endswith(trim_suffix):
                raise SystemExit(f"{control} expected unsupported terminal residue is absent")
            rdrp_supported = rdrp_raw[:-len(trim_suffix)]
        else:
            rdrp_supported = rdrp_raw
        if len(rdrp_supported) != expected["rdrp_supported_length"] or sha(rdrp_supported) != expected["rdrp_supported_sha256"]:
            raise SystemExit(f"{control} contig-supported screened RdRP sequence changed")
        offset = int(expected["forward_offset"])
        translated = translate(nt, offset)
        segments = [(len(segment), start, segment) for start,segment in (
            (sum(len(x)+1 for x in translated.split('*')[:i]), x)
            for i,x in enumerate(translated.split('*'))
        ) if segment]
        _, aa_start, orf = max(segments)
        if len(orf) != expected["orf_length"] or sha(orf) != expected["orf_sha256"]:
            raise SystemExit(f"{control} exact positive-control ORF changed")
        if rdrp_supported not in orf:
            raise SystemExit(f"{control} contig-supported screened RdRP is not contained in its reconstructed ORF")
        observed_start = offset + aa_start * 3 + 1
        observed_end = offset + (aa_start + len(orf)) * 3
        observed_coordinates = f"{observed_start}-{observed_end}"
        if observed_coordinates != expected["coordinates"]:
            raise SystemExit(f"{control} reconstructed coordinates changed: {observed_coordinates}")
        control_nt_records.append((control, f"source={expected['source_id']}; Durnavirales_like_pipeline_positive_control", nt))
        control_orf_records.append((control, f"source={expected['source_id']}; Durnavirales_like_pipeline_positive_control_ORF; coordinates={observed_coordinates}", orf))
        control_manifest.append({
            "control":control, "source_id":expected["source_id"], "source_rdrp":expected["rdrp_id"],
            "strand":"+", "frame":offset+1, "original_coordinates_1based":observed_coordinates, "nt_length":len(nt),
            "orf_aa_length":len(orf), "nt_sha256":sha(nt), "orf_sha256":sha(orf),
            "source_rdrp_raw_aa_length":len(rdrp_raw), "contig_supported_rdrp_aa_length":len(rdrp_supported),
            "removed_unsupported_terminal_suffix":trim_suffix,
            "sequence_cleaning_reason":expected["rdrp_cleaning_reason"],
            "source_rdrp_raw_sha256":sha(rdrp_raw), "contig_supported_rdrp_sha256":sha(rdrp_supported),
            "role":"Durnavirales-like pipeline positive control; exact known-virus identity not asserted; excluded from novelty candidates",
        })

    write_fasta(args.out / "panax_three_contigs.fna", nt_records)
    write_fasta(args.out / "panax_three_partial_orfs.fna", orf_nt_records)
    write_fasta(args.out / "panax_three_partial_orfs.faa", orf_records)
    write_fasta(args.out / "panax_three_rdrp_segments.faa", rdrp_records)
    write_fasta(args.out / "panax_pipeline_controls_contigs.fna", control_nt_records)
    write_fasta(args.out / "panax_pipeline_controls_orfs.faa", control_orf_records)
    write_fasta(args.out / "panax_candidates_plus_controls_contigs.fna", nt_records + control_nt_records)
    write_fasta(args.out / "panax_candidates_plus_controls_orfs.faa", orf_records + control_orf_records)
    (args.out / "QUERY_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (args.out / "QUERY_MANIFEST.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(manifest)
    (args.out / "PIPELINE_CONTROL_MANIFEST.json").write_text(json.dumps(control_manifest, indent=2) + "\n")
    with (args.out / "PIPELINE_CONTROL_MANIFEST.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(control_manifest[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(control_manifest)

    separability = []
    names = list(CANDIDATES)
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            a, b = sequences[left]["nt"], sequences[right]["nt"]
            k = 31
            ka = {a[x:x + k] for x in range(len(a) - k + 1)}
            kb = {b[x:x + k] for x in range(len(b) - k + 1)}
            aa, bb = sequences[left]["orf"], sequences[right]["orf"]
            positional = sum(x == y for x, y in zip(aa, bb)) / max(len(aa), len(bb))
            separability.append({
                "candidate_1": left, "candidate_2": right, "k": k,
                "shared_exact_kmers": len(ka & kb), "candidate_1_unique_kmers": len(ka - kb),
                "candidate_2_unique_kmers": len(kb - ka),
                "unaligned_positional_orf_identity": f"{positional:.6f}",
                "note": "diagnostic sequence separability only; not a species-delimitation test",
            })
    with (args.out / "PAIRWISE_SEPARABILITY.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(separability[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(separability)

    # Preserve how the old coarse Picorna A/B lineage groups were split into
    # exact A1/A2/B candidates and retain the omitted cross-run members.
    crosswalk_spec = {
        "PNX_Picorna_A1": (
            "PNX_Picorna_A", "DRR853912_10707_frame=-1",
            {
                "DRR853912_10707_frame=-1": "selected_source",
                "DRR853910_26168_frame=-1": "independent_run_support",
                "DRR853912_79526_frame=-1": "separated_sibling_candidate",
            },
            "long exact A1 contig selected; the former coarse A group is split because A1 and A2 are sequence-distinguishable",
        ),
        "PNX_Picorna_A2": (
            "PNX_Picorna_A", "DRR853912_79526_frame=-1",
            {
                "DRR853912_79526_frame=-1": "selected_source",
                "DRR853910_26168_frame=-1": "related_parent_member_not_independent_A2_support",
                "DRR853912_10707_frame=-1": "separated_sibling_candidate",
            },
            "distinct long A2 contig retained separately rather than merged into coarse Picorna A",
        ),
        "PNX_Picorna_B": (
            "PNX_Picorna_B", "DRR853910_21434_frame=-3",
            {
                "DRR853910_21434_frame=-3": "selected_source",
                "DRR853912_43327_frame=-1": "independent_run_support",
            },
            "longest contig-supported B ORF selected; shorter independent-run member retained in provenance",
        ),
    }
    crosswalk: list[dict[str, object]] = []
    for new_candidate, (parent, selected_member, roles, reason) in crosswalk_spec.items():
        selected = sequences[new_candidate]["rdrp"]
        for member in parent_manifest[parent]["members"]:
            member_sequence = rdrps[member]
            if selected in member_sequence:
                aligned = len(selected); matches = len(selected)
                relation = "exact_selected_sequence_within_member"
            elif member_sequence in selected:
                aligned = len(member_sequence); matches = len(member_sequence)
                relation = "exact_member_sequence_within_selected"
            elif len(selected) == len(member_sequence):
                aligned = len(selected); matches = sum(a == b for a, b in zip(selected, member_sequence))
                relation = "ungapped_equal_length_comparison"
            else:
                raise SystemExit(f"unsupported crosswalk relation for {new_candidate}/{member}")
            crosswalk.append({
                "new_candidate": new_candidate, "parent_lineage": parent,
                "selected_member": selected_member, "parent_member": member,
                "member_role": roles[member], "member_run": strict_by_query[member]["accession"],
                "selected_rdrp_length": len(selected), "member_rdrp_length": len(member_sequence),
                "relation": relation, "aligned_aa": aligned, "identical_aa": matches,
                "identity_over_aligned_pct": f"{100 * matches / aligned:.6f}",
                "selected_rdrp_sha256": sha(selected), "member_raw_rdrp_sha256": sha(member_sequence),
                "selection_reason": reason,
            })
    with (args.out / "LINEAGE_CROSSWALK.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(crosswalk[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(crosswalk)
    (args.out / "LINEAGE_CROSSWALK.json").write_text(json.dumps(crosswalk, indent=2) + "\n")

    provenance = {
        "source_contig_files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.contig_part},
        "source_rdrp_file": {args.rdrp.name: hashlib.sha256(args.rdrp.read_bytes()).hexdigest()},
        "source_control_contig_file": {args.control_contigs.name: hashlib.sha256(args.control_contigs.read_bytes()).hexdigest()},
        "source_candidate_manifest": {args.candidate_manifest.name: hashlib.sha256(args.candidate_manifest.read_bytes()).hexdigest()},
        "source_strict_members": {args.strict_members.name: hashlib.sha256(args.strict_members.read_bytes()).hexdigest()},
        "candidate_count": len(manifest),
        "pipeline_control_count":len(control_manifest),
        "interpretation_boundary": "These immutable queries support audits of partial sequence candidates only.",
    }
    (args.out / "QUERY_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    lines = []
    for path in sorted(args.out.glob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            lines.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (args.out / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
