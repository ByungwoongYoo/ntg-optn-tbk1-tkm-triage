#!/usr/bin/env python3
"""Corpus-scale audit of disease specificity in herbal network-pharmacology target lists.

The script is designed to run from a clean GitHub Actions runner. It harvests a frozen
Europe PMC Open Access corpus, downloads full-text XML, extracts explicit author-reported
hub/core/key target lists, maps title diseases to a predeclared lexicon, and runs
label-permutation, frequency-preserving null, and disease-recoverability analyses.

It never treats the extracted lists as proof of efficacy or causal mechanism.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import gzip
import hashlib
import html
import io
import json
import math
import os
import random
import re
import shutil
import statistics
import sys
import tarfile
import textwrap
import time
import traceback
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np
import pandas as pd
from lxml import etree
from scipy import sparse
from scipy.stats import entropy
from sklearn.dummy import DummyClassifier
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Frozen design constants
# ---------------------------------------------------------------------------
SNAPSHOT_DATE = "2026-08-16"
BASE_SEED = 20260816
EUROPE_PMC_QUERY = (
    'OPEN_ACCESS:Y AND FIRST_PDATE:[2015-01-01 TO 2026-08-16] '
    'AND TITLE:"network pharmacology" '
    'AND (decoction OR formula OR herb OR herbal OR phytochemical '
    'OR "traditional Chinese medicine" OR "traditional medicine" '
    'OR natural OR plant)'
)
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
NCBI_GENE_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
USER_AGENT = (
    "HerbalDiseaseSpecificityAudit/1.0 "
    "(independent reproducibility research; contact yoonge3@gmail.com)"
)

# Fixed exclusions before outcome inspection.
TITLE_EXCLUDE_PATTERNS = [
    r"\bsystematic review\b",
    r"\bmeta-analysis\b",
    r"\bmeta analysis\b",
    r"\bbibliometric\b",
    r"\bscoping review\b",
    r"\bnarrative review\b",
    r"\breview\b",
    r"\beditorial\b",
    r"\bprotocol\b",
    r"\boverview\b",
    r"\bperspective\b",
    r"\bcommentary\b",
    r"\bretraction\b",
]
INTERVENTION_PATTERNS = [
    r"\bdecoction\b", r"\bformula\b", r"\bherbal\b", r"\bherb\b",
    r"\btraditional chinese medicine\b", r"\btcm\b", r"\btraditional medicine\b",
    r"\bphytochemical\b", r"\bmedicinal plant\b", r"\bnatural product\b",
    r"\bextract\b", r"\bcompound\b", r"\bflavonoid\b", r"\balkaloid\b",
    r"\bpolyphenol\b", r"\bsaponin\b", r"\bterpenoid\b", r"\bpolysaccharide\b",
    r"\bginseng\b", r"\bberberine\b", r"\bcurcumin\b", r"\bquercetin\b",
    r"\bresveratrol\b", r"\bastragalus\b", r"\bsalvia\b", r"\brhubarb\b",
]

ANCHOR_PATTERNS: list[tuple[str, str, int]] = [
    ("hub", r"\bhub\s+(?:genes?|targets?|proteins?|nodes?)\b", 8),
    ("core", r"\bcore\s+(?:genes?|targets?|proteins?|nodes?)\b", 7),
    ("key", r"\bkey\s+(?:genes?|targets?|proteins?|nodes?)\b", 6),
    ("central", r"\bcentral\s+(?:genes?|targets?|proteins?|nodes?)\b", 5),
    ("crucial", r"\bcrucial\s+(?:genes?|targets?|proteins?|nodes?)\b", 5),
    ("important", r"\bimportant\s+(?:genes?|targets?|proteins?|nodes?)\b", 4),
    ("top", r"\btop\s+\d{0,2}\s*(?:genes?|targets?|proteins?|nodes?)\b", 5),
    ("intersection", r"\bintersection\s+(?:genes?|targets?|proteins?)\b", 4),
    ("overlap", r"\boverlapping\s+(?:genes?|targets?|proteins?)\b", 4),
    ("candidate", r"\bcandidate\s+(?:genes?|targets?|proteins?)\b", 3),
    ("potential", r"\bpotential\s+(?:genes?|targets?|proteins?)\b", 2),
]
ANCHOR_REGEX = [(name, re.compile(pattern, re.I), weight) for name, pattern, weight in ANCHOR_PATTERNS]

GENE_STOPWORDS = {
    "AND", "OR", "THE", "WITH", "FROM", "THIS", "THAT", "WERE", "WAS", "ARE",
    "GENE", "GENES", "TARGET", "TARGETS", "PROTEIN", "PROTEINS", "HUB", "CORE",
    "KEY", "TOP", "PPI", "GO", "KEGG", "STRING", "TCM", "RNA", "DNA", "MRNA",
    "LNCRNA", "MIRNA", "ROS", "NO", "NOS", "MAPK", "AKT", "ERK", "JNK", "PI3K",
    "NF", "NFKB", "NF-KB", "HIF", "VEGF", "TGF", "IL", "TNF-A", "COX", "ACE",
    "AUC", "ROC", "CI", "SD", "SEM", "PCA", "GSEA", "GEO", "CTD", "OMIM",
    "DRUGBANK", "GENECARDS", "DISGENET", "CYP", "ADME", "OB", "DL", "FDR",
}

# Predeclared title lexicon. Specific phrases must precede broad phrases.
# label, broad category, regex patterns
DISEASE_LEXICON: list[tuple[str, str, list[str]]] = [
    ("breast_cancer", "cancer", [r"\bbreast cancer\b", r"\bbreast carcinoma\b"]),
    ("lung_cancer", "cancer", [r"\blung cancer\b", r"\bnon[- ]small cell lung cancer\b", r"\bnsclc\b"]),
    ("hepatocellular_carcinoma", "cancer", [r"\bhepatocellular carcinoma\b", r"\bliver cancer\b", r"\bhcc\b"]),
    ("colorectal_cancer", "cancer", [r"\bcolorectal cancer\b", r"\bcolon cancer\b", r"\brectal cancer\b"]),
    ("gastric_cancer", "cancer", [r"\bgastric cancer\b", r"\bstomach cancer\b", r"\bgastric carcinoma\b"]),
    ("prostate_cancer", "cancer", [r"\bprostate cancer\b", r"\bprostatic carcinoma\b"]),
    ("pancreatic_cancer", "cancer", [r"\bpancreatic cancer\b", r"\bpancreatic carcinoma\b"]),
    ("ovarian_cancer", "cancer", [r"\bovarian cancer\b", r"\bovarian carcinoma\b"]),
    ("cervical_cancer", "cancer", [r"\bcervical cancer\b", r"\bcervical carcinoma\b"]),
    ("endometrial_cancer", "cancer", [r"\bendometrial cancer\b", r"\bendometrial carcinoma\b"]),
    ("esophageal_cancer", "cancer", [r"\besophageal cancer\b", r"\boesophageal cancer\b"]),
    ("nasopharyngeal_carcinoma", "cancer", [r"\bnasopharyngeal carcinoma\b"]),
    ("oral_cancer", "cancer", [r"\boral cancer\b", r"\boral squamous cell carcinoma\b"]),
    ("thyroid_cancer", "cancer", [r"\bthyroid cancer\b", r"\bthyroid carcinoma\b"]),
    ("bladder_cancer", "cancer", [r"\bbladder cancer\b", r"\burothelial carcinoma\b"]),
    ("renal_cell_carcinoma", "cancer", [r"\brenal cell carcinoma\b", r"\bkidney cancer\b"]),
    ("glioma", "cancer", [r"\bglioblastoma\b", r"\bglioma\b"]),
    ("leukemia", "cancer", [r"\bleukemia\b", r"\bleukaemia\b"]),
    ("multiple_myeloma", "cancer", [r"\bmultiple myeloma\b"]),
    ("lymphoma", "cancer", [r"\blymphoma\b"]),
    ("cancer_unspecified", "cancer", [r"\bcancer\b", r"\bcarcinoma\b", r"\btumou?r\b"]),

    ("rheumatoid_arthritis", "autoimmune", [r"\brheumatoid arthritis\b"]),
    ("osteoarthritis", "musculoskeletal", [r"\bosteoarthritis\b", r"\bknee osteoarthritis\b"]),
    ("gout", "musculoskeletal", [r"\bgout(?:y arthritis)?\b", r"\bhyperuricemia\b"]),
    ("osteoporosis", "musculoskeletal", [r"\bosteoporosis\b"]),
    ("intervertebral_disc_degeneration", "musculoskeletal", [r"\bintervertebral disc degeneration\b", r"\bdisc degeneration\b"]),
    ("ankylosing_spondylitis", "autoimmune", [r"\bankylosing spondylitis\b"]),
    ("fibromyalgia", "musculoskeletal", [r"\bfibromyalgia\b"]),

    ("ulcerative_colitis", "gastrointestinal", [r"\bulcerative colitis\b"]),
    ("crohn_disease", "gastrointestinal", [r"\bcrohn'?s disease\b"]),
    ("irritable_bowel_syndrome", "gastrointestinal", [r"\birritable bowel syndrome\b", r"\bibs\b"]),
    ("gastritis", "gastrointestinal", [r"\bchronic atrophic gastritis\b", r"\bgastritis\b"]),
    ("gastric_ulcer", "gastrointestinal", [r"\bgastric ulcer\b", r"\bpeptic ulcer\b"]),
    ("constipation", "gastrointestinal", [r"\bconstipation\b"]),
    ("acute_pancreatitis", "gastrointestinal", [r"\bacute pancreatitis\b", r"\bpancreatitis\b"]),
    ("inflammatory_bowel_disease", "gastrointestinal", [r"\binflammatory bowel disease\b"]),

    ("diabetic_nephropathy", "renal_metabolic", [r"\bdiabetic nephropathy\b", r"\bdiabetic kidney disease\b"]),
    ("type_2_diabetes", "metabolic", [r"\btype 2 diabetes(?: mellitus)?\b", r"\bt2dm\b"]),
    ("diabetes", "metabolic", [r"\bdiabetes mellitus\b", r"\bdiabetes\b"]),
    ("obesity", "metabolic", [r"\bobesity\b"]),
    ("metabolic_syndrome", "metabolic", [r"\bmetabolic syndrome\b"]),
    ("nonalcoholic_fatty_liver", "hepatic_metabolic", [r"\bnon[- ]alcoholic fatty liver disease\b", r"\bnafld\b", r"\bmetabolic dysfunction-associated steatotic liver disease\b", r"\bmasld\b"]),
    ("liver_fibrosis", "hepatic", [r"\bliver fibrosis\b", r"\bhepatic fibrosis\b"]),
    ("acute_liver_injury", "hepatic", [r"\bacute liver injury\b", r"\bliver injury\b"]),
    ("chronic_kidney_disease", "renal", [r"\bchronic kidney disease\b", r"\brenal fibrosis\b"]),

    ("alzheimer_disease", "neurological", [r"\balzheimer'?s disease\b", r"\balzheimer disease\b"]),
    ("parkinson_disease", "neurological", [r"\bparkinson'?s disease\b", r"\bparkinson disease\b"]),
    ("ischemic_stroke", "neurological", [r"\bischemic stroke\b", r"\bischaemic stroke\b", r"\bcerebral ischemia\b", r"\bcerebral ischaemia\b"]),
    ("intracerebral_hemorrhage", "neurological", [r"\bintracerebral hemorrhage\b", r"\bcerebral hemorrhage\b"]),
    ("epilepsy", "neurological", [r"\bepilepsy\b"]),
    ("neuropathic_pain", "neurological", [r"\bneuropathic pain\b"]),
    ("spinal_cord_injury", "neurological", [r"\bspinal cord injury\b"]),
    ("vascular_dementia", "neurological", [r"\bvascular dementia\b"]),
    ("depression", "psychiatric", [r"\bmajor depressive disorder\b", r"\bdepression\b", r"\bantidepressant\b"]),
    ("anxiety", "psychiatric", [r"\banxiety\b"]),
    ("insomnia", "psychiatric", [r"\binsomnia\b"]),

    ("atherosclerosis", "cardiovascular", [r"\batherosclerosis\b"]),
    ("coronary_heart_disease", "cardiovascular", [r"\bcoronary heart disease\b", r"\bcoronary artery disease\b", r"\bangina pectoris\b"]),
    ("myocardial_infarction", "cardiovascular", [r"\bmyocardial infarction\b", r"\bmyocardial ischemia\b", r"\bmyocardial ischaemia\b"]),
    ("heart_failure", "cardiovascular", [r"\bheart failure\b"]),
    ("hypertension", "cardiovascular", [r"\bhypertension\b"]),
    ("cardiac_fibrosis", "cardiovascular", [r"\bcardiac fibrosis\b", r"\bmyocardial fibrosis\b"]),
    ("arrhythmia", "cardiovascular", [r"\barrhythmia\b", r"\batrial fibrillation\b"]),

    ("asthma", "respiratory", [r"\basthma\b"]),
    ("copd", "respiratory", [r"\bchronic obstructive pulmonary disease\b", r"\bcopd\b"]),
    ("pulmonary_fibrosis", "respiratory", [r"\bpulmonary fibrosis\b"]),
    ("acute_lung_injury", "respiratory", [r"\bacute lung injury\b", r"\bards\b", r"\bacute respiratory distress syndrome\b"]),
    ("pneumonia", "respiratory_infectious", [r"\bpneumonia\b"]),

    ("sepsis", "infectious_inflammatory", [r"\bsepsis\b"]),
    ("covid_19", "infectious", [r"\bcovid[- ]?19\b", r"\bsars[- ]cov[- ]2\b"]),
    ("influenza", "infectious", [r"\binfluenza\b"]),

    ("psoriasis", "dermatologic_autoimmune", [r"\bpsoriasis\b"]),
    ("atopic_dermatitis", "dermatologic", [r"\batopic dermatitis\b", r"\beczema\b"]),
    ("acne", "dermatologic", [r"\bacne\b"]),
    ("alopecia_areata", "dermatologic_autoimmune", [r"\balopecia areata\b"]),
    ("vitiligo", "dermatologic_autoimmune", [r"\bvitiligo\b"]),

    ("endometriosis", "gynecologic", [r"\bendometriosis\b"]),
    ("polycystic_ovary_syndrome", "gynecologic_metabolic", [r"\bpolycystic ovar(?:y|ian) syndrome\b", r"\bpcos\b"]),
    ("preeclampsia", "gynecologic", [r"\bpreeclampsia\b", r"\bpre-eclampsia\b"]),
    ("uterine_fibroids", "gynecologic", [r"\buterine fibroids?\b", r"\buterine leiomyoma\b"]),
    ("male_infertility", "reproductive", [r"\bmale infertility\b"]),
    ("female_infertility", "reproductive", [r"\bfemale infertility\b", r"\binfertility\b"]),

    ("periodontitis", "oral", [r"\bperiodontitis\b"]),
    ("dry_eye_disease", "ophthalmic", [r"\bdry eye disease\b"]),
    ("age_related_macular_degeneration", "ophthalmic", [r"\bage-related macular degeneration\b", r"\bmacular degeneration\b"]),
    ("glaucoma", "ophthalmic", [r"\bglaucoma\b"]),
]
DISEASE_COMPILED = [
    (label, broad, [re.compile(p, re.I) for p in patterns])
    for label, broad, patterns in DISEASE_LEXICON
]

SECTION_WEIGHTS = {
    "abstract": 7,
    "results": 6,
    "result": 6,
    "conclusion": 5,
    "discussion": 3,
    "methods": 0,
    "materials and methods": 0,
    "introduction": -1,
    "table": 5,
    "figure": 4,
    "unknown": 1,
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def stable_int(text: str, seed: int = BASE_SEED) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}|{text}".encode()).digest()[:8], "big")


def norm_space(text: str) -> str:
    text = html.unescape(text or "")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def get_url(url: str, *, timeout: int = 90, retries: int = 5) -> bytes:
    last: Exception | None = None
    headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(min(15, 1.5 ** attempt + random.random()))
    raise RuntimeError(f"GET failed after {retries} attempts: {url}: {last}")


def get_json(url: str, **kwargs: Any) -> dict[str, Any]:
    return json.loads(get_url(url, **kwargs).decode("utf-8"))


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        keys: list[str] = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def title_exclusion_reason(title: str, pub_types: Sequence[str] | None = None) -> str | None:
    low = norm_space(title).lower()
    for pattern in TITLE_EXCLUDE_PATTERNS:
        if re.search(pattern, low, re.I):
            return f"title_pattern:{pattern}"
    pts = " ".join(pub_types or []).lower()
    if any(x in pts for x in ["review", "meta-analysis", "editorial", "letter", "comment"]):
        return f"publication_type:{pts}"
    return None


def has_intervention_signal(text: str) -> bool:
    return any(re.search(p, text, re.I) for p in INTERVENTION_PATTERNS)


def disease_from_title(title: str) -> tuple[str | None, str | None, str | None]:
    title_clean = norm_space(title)
    matches: list[tuple[int, int, str, str, str]] = []
    for order, (label, broad, patterns) in enumerate(DISEASE_COMPILED):
        for pattern in patterns:
            m = pattern.search(title_clean)
            if m:
                matches.append((len(m.group(0)), -order, label, broad, m.group(0)))
    if not matches:
        return None, None, None
    matches.sort(reverse=True)
    _, _, label, broad, raw = matches[0]
    return label, broad, raw


def gini(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or np.all(arr == 0):
        return float("nan")
    arr = np.sort(arr)
    n = arr.size
    return float((2 * np.sum((np.arange(1, n + 1)) * arr) / (n * np.sum(arr))) - (n + 1) / n)


def percentile_ci(values: Sequence[float], alpha: float = 0.05) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    return float(np.quantile(arr, alpha / 2)), float(np.quantile(arr, 1 - alpha / 2))


# ---------------------------------------------------------------------------
# Corpus harvest
# ---------------------------------------------------------------------------
def harvest_metadata(out_dir: Path, query: str = EUROPE_PMC_QUERY) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cursor = "*"
    page = 0
    hit_count: int | None = None
    while True:
        params = {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": 1000,
            "cursorMark": cursor,
        }
        url = EPMC_SEARCH + "?" + urllib.parse.urlencode(params)
        payload = get_json(url, timeout=120)
        if hit_count is None:
            hit_count = int(payload.get("hitCount", 0))
        result = payload.get("resultList", {}).get("result", [])
        if not result:
            break
        rows.extend(result)
        page += 1
        next_cursor = payload.get("nextCursorMark")
        print(f"metadata page={page} rows={len(rows)}/{hit_count}", flush=True)
        if not next_cursor or next_cursor == cursor or len(rows) >= (hit_count or 0):
            cursor = next_cursor or cursor
            break
        cursor = next_cursor
        time.sleep(0.15)

    # Deterministic deduplication by PMCID, then DOI/PMID fallback.
    dedup: dict[str, dict[str, Any]] = {}
    for row in rows:
        key = row.get("pmcid") or row.get("doi") or row.get("pmid") or f"anon:{len(dedup)}"
        dedup.setdefault(str(key), row)
    rows = list(dedup.values())
    manifest = {
        "query": query,
        "snapshot_date": SNAPSHOT_DATE,
        "harvested_at_utc": utc_now(),
        "reported_hit_count": hit_count,
        "records_after_deduplication": len(rows),
        "last_cursor_mark": cursor,
        "endpoint": EPMC_SEARCH,
    }
    write_json(out_dir / "raw" / "metadata.json", rows)
    write_json(out_dir / "raw" / "metadata_manifest.json", manifest)
    return rows, manifest


def select_metadata_for_xml(
    metadata: Sequence[dict[str, Any]],
    *,
    max_xml: int | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    selected: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for row in metadata:
        title = row.get("title") or ""
        pub_types = (row.get("pubTypeList") or {}).get("pubType", []) if isinstance(row.get("pubTypeList"), dict) else []
        reason = title_exclusion_reason(title, pub_types)
        if reason:
            excluded.append({"pmcid": row.get("pmcid"), "title": title, "stage": "metadata", "reason": reason})
            continue
        if not row.get("pmcid"):
            excluded.append({"pmcid": None, "title": title, "stage": "metadata", "reason": "no_pmcid"})
            continue
        title_abs = f"{title} {row.get('abstractText') or ''}"
        if not has_intervention_signal(title_abs):
            excluded.append({"pmcid": row.get("pmcid"), "title": title, "stage": "metadata", "reason": "no_intervention_signal"})
            continue
        selected.append(row)

    # A pilot uses a deterministic hash sample across the whole eligible metadata frame,
    # not the API relevance order. Full analysis sets max_xml=None.
    selected.sort(key=lambda r: stable_int(str(r.get("pmcid"))))
    if max_xml is not None:
        selected = selected[:max_xml]
    return selected, excluded


def download_xml_one(row: Mapping[str, Any], xml_dir: Path) -> dict[str, Any]:
    pmcid = str(row["pmcid"])
    path = xml_dir / f"{pmcid}.xml"
    if path.exists() and path.stat().st_size > 100:
        data = path.read_bytes()
        return {"pmcid": pmcid, "status": "cached", "bytes": len(data), "sha256": sha256_bytes(data)}
    url = EPMC_FULLTEXT.format(pmcid=urllib.parse.quote(pmcid))
    try:
        data = get_url(url, timeout=120, retries=5)
        if b"<article" not in data[:5000] and b"<book-part" not in data[:5000]:
            raise RuntimeError("response_does_not_look_like_article_xml")
        path.write_bytes(data)
        return {"pmcid": pmcid, "status": "downloaded", "bytes": len(data), "sha256": sha256_bytes(data)}
    except Exception as exc:  # noqa: BLE001
        return {"pmcid": pmcid, "status": "failed", "error": repr(exc)}


def download_xmls(rows: Sequence[dict[str, Any]], out_dir: Path, workers: int = 8) -> list[dict[str, Any]]:
    xml_dir = out_dir / "raw" / "xml"
    xml_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_xml_one, row, xml_dir): row for row in rows}
        for idx, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            if idx % 50 == 0 or idx == len(futures):
                ok = sum(r.get("status") != "failed" for r in results)
                print(f"fulltext {idx}/{len(futures)} ok={ok}", flush=True)
    results.sort(key=lambda r: r["pmcid"])
    write_json(out_dir / "raw" / "xml_manifest.json", results)
    return results


# ---------------------------------------------------------------------------
# Gene dictionary
# ---------------------------------------------------------------------------
@dataclass
class GeneDictionary:
    official_symbols: set[str]
    unique_alias_to_symbol: dict[str, str]
    ambiguous_aliases: set[str]
    source_sha256: str
    source_url: str


def normalize_gene_token(token: str) -> str:
    token = token.strip().upper().replace("−", "-").replace("–", "-")
    token = re.sub(r"[^A-Z0-9-]", "", token)
    return token


def build_gene_dictionary(out_dir: Path) -> GeneDictionary:
    raw_path = out_dir / "raw" / "Homo_sapiens.gene_info.gz"
    if not raw_path.exists():
        raw_path.write_bytes(get_url(NCBI_GENE_INFO, timeout=180))
    data = raw_path.read_bytes()
    official: set[str] = set()
    alias_candidates: defaultdict[str, set[str]] = defaultdict(set)
    with gzip.open(io.BytesIO(data), "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("#tax_id") != "9606":
                continue
            symbol = normalize_gene_token(row.get("Symbol") or "")
            if not symbol or symbol == "-":
                continue
            official.add(symbol)
            alias_candidates[symbol].add(symbol)
            for field in [row.get("Synonyms") or "", row.get("Symbol_from_nomenclature_authority") or ""]:
                for alias in field.split("|"):
                    alias = normalize_gene_token(alias)
                    if len(alias) >= 2 and alias not in GENE_STOPWORDS and alias != "-":
                        alias_candidates[alias].add(symbol)
    unique_alias: dict[str, str] = {}
    ambiguous: set[str] = set()
    for alias, symbols in alias_candidates.items():
        if len(symbols) == 1:
            unique_alias[alias] = next(iter(symbols))
        else:
            ambiguous.add(alias)
    manifest = {
        "source_url": NCBI_GENE_INFO,
        "source_sha256": sha256_bytes(data),
        "official_symbol_count": len(official),
        "unique_alias_count": len(unique_alias),
        "ambiguous_alias_count": len(ambiguous),
    }
    write_json(out_dir / "raw" / "gene_dictionary_manifest.json", manifest)
    return GeneDictionary(official, unique_alias, ambiguous, manifest["source_sha256"], NCBI_GENE_INFO)


# ---------------------------------------------------------------------------
# XML parsing and target extraction
# ---------------------------------------------------------------------------
@dataclass
class TextChunk:
    section: str
    kind: str
    order: int
    text: str


@dataclass
class ExtractionCandidate:
    anchor: str
    anchor_weight: int
    section: str
    kind: str
    order: int
    span: str
    genes: list[str]
    raw_tokens: list[str]
    explicit_count: int | None
    score: float
    density: float
    count_match: bool
    confidence: str


@dataclass
class ArticleExtraction:
    pmcid: str
    pmid: str | None
    doi: str | None
    title: str
    year: int | None
    journal: str | None
    author_string: str | None
    article_type: str | None
    disease_label: str | None
    broad_disease: str | None
    disease_raw: str | None
    genes: list[str]
    extraction_confidence: str
    extraction_score: float
    anchor: str
    source_section: str
    source_kind: str
    source_order: int
    source_span: str
    explicit_count: int | None
    count_match: bool
    candidate_count: int
    xml_sha256: str


def element_text(element: etree._Element) -> str:
    return norm_space(" ".join(element.itertext()))


def section_name(element: etree._Element) -> str:
    # Walk to nearest sec and use its title; otherwise infer abstract/table/figure.
    node: etree._Element | None = element
    while node is not None:
        tag = etree.QName(node).localname.lower() if isinstance(node.tag, str) else ""
        if tag == "abstract":
            return "abstract"
        if tag == "table-wrap":
            return "table"
        if tag == "fig":
            return "figure"
        if tag == "sec":
            title = node.find("./title")
            if title is not None:
                return element_text(title).lower() or "unknown"
            return "unknown"
        node = node.getparent()
    return "unknown"


def parse_article_xml(path: Path) -> tuple[dict[str, Any], list[TextChunk]]:
    parser = etree.XMLParser(recover=True, huge_tree=True, resolve_entities=False)
    root = etree.fromstring(path.read_bytes(), parser=parser)
    article = root if etree.QName(root).localname == "article" else root.find(".//article")
    if article is None:
        article = root
    article_type = article.get("article-type")
    title_el = article.find(".//article-title")
    title = element_text(title_el) if title_el is not None else ""
    pmcid = None
    pmid = None
    doi = None
    for aid in article.findall(".//article-id"):
        typ = aid.get("pub-id-type")
        value = element_text(aid)
        if typ == "pmc":
            pmcid = value if value.upper().startswith("PMC") else f"PMC{value}"
        elif typ == "pmid":
            pmid = value
        elif typ == "doi":
            doi = value
    journal_el = article.find(".//journal-title")
    journal = element_text(journal_el) if journal_el is not None else None
    year = None
    for y in article.findall(".//pub-date/year"):
        txt = element_text(y)
        if txt.isdigit():
            year = int(txt)
            break
    authors = []
    for contrib in article.findall(".//contrib[@contrib-type='author']"):
        surname = contrib.find(".//surname")
        given = contrib.find(".//given-names")
        name = " ".join(x for x in [element_text(given) if given is not None else "", element_text(surname) if surname is not None else ""] if x)
        if name:
            authors.append(name)
    metadata = {
        "pmcid": pmcid,
        "pmid": pmid,
        "doi": doi,
        "title": title,
        "journal": journal,
        "year": year,
        "author_string": "; ".join(authors) if authors else None,
        "article_type": article_type,
    }

    chunks: list[TextChunk] = []
    order = 0
    # Paragraphs, table cells/captions, and figure captions; exclude reference list.
    xpath = (
        ".//abstract//p | .//body//p | .//table-wrap//caption | .//table-wrap//td | "
        ".//table-wrap//th | .//fig//caption"
    )
    for el in article.xpath(xpath):
        if el.xpath("ancestor::ref-list"):
            continue
        text = element_text(el)
        if len(text) < 20:
            continue
        tag = etree.QName(el).localname.lower() if isinstance(el.tag, str) else "unknown"
        kind = "paragraph"
        if tag in {"td", "th"}:
            kind = "table_cell"
        elif tag == "caption":
            kind = "caption"
        section = section_name(el)
        chunks.append(TextChunk(section=section, kind=kind, order=order, text=text))
        order += 1
    return metadata, chunks


def section_weight(section: str, kind: str) -> int:
    sec = section.lower()
    if sec in SECTION_WEIGHTS:
        base = SECTION_WEIGHTS[sec]
    else:
        base = 1
        for key, value in SECTION_WEIGHTS.items():
            if key != "unknown" and key in sec:
                base = value
                break
    if kind == "table_cell":
        base += 1
    return base


def explicit_count_near(text: str, anchor_start: int, anchor_end: int) -> int | None:
    local = text[max(0, anchor_start - 80): min(len(text), anchor_end + 80)]
    patterns = [
        r"\b(?:top\s+)?(\d{1,2})\s+(?:hub|core|key|central|crucial|important|candidate|potential|intersection|overlapping)\b",
        r"\b(\d{1,2})\s+(?:genes?|targets?|proteins?)\s+(?:were|was|are|is)\s+(?:identified|selected|screened|obtained|considered)\b",
        r"\b(?:identified|selected|screened|obtained)\s+(\d{1,2})\s+(?:hub|core|key|central|crucial|important|candidate|potential)?\s*(?:genes?|targets?|proteins?)\b",
    ]
    for pattern in patterns:
        m = re.search(pattern, local, re.I)
        if m:
            value = int(m.group(1))
            if 2 <= value <= 60:
                return value
    return None


def span_around_anchor(text: str, start: int, end: int) -> str:
    # Prefer sentence-level extraction, but include an adjacent sentence after a colon/list lead-in.
    left_candidates = [text.rfind(x, 0, start) for x in [". ", "? ", "! ", "; "]]
    left = max(left_candidates)
    left = 0 if left < 0 else left + 2
    right_positions = [p for x in [". ", "? ", "! "] if (p := text.find(x, end)) >= 0]
    right = min(right_positions) + 1 if right_positions else min(len(text), end + 700)
    span = text[left:right]
    # If the anchor sentence contains very few symbols, append one following sentence.
    if len(span) < 180 and right < len(text):
        next_positions = [p for x in [". ", "? ", "! "] if (p := text.find(x, right + 1)) >= 0]
        next_right = min(next_positions) + 1 if next_positions else min(len(text), right + 500)
        span = text[left:next_right]
    # Cap to prevent whole-paragraph extraction.
    if len(span) > 1200:
        anchor_rel = start - left
        a = max(0, anchor_rel - 250)
        b = min(len(span), anchor_rel + 850)
        span = span[a:b]
    return norm_space(span)


def gene_tokens_from_text(text: str, gene_dict: GeneDictionary) -> tuple[list[str], list[str]]:
    # Gene-like tokens, including hyphenated aliases. We require either all caps or at least one digit
    # to avoid mapping ordinary title-case words through aliases.
    raw = re.findall(r"(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9]{1,14}(?:-[A-Za-z0-9]{1,8})?)(?![A-Za-z0-9])", text)
    genes: list[str] = []
    raw_kept: list[str] = []
    seen: set[str] = set()
    for token in raw:
        norm = normalize_gene_token(token)
        if not norm or norm in GENE_STOPWORDS:
            continue
        is_gene_style = token.isupper() or any(ch.isdigit() for ch in token) or (token[:1].isupper() and token[1:].islower() is False)
        if not is_gene_style:
            continue
        candidates = [norm, norm.replace("-", "")]
        symbol = None
        for candidate in candidates:
            if candidate in gene_dict.official_symbols:
                symbol = candidate
                break
            if candidate in gene_dict.unique_alias_to_symbol:
                symbol = gene_dict.unique_alias_to_symbol[candidate]
                break
        if symbol and symbol not in seen:
            seen.add(symbol)
            genes.append(symbol)
            raw_kept.append(token)
    return genes, raw_kept


def candidate_confidence(
    *,
    genes: Sequence[str],
    density: float,
    explicit_count: int | None,
    count_match: bool,
    section: str,
    anchor: str,
) -> str:
    sec = section.lower()
    if count_match and 3 <= len(genes) <= 40:
        return "high"
    if anchor in {"hub", "core", "key", "top"} and density >= 0.10 and any(k in sec for k in ["abstract", "result", "conclusion", "table"]):
        return "high"
    if anchor in {"hub", "core", "key", "central", "crucial", "top"} and density >= 0.06 and 3 <= len(genes) <= 40:
        return "medium"
    return "low"


def extract_candidates(chunks: Sequence[TextChunk], gene_dict: GeneDictionary) -> list[ExtractionCandidate]:
    candidates: list[ExtractionCandidate] = []
    for chunk in chunks:
        for anchor, regex, anchor_weight in ANCHOR_REGEX:
            for match in regex.finditer(chunk.text):
                span = span_around_anchor(chunk.text, match.start(), match.end())
                genes, raw_tokens = gene_tokens_from_text(span, gene_dict)
                if not (3 <= len(genes) <= 60):
                    continue
                word_count = max(1, len(re.findall(r"\b\w+\b", span)))
                density = len(genes) / word_count
                explicit_count = explicit_count_near(chunk.text, match.start(), match.end())
                count_match = explicit_count is not None and abs(explicit_count - len(genes)) <= 1
                score = float(anchor_weight + section_weight(chunk.section, chunk.kind))
                score += min(4.0, density * 20)
                if count_match:
                    score += 5
                elif explicit_count is not None:
                    score -= min(4, abs(explicit_count - len(genes)) / 2)
                if ":" in span or ";" in span:
                    score += 1
                if 5 <= len(genes) <= 30:
                    score += 1
                if re.search(r"\b(?:for example|such as|including but not limited to)\b", span, re.I):
                    score -= 1
                confidence = candidate_confidence(
                    genes=genes,
                    density=density,
                    explicit_count=explicit_count,
                    count_match=count_match,
                    section=chunk.section,
                    anchor=anchor,
                )
                candidates.append(
                    ExtractionCandidate(
                        anchor=anchor,
                        anchor_weight=anchor_weight,
                        section=chunk.section,
                        kind=chunk.kind,
                        order=chunk.order,
                        span=span,
                        genes=genes,
                        raw_tokens=raw_tokens,
                        explicit_count=explicit_count,
                        score=score,
                        density=density,
                        count_match=count_match,
                        confidence=confidence,
                    )
                )
    # Deduplicate identical sets; retain best-scoring source.
    best_by_set: dict[tuple[str, ...], ExtractionCandidate] = {}
    for candidate in candidates:
        key = tuple(sorted(candidate.genes))
        old = best_by_set.get(key)
        if old is None or candidate.score > old.score:
            best_by_set[key] = candidate
    return sorted(best_by_set.values(), key=lambda c: (c.score, c.confidence == "high", -c.order), reverse=True)


def choose_primary_candidate(candidates: Sequence[ExtractionCandidate]) -> ExtractionCandidate | None:
    if not candidates:
        return None
    # High confidence is primary even if a medium candidate scores slightly higher.
    high = [c for c in candidates if c.confidence == "high"]
    if high:
        return max(high, key=lambda c: c.score)
    medium = [c for c in candidates if c.confidence == "medium"]
    if medium:
        return max(medium, key=lambda c: c.score)
    return None


def article_extraction_from_xml(
    row: Mapping[str, Any],
    xml_path: Path,
    gene_dict: GeneDictionary,
) -> tuple[ArticleExtraction | None, dict[str, Any] | None, list[dict[str, Any]]]:
    try:
        parsed_meta, chunks = parse_article_xml(xml_path)
    except Exception as exc:  # noqa: BLE001
        return None, {"pmcid": row.get("pmcid"), "title": row.get("title"), "stage": "xml_parse", "reason": repr(exc)}, []

    title = parsed_meta.get("title") or row.get("title") or ""
    article_type = (parsed_meta.get("article_type") or "").lower()
    if article_type and any(x in article_type for x in ["review", "editorial", "letter", "commentary", "correction"]):
        return None, {"pmcid": row.get("pmcid"), "title": title, "stage": "xml", "reason": f"article_type:{article_type}"}, []
    reason = title_exclusion_reason(title, (row.get("pubTypeList") or {}).get("pubType", []) if isinstance(row.get("pubTypeList"), dict) else [])
    if reason:
        return None, {"pmcid": row.get("pmcid"), "title": title, "stage": "xml", "reason": reason}, []
    if not has_intervention_signal(f"{title} {row.get('abstractText') or ''}"):
        return None, {"pmcid": row.get("pmcid"), "title": title, "stage": "xml", "reason": "no_intervention_signal"}, []

    candidates = extract_candidates(chunks, gene_dict)
    candidate_rows = []
    for rank, c in enumerate(candidates[:10], start=1):
        candidate_rows.append({
            "pmcid": row.get("pmcid"),
            "rank": rank,
            "anchor": c.anchor,
            "section": c.section,
            "kind": c.kind,
            "order": c.order,
            "gene_count": len(c.genes),
            "genes": "|".join(c.genes),
            "explicit_count": c.explicit_count,
            "count_match": c.count_match,
            "density": c.density,
            "score": c.score,
            "confidence": c.confidence,
            "source_span": c.span,
        })
    primary = choose_primary_candidate(candidates)
    if primary is None:
        return None, {"pmcid": row.get("pmcid"), "title": title, "stage": "extraction", "reason": "no_high_or_medium_explicit_target_list"}, candidate_rows

    disease_label, broad_disease, disease_raw = disease_from_title(title)
    extraction = ArticleExtraction(
        pmcid=str(row.get("pmcid") or parsed_meta.get("pmcid")),
        pmid=str(row.get("pmid") or parsed_meta.get("pmid")) if (row.get("pmid") or parsed_meta.get("pmid")) else None,
        doi=str(row.get("doi") or parsed_meta.get("doi")) if (row.get("doi") or parsed_meta.get("doi")) else None,
        title=title,
        year=int(row.get("pubYear") or parsed_meta.get("year")) if str(row.get("pubYear") or parsed_meta.get("year") or "").isdigit() else None,
        journal=row.get("journalTitle") or parsed_meta.get("journal"),
        author_string=row.get("authorString") or parsed_meta.get("author_string"),
        article_type=parsed_meta.get("article_type"),
        disease_label=disease_label,
        broad_disease=broad_disease,
        disease_raw=disease_raw,
        genes=primary.genes,
        extraction_confidence=primary.confidence,
        extraction_score=primary.score,
        anchor=primary.anchor,
        source_section=primary.section,
        source_kind=primary.kind,
        source_order=primary.order,
        source_span=primary.span,
        explicit_count=primary.explicit_count,
        count_match=primary.count_match,
        candidate_count=len(candidates),
        xml_sha256=sha256_bytes(xml_path.read_bytes()),
    )
    return extraction, None, candidate_rows


def extract_all(
    rows: Sequence[dict[str, Any]],
    out_dir: Path,
    gene_dict: GeneDictionary,
) -> tuple[list[ArticleExtraction], list[dict[str, Any]], list[dict[str, Any]]]:
    xml_dir = out_dir / "raw" / "xml"
    extractions: list[ArticleExtraction] = []
    exclusions: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []
    for idx, row in enumerate(rows, start=1):
        xml_path = xml_dir / f"{row['pmcid']}.xml"
        if not xml_path.exists():
            exclusions.append({"pmcid": row.get("pmcid"), "title": row.get("title"), "stage": "download", "reason": "xml_missing"})
            continue
        extraction, exclusion, candidate_rows = article_extraction_from_xml(row, xml_path, gene_dict)
        candidates.extend(candidate_rows)
        if extraction is not None:
            extractions.append(extraction)
        if exclusion is not None:
            exclusions.append(exclusion)
        if idx % 100 == 0 or idx == len(rows):
            print(f"extract {idx}/{len(rows)} included={len(extractions)}", flush=True)
    return extractions, exclusions, candidates


# ---------------------------------------------------------------------------
# Statistical analysis
# ---------------------------------------------------------------------------
def build_incidence(articles: Sequence[ArticleExtraction], genes_override: Sequence[Sequence[str]] | None = None) -> tuple[sparse.csr_matrix, list[str]]:
    genesets = [list(gs) for gs in genes_override] if genes_override is not None else [a.genes for a in articles]
    vocabulary = sorted({g for genes in genesets for g in genes})
    index = {g: j for j, g in enumerate(vocabulary)}
    rows: list[int] = []
    cols: list[int] = []
    for i, genes in enumerate(genesets):
        for g in set(genes):
            if g in index:
                rows.append(i)
                cols.append(index[g])
    data = np.ones(len(rows), dtype=np.uint8)
    matrix = sparse.csr_matrix((data, (rows, cols)), shape=(len(genesets), len(vocabulary)), dtype=np.uint8)
    return matrix, vocabulary


def jaccard_matrix(matrix: sparse.csr_matrix) -> np.ndarray:
    matrix = matrix.astype(np.float64)
    intersections = (matrix @ matrix.T).toarray()
    sizes = np.asarray(matrix.sum(axis=1)).ravel()
    unions = sizes[:, None] + sizes[None, :] - intersections
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.divide(intersections, unions, out=np.zeros_like(intersections), where=unions > 0)
    np.fill_diagonal(jac, 1.0)
    return jac


def pairwise_delta(sim: np.ndarray, labels: np.ndarray) -> tuple[float, float, float, int, int]:
    iu = np.triu_indices(len(labels), k=1)
    same = labels[iu[0]] == labels[iu[1]]
    values = sim[iu]
    if same.sum() == 0 or (~same).sum() == 0:
        return float("nan"), float("nan"), float("nan"), int(same.sum()), int((~same).sum())
    within = float(values[same].mean())
    between = float(values[~same].mean())
    return within, between, within - between, int(same.sum()), int((~same).sum())


def stratified_label_permutation(
    labels: np.ndarray,
    list_sizes: np.ndarray,
    years: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    # Permute within list-size quartiles and publication-era bins. Singleton strata are pooled.
    size_q = pd.qcut(pd.Series(list_sizes), q=min(4, len(np.unique(list_sizes))), labels=False, duplicates="drop").to_numpy()
    year_bins = np.digitize(years, bins=[2019, 2021, 2023, 2025])
    strata = np.array([f"{a}|{b}" for a, b in zip(size_q, year_bins)], dtype=object)
    perm = labels.copy()
    for stratum in np.unique(strata):
        idx = np.flatnonzero(strata == stratum)
        if len(idx) >= 2:
            perm[idx] = rng.permutation(perm[idx])
    return perm


def permutation_test_delta(
    sim: np.ndarray,
    labels: np.ndarray,
    *,
    n_perm: int,
    seed: int,
    stratified: bool = False,
    list_sizes: np.ndarray | None = None,
    years: np.ndarray | None = None,
) -> dict[str, Any]:
    within, between, observed, n_within, n_between = pairwise_delta(sim, labels)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        if stratified:
            assert list_sizes is not None and years is not None
            perm = stratified_label_permutation(labels, list_sizes, years, rng)
        else:
            perm = rng.permutation(labels)
        null[k] = pairwise_delta(sim, perm)[2]
    p_greater = float((1 + np.count_nonzero(null >= observed)) / (n_perm + 1))
    p_two = float((1 + np.count_nonzero(np.abs(null - null.mean()) >= abs(observed - null.mean()))) / (n_perm + 1))
    lo, hi = percentile_ci(null)
    return {
        "within_mean_jaccard": within,
        "between_mean_jaccard": between,
        "delta": observed,
        "n_within_pairs": n_within,
        "n_between_pairs": n_between,
        "n_permutations": n_perm,
        "p_one_sided_greater": p_greater,
        "p_two_sided": p_two,
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
        "null_95_interval": [lo, hi],
        "stratified": stratified,
    }


def bootstrap_delta(
    genesets: Sequence[Sequence[str]],
    labels: Sequence[str],
    *,
    n_boot: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    labels_arr = np.asarray(labels, dtype=object)
    by_class = {label: np.flatnonzero(labels_arr == label) for label in np.unique(labels_arr)}
    values: list[float] = []
    for _ in range(n_boot):
        idx: list[int] = []
        for label in sorted(by_class):
            group = by_class[label]
            idx.extend(rng.choice(group, size=len(group), replace=True).tolist())
        sampled_sets = [genesets[i] for i in idx]
        sampled_labels = labels_arr[idx]
        mat, _ = build_incidence([], genes_override=sampled_sets)
        sim = jaccard_matrix(mat)
        delta = pairwise_delta(sim, sampled_labels)[2]
        if math.isfinite(delta):
            values.append(delta)
    return percentile_ci(values)


def filter_articles_for_label(
    articles: Sequence[ArticleExtraction],
    label_attr: str,
    min_class_n: int,
    confidence: set[str],
) -> list[ArticleExtraction]:
    eligible = [a for a in articles if getattr(a, label_attr) and a.extraction_confidence in confidence]
    counts = Counter(getattr(a, label_attr) for a in eligible)
    return [a for a in eligible if counts[getattr(a, label_attr)] >= min_class_n]


def remove_top_genes(articles: Sequence[ArticleExtraction], n: int) -> tuple[list[list[str]], list[str]]:
    counter = Counter(g for a in articles for g in set(a.genes))
    removed = [g for g, _ in counter.most_common(n)]
    removed_set = set(removed)
    genesets = [[g for g in a.genes if g not in removed_set] for a in articles]
    return genesets, removed


def double_edge_swap_gene_sets(
    genesets: Sequence[Sequence[str]],
    *,
    rng: np.random.Generator,
    n_swaps: int,
) -> list[list[str]]:
    rows = [set(gs) for gs in genesets]
    nonempty = [i for i, s in enumerate(rows) if s]
    if len(nonempty) < 2:
        return [sorted(s) for s in rows]
    accepted = 0
    attempts = 0
    max_attempts = max(n_swaps * 20, 1000)
    while accepted < n_swaps and attempts < max_attempts:
        attempts += 1
        i, j = rng.choice(nonempty, size=2, replace=False)
        only_i = tuple(rows[i] - rows[j])
        only_j = tuple(rows[j] - rows[i])
        if not only_i or not only_j:
            continue
        gi = only_i[int(rng.integers(len(only_i)))]
        gj = only_j[int(rng.integers(len(only_j)))]
        rows[i].remove(gi); rows[j].remove(gj)
        rows[i].add(gj); rows[j].add(gi)
        accepted += 1
    return [sorted(s) for s in rows]


def frequency_preserving_null(
    articles: Sequence[ArticleExtraction],
    label_attr: str,
    *,
    n_random: int,
    swaps_per_edge: int,
    seed: int,
) -> dict[str, Any]:
    labels = np.asarray([getattr(a, label_attr) for a in articles], dtype=object)
    genesets = [a.genes for a in articles]
    base_mat, _ = build_incidence(articles)
    base_sim = jaccard_matrix(base_mat)
    observed = pairwise_delta(base_sim, labels)[2]
    edges = sum(len(set(gs)) for gs in genesets)
    rng = np.random.default_rng(seed)
    null = []
    for k in range(n_random):
        randomized = double_edge_swap_gene_sets(
            genesets,
            rng=rng,
            n_swaps=max(1, edges * swaps_per_edge),
        )
        mat, _ = build_incidence([], genes_override=randomized)
        sim = jaccard_matrix(mat)
        null.append(pairwise_delta(sim, labels)[2])
    arr = np.asarray(null)
    lo, hi = percentile_ci(arr)
    return {
        "observed_delta": observed,
        "n_random_matrices": n_random,
        "swaps_per_edge": swaps_per_edge,
        "null_mean": float(arr.mean()),
        "null_sd": float(arr.std(ddof=1)),
        "null_95_interval": [lo, hi],
        "p_one_sided_greater": float((1 + np.count_nonzero(arr >= observed)) / (len(arr) + 1)),
        "p_two_sided": float((1 + np.count_nonzero(np.abs(arr - arr.mean()) >= abs(observed - arr.mean()))) / (len(arr) + 1)),
    }


def disease_recoverability(
    articles: Sequence[ArticleExtraction],
    label_attr: str,
    *,
    repeats: int,
    splits: int,
    n_label_permutations: int,
    seed: int,
    remove_genes: set[str] | None = None,
) -> dict[str, Any]:
    labels = np.asarray([getattr(a, label_attr) for a in articles], dtype=object)
    genesets = [[g for g in a.genes if not remove_genes or g not in remove_genes] for a in articles]
    X, vocabulary = build_incidence([], genes_override=genesets)
    # TF-IDF reduces dominance of universally recycled genes without using outcomes.
    model = make_pipeline(TfidfTransformer(), LinearSVC(class_weight="balanced", C=1.0, random_state=seed))
    cv = RepeatedStratifiedKFold(n_splits=splits, n_repeats=repeats, random_state=seed)
    observed_rows = []
    for fold, (train, test) in enumerate(cv.split(X, labels), start=1):
        model.fit(X[train], labels[train])
        pred = model.predict(X[test])
        observed_rows.append({
            "fold": fold,
            "accuracy": accuracy_score(labels[test], pred),
            "balanced_accuracy": balanced_accuracy_score(labels[test], pred),
            "macro_f1": f1_score(labels[test], pred, average="macro", zero_division=0),
        })
    observed = {metric: float(np.mean([r[metric] for r in observed_rows])) for metric in ["accuracy", "balanced_accuracy", "macro_f1"]}
    observed_sd = {metric: float(np.std([r[metric] for r in observed_rows], ddof=1)) for metric in ["accuracy", "balanced_accuracy", "macro_f1"]}

    rng = np.random.default_rng(seed + 991)
    # Use a single stratified split set for each permuted label vector to keep runtime bounded.
    perm_metrics = []
    for p in range(n_label_permutations):
        y_perm = rng.permutation(labels)
        perm_cv = RepeatedStratifiedKFold(n_splits=splits, n_repeats=1, random_state=seed + p + 1)
        scores = []
        for train, test in perm_cv.split(X, y_perm):
            model.fit(X[train], y_perm[train])
            pred = model.predict(X[test])
            scores.append(balanced_accuracy_score(y_perm[test], pred))
        perm_metrics.append(float(np.mean(scores)))
    perm_arr = np.asarray(perm_metrics)

    class_counts = Counter(labels.tolist())
    majority = max(class_counts.values()) / len(labels)
    chance_balanced = 1 / len(class_counts)
    return {
        "n_articles": len(articles),
        "n_classes": len(class_counts),
        "class_counts": dict(sorted(class_counts.items(), key=lambda kv: (-kv[1], kv[0]))),
        "n_features": len(vocabulary),
        "cv_splits": splits,
        "cv_repeats": repeats,
        "fold_count": len(observed_rows),
        "observed_mean": observed,
        "observed_sd": observed_sd,
        "majority_accuracy": majority,
        "balanced_accuracy_chance": chance_balanced,
        "label_permutations": n_label_permutations,
        "permuted_balanced_accuracy_mean": float(perm_arr.mean()),
        "permuted_balanced_accuracy_sd": float(perm_arr.std(ddof=1)),
        "permutation_p_balanced_accuracy": float((1 + np.count_nonzero(perm_arr >= observed["balanced_accuracy"])) / (len(perm_arr) + 1)),
        "top_genes_removed": sorted(remove_genes) if remove_genes else [],
        "fold_metrics": observed_rows,
    }


def temporal_holdout(
    articles: Sequence[ArticleExtraction],
    label_attr: str,
    *,
    cutoff_year: int,
    seed: int,
) -> dict[str, Any]:
    labels = np.asarray([getattr(a, label_attr) for a in articles], dtype=object)
    years = np.asarray([a.year or 0 for a in articles], dtype=int)
    train_mask = years <= cutoff_year
    test_mask = years > cutoff_year
    train_classes = Counter(labels[train_mask].tolist())
    test_classes = Counter(labels[test_mask].tolist())
    keep_classes = {c for c in train_classes if train_classes[c] >= 3 and test_classes.get(c, 0) >= 2}
    keep = np.array([label in keep_classes for label in labels])
    train = np.flatnonzero(keep & train_mask)
    test = np.flatnonzero(keep & test_mask)
    if len(keep_classes) < 2 or len(train) < 10 or len(test) < 5:
        return {"status": "insufficient", "cutoff_year": cutoff_year, "n_classes": len(keep_classes), "n_train": len(train), "n_test": len(test)}
    X, vocabulary = build_incidence(articles)
    model = make_pipeline(TfidfTransformer(), LinearSVC(class_weight="balanced", C=1.0, random_state=seed))
    model.fit(X[train], labels[train])
    pred = model.predict(X[test])
    return {
        "status": "completed",
        "cutoff_year": cutoff_year,
        "n_classes": len(keep_classes),
        "classes": sorted(keep_classes),
        "n_train": len(train),
        "n_test": len(test),
        "accuracy": accuracy_score(labels[test], pred),
        "balanced_accuracy": balanced_accuracy_score(labels[test], pred),
        "macro_f1": f1_score(labels[test], pred, average="macro", zero_division=0),
        "majority_accuracy_test": max(Counter(labels[test]).values()) / len(test),
        "n_features": len(vocabulary),
    }


def hub_statistics(articles: Sequence[ArticleExtraction], label_attr: str = "disease_label") -> dict[str, Any]:
    edge_counter = Counter(g for a in articles for g in set(a.genes))
    total_edges = sum(edge_counter.values())
    labels = [getattr(a, label_attr) for a in articles]
    class_values = sorted({x for x in labels if x})
    gene_disease_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for a in articles:
        label = getattr(a, label_attr)
        if not label:
            continue
        for g in set(a.genes):
            gene_disease_counts[g][label] += 1
    top_rows = []
    for gene, count in edge_counter.most_common(100):
        dist = np.asarray([gene_disease_counts[gene].get(c, 0) for c in class_values], dtype=float)
        nz = dist[dist > 0]
        ent = float(entropy(nz, base=2)) if len(nz) > 1 else 0.0
        max_ent = math.log2(len(class_values)) if len(class_values) > 1 else 1.0
        top_rows.append({
            "gene": gene,
            "article_count": count,
            "article_share": count / len(articles) if articles else float("nan"),
            "edge_share": count / total_edges if total_edges else float("nan"),
            "disease_count": int(np.count_nonzero(dist)),
            "disease_entropy": ent,
            "normalized_disease_entropy": ent / max_ent if max_ent else 0.0,
        })
    return {
        "n_articles": len(articles),
        "total_article_gene_edges": total_edges,
        "unique_genes": len(edge_counter),
        "top_10_edge_share": sum(v for _, v in edge_counter.most_common(10)) / total_edges if total_edges else float("nan"),
        "top_20_edge_share": sum(v for _, v in edge_counter.most_common(20)) / total_edges if total_edges else float("nan"),
        "gene_frequency_gini": gini(list(edge_counter.values())),
        "genes_in_at_least_10_percent": sum(v >= math.ceil(0.10 * len(articles)) for v in edge_counter.values()),
        "top_genes": top_rows,
    }


def analyze_label_level(
    all_articles: Sequence[ArticleExtraction],
    *,
    label_attr: str,
    min_class_n: int,
    out_dir: Path,
    full: bool,
) -> dict[str, Any]:
    articles = filter_articles_for_label(
        all_articles,
        label_attr=label_attr,
        min_class_n=min_class_n,
        confidence={"high"},
    )
    if len(articles) < 20 or len({getattr(a, label_attr) for a in articles}) < 2:
        return {"status": "insufficient", "n_articles": len(articles), "label_attr": label_attr}
    labels = np.asarray([getattr(a, label_attr) for a in articles], dtype=object)
    list_sizes = np.asarray([len(a.genes) for a in articles], dtype=int)
    years = np.asarray([a.year or 0 for a in articles], dtype=int)
    matrix, vocabulary = build_incidence(articles)
    sim = jaccard_matrix(matrix)
    n_perm = 5000 if full else 1000
    n_boot = 1000 if full else 300
    unrestricted = permutation_test_delta(sim, labels, n_perm=n_perm, seed=BASE_SEED)
    stratified = permutation_test_delta(
        sim,
        labels,
        n_perm=n_perm,
        seed=BASE_SEED + 1,
        stratified=True,
        list_sizes=list_sizes,
        years=years,
    )
    boot_ci = bootstrap_delta([a.genes for a in articles], labels.tolist(), n_boot=n_boot, seed=BASE_SEED + 2)

    sensitivities = []
    for remove_n in [0, 5, 10, 20]:
        if remove_n == 0:
            genesets = [a.genes for a in articles]
            removed = []
        else:
            genesets, removed = remove_top_genes(articles, remove_n)
        keep_idx = [i for i, gs in enumerate(genesets) if len(gs) >= 2]
        mat, _ = build_incidence([], genes_override=[genesets[i] for i in keep_idx])
        sm = jaccard_matrix(mat)
        lab = labels[keep_idx]
        result = permutation_test_delta(sm, lab, n_perm=(2000 if full else 500), seed=BASE_SEED + 100 + remove_n)
        result.update({"remove_top_n": remove_n, "removed_genes": removed, "n_articles": len(keep_idx)})
        sensitivities.append(result)

    freq_null = frequency_preserving_null(
        articles,
        label_attr,
        n_random=(300 if full else 50),
        swaps_per_edge=5,
        seed=BASE_SEED + 3,
    )

    class_counts = Counter(labels.tolist())
    classification_articles = [a for a in articles if class_counts[getattr(a, label_attr)] >= (10 if full else 6)]
    recoverability = None
    recoverability_no_hubs = None
    if len(classification_articles) >= 30 and len({getattr(a, label_attr) for a in classification_articles}) >= 3:
        top10 = set(g for g, _ in Counter(g for a in classification_articles for g in set(a.genes)).most_common(10))
        min_count = min(Counter(getattr(a, label_attr) for a in classification_articles).values())
        splits = max(2, min(5, min_count))
        recoverability = disease_recoverability(
            classification_articles,
            label_attr,
            repeats=(10 if full else 3),
            splits=splits,
            n_label_permutations=(200 if full else 30),
            seed=BASE_SEED + 4,
        )
        recoverability_no_hubs = disease_recoverability(
            classification_articles,
            label_attr,
            repeats=(10 if full else 3),
            splits=splits,
            n_label_permutations=(200 if full else 30),
            seed=BASE_SEED + 5,
            remove_genes=top10,
        )

    temporal = temporal_holdout(articles, label_attr, cutoff_year=2022, seed=BASE_SEED + 6)

    # Pairwise table, sampled for artifact size when very large.
    iu = np.triu_indices(len(articles), k=1)
    pair_rows = []
    for i, j, value in zip(iu[0], iu[1], sim[iu]):
        pair_rows.append({
            "pmcid_a": articles[int(i)].pmcid,
            "pmcid_b": articles[int(j)].pmcid,
            "label_a": labels[int(i)],
            "label_b": labels[int(j)],
            "same_label": bool(labels[int(i)] == labels[int(j)]),
            "jaccard": float(value),
            "n_genes_a": len(articles[int(i)].genes),
            "n_genes_b": len(articles[int(j)].genes),
        })
    if len(pair_rows) > 500_000:
        pair_rows.sort(key=lambda r: stable_int(f"{r['pmcid_a']}|{r['pmcid_b']}"))
        pair_rows = pair_rows[:500_000]
    write_csv(out_dir / "tables" / f"pairwise_{label_attr}.csv", pair_rows)

    return {
        "status": "completed",
        "label_attr": label_attr,
        "min_class_n": min_class_n,
        "n_articles": len(articles),
        "n_classes": len(class_counts),
        "class_counts": dict(sorted(class_counts.items(), key=lambda kv: (-kv[1], kv[0]))),
        "n_unique_genes": len(vocabulary),
        "primary_unrestricted": unrestricted,
        "primary_bootstrap_95_ci": list(boot_ci),
        "stratified_list_size_year": stratified,
        "sensitivities_remove_hubs": sensitivities,
        "frequency_preserving_null": freq_null,
        "recoverability": recoverability,
        "recoverability_remove_top10": recoverability_no_hubs,
        "temporal_holdout": temporal,
        "article_pmcids": [a.pmcid for a in articles],
    }


# ---------------------------------------------------------------------------
# Figures and report
# ---------------------------------------------------------------------------
def save_figures(
    all_articles: Sequence[ArticleExtraction],
    specific: Mapping[str, Any],
    broad: Mapping[str, Any],
    hub: Mapping[str, Any],
    out_dir: Path,
) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Top genes.
    top = hub.get("top_genes", [])[:20]
    if top:
        fig, ax = plt.subplots(figsize=(9, 6))
        names = [r["gene"] for r in top][::-1]
        counts = [r["article_count"] for r in top][::-1]
        ax.barh(names, counts)
        ax.set_xlabel("Number of included articles")
        ax.set_title("Most frequently reported hub/core/key genes")
        fig.tight_layout()
        fig.savefig(fig_dir / "figure_2_top_genes.png", dpi=220)
        plt.close(fig)

    # Flow diagram as text boxes.
    flow = json.loads((out_dir / "results" / "flow.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(figsize=(8.5, 8))
    ax.axis("off")
    items = [
        ("Europe PMC query", flow.get("query_hits", 0)),
        ("Metadata after deduplication", flow.get("metadata_deduplicated", 0)),
        ("Selected for full-text XML", flow.get("selected_for_xml", 0)),
        ("XML downloaded", flow.get("xml_downloaded", 0)),
        ("Explicit target list extracted", flow.get("target_lists_extracted", 0)),
        ("High-confidence target list", flow.get("high_confidence", 0)),
        ("Specific disease label mapped", flow.get("specific_disease_mapped", 0)),
        ("Primary specific-label analysis", specific.get("n_articles", 0)),
    ]
    ys = np.linspace(0.93, 0.08, len(items))
    for (label, value), y in zip(items, ys):
        ax.text(0.5, y, f"{label}\nN = {value:,}", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.5", edgecolor="black", facecolor="white"), fontsize=11)
        if y != ys[-1]:
            ax.annotate("", xy=(0.5, y - 0.065), xytext=(0.5, y - 0.025), arrowprops=dict(arrowstyle="->"))
    ax.set_title("Corpus construction and analysis flow", fontsize=14)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure_1_flow.png", dpi=220)
    plt.close(fig)

    # Sensitivity deltas.
    if specific.get("status") == "completed":
        sens = specific["sensitivities_remove_hubs"]
        labels = [f"Remove top {r['remove_top_n']}" for r in sens]
        deltas = [r["delta"] for r in sens]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(labels, deltas)
        ax.axhline(0, linewidth=1)
        ax.set_ylabel("Within-disease minus between-disease Jaccard")
        ax.set_title("Disease-specificity sensitivity to generic hub removal")
        ax.tick_params(axis="x", rotation=25)
        fig.tight_layout()
        fig.savefig(fig_dir / "figure_3_sensitivity.png", dpi=220)
        plt.close(fig)

    # Hub frequency vs entropy.
    if top:
        rows = hub.get("top_genes", [])
        fig, ax = plt.subplots(figsize=(8, 6))
        x = [r["article_share"] for r in rows]
        y = [r["normalized_disease_entropy"] for r in rows]
        ax.scatter(x, y, alpha=0.7)
        for r in rows[:12]:
            ax.annotate(r["gene"], (r["article_share"], r["normalized_disease_entropy"]), fontsize=8)
        ax.set_xlabel("Article prevalence")
        ax.set_ylabel("Normalized disease entropy")
        ax.set_title("Frequently reused genes spread across disease labels")
        fig.tight_layout()
        fig.savefig(fig_dir / "figure_4_gene_generality.png", dpi=220)
        plt.close(fig)


def fmt(x: Any, digits: int = 4) -> str:
    if x is None:
        return "NA"
    if isinstance(x, (float, np.floating)):
        if math.isnan(float(x)):
            return "NA"
        return f"{float(x):.{digits}f}"
    return str(x)


def generate_report(
    *,
    out_dir: Path,
    flow: Mapping[str, Any],
    specific: Mapping[str, Any],
    broad: Mapping[str, Any],
    hub: Mapping[str, Any],
    all_articles: Sequence[ArticleExtraction],
    full: bool,
) -> str:
    mode = "full frozen-corpus analysis" if full else "deterministic pilot"
    sp = specific.get("primary_unrestricted", {}) if specific.get("status") == "completed" else {}
    st = specific.get("stratified_list_size_year", {}) if specific.get("status") == "completed" else {}
    fn = specific.get("frequency_preserving_null", {}) if specific.get("status") == "completed" else {}
    rec = specific.get("recoverability") or {}
    rec_no = specific.get("recoverability_remove_top10") or {}
    boot = specific.get("primary_bootstrap_95_ci", [None, None])
    top_genes = hub.get("top_genes", [])[:20]
    top_table = "\n".join(
        f"| {i} | {r['gene']} | {r['article_count']} | {100*r['article_share']:.1f}% | {r['disease_count']} | {r['normalized_disease_entropy']:.3f} |"
        for i, r in enumerate(top_genes, start=1)
    ) or "| – | – | – | – | – | – |"
    class_table = "\n".join(
        f"| {label} | {n} |" for label, n in list((specific.get("class_counts") or {}).items())[:40]
    ) or "| – | – |"

    conclusion = "The primary analysis could not be completed because too few repeatedly represented disease labels were extracted."
    if specific.get("status") == "completed":
        if sp.get("p_one_sided_greater", 1) <= 0.05 and sp.get("delta", 0) > 0 and fn.get("p_one_sided_greater", 1) <= 0.05:
            conclusion = (
                "Author-reported target sets retained a statistically detectable disease-specific signal under both "
                "label permutation and the frequency-preserving incidence null. The magnitude and recoverability results "
                "determine whether that signal is practically useful."
            )
        elif sp.get("delta", 0) <= 0 or sp.get("p_one_sided_greater", 1) > 0.05:
            conclusion = (
                "Author-reported hub/core/key target sets did not show a reliable within-disease similarity advantage "
                "over between-disease pairs under the primary label-permutation test. Under the tested extraction and "
                "label definitions, the published target lists therefore do not encode a reproducible disease-specific fingerprint."
            )
        else:
            conclusion = (
                "A small within-disease signal was observed, but it did not survive all predeclared null controls; "
                "the evidence is insufficient to call the reported target sets disease-specific."
            )

    report = f"""# Does the disease label matter in herbal network pharmacology?

## Final computational audit report

**Analysis date:** {SNAPSHOT_DATE}  
**Mode:** {mode}  
**Data source:** Europe PMC Open Access full-text XML  
**Primary unit:** one article-level, author-reported hub/core/key target-gene set  
**Author:** Byungwoong Yoo, Independent Researcher

## Executive finding

{conclusion}

The conclusion is deliberately bounded. This study evaluates whether the **reported computational target lists** carry recoverable information about the named disease. It does not determine whether an herb is clinically effective, whether a listed target is causal, or whether any individual formula is safe.

## Frozen question and estimand

For article `i`, let `G_i` be the explicit hub/core/key gene set and `Y_i` its disease label. The primary statistic was:

```text
Delta_J = mean[J(G_i,G_j) | Y_i = Y_j] - mean[J(G_i,G_j) | Y_i != Y_j]
J(A,B) = |A intersect B| / |A union B|
```

A positive, non-random `Delta_J` is the minimum expected signature of disease specificity.

## Corpus flow

| Stage | Count |
|---|---:|
| Europe PMC query hits | {flow.get('query_hits', 0):,} |
| Metadata after deduplication | {flow.get('metadata_deduplicated', 0):,} |
| Selected for XML | {flow.get('selected_for_xml', 0):,} |
| XML downloaded | {flow.get('xml_downloaded', 0):,} |
| Explicit target lists extracted | {flow.get('target_lists_extracted', 0):,} |
| High-confidence lists | {flow.get('high_confidence', 0):,} |
| Specific disease mapped | {flow.get('specific_disease_mapped', 0):,} |
| Primary specific-label analysis | {specific.get('n_articles', 0):,} |
| Specific disease classes in primary analysis | {specific.get('n_classes', 0):,} |

The frozen query was:

```text
{EUROPE_PMC_QUERY}
```

Reviews, meta-analyses, bibliometric studies, editorials, protocols, and articles without an explicit finite target list were excluded by predeclared rules.

## Primary result: specific disease labels

| Quantity | Result |
|---|---:|
| Articles | {specific.get('n_articles', 'NA')} |
| Disease classes | {specific.get('n_classes', 'NA')} |
| Within-disease mean Jaccard | {fmt(sp.get('within_mean_jaccard'))} |
| Between-disease mean Jaccard | {fmt(sp.get('between_mean_jaccard'))} |
| `Delta_J` | {fmt(sp.get('delta'))} |
| Stratified-bootstrap 95% CI | [{fmt(boot[0])}, {fmt(boot[1])}] |
| Label permutations | {sp.get('n_permutations', 'NA')} |
| One-sided permutation p | {fmt(sp.get('p_one_sided_greater'))} |
| Two-sided permutation p | {fmt(sp.get('p_two_sided'))} |

### List-size and publication-era stratified permutation

| Quantity | Result |
|---|---:|
| Stratified `Delta_J` | {fmt(st.get('delta'))} |
| One-sided p | {fmt(st.get('p_one_sided_greater'))} |
| Null 95% interval | {st.get('null_95_interval', 'NA')} |

### Frequency-preserving incidence null

This null randomizes article–gene edges while preserving every article's list length and every gene's total article frequency. It therefore asks whether disease grouping explains more overlap than the generic hub-frequency structure alone.

| Quantity | Result |
|---|---:|
| Observed `Delta_J` | {fmt(fn.get('observed_delta'))} |
| Randomized incidence matrices | {fn.get('n_random_matrices', 'NA')} |
| Null mean | {fmt(fn.get('null_mean'))} |
| Null 95% interval | {fn.get('null_95_interval', 'NA')} |
| One-sided p | {fmt(fn.get('p_one_sided_greater'))} |

## Can an AI recover the disease from the target genes alone?

A linear support-vector classifier received only the article-by-gene binary matrix. TF–IDF weighting was fitted inside the cross-validation pipeline. No title, herb name, journal, year, pathway name, or abstract text was supplied.

| Metric | All genes | After removing top 10 recurrent genes |
|---|---:|---:|
| Articles | {rec.get('n_articles', 'NA')} | {rec_no.get('n_articles', 'NA')} |
| Classes | {rec.get('n_classes', 'NA')} | {rec_no.get('n_classes', 'NA')} |
| Accuracy | {fmt((rec.get('observed_mean') or {}).get('accuracy'))} | {fmt((rec_no.get('observed_mean') or {}).get('accuracy'))} |
| Balanced accuracy | {fmt((rec.get('observed_mean') or {}).get('balanced_accuracy'))} | {fmt((rec_no.get('observed_mean') or {}).get('balanced_accuracy'))} |
| Macro-F1 | {fmt((rec.get('observed_mean') or {}).get('macro_f1'))} | {fmt((rec_no.get('observed_mean') or {}).get('macro_f1'))} |
| Balanced-accuracy chance | {fmt(rec.get('balanced_accuracy_chance'))} | {fmt(rec_no.get('balanced_accuracy_chance'))} |
| Label-permutation p | {fmt(rec.get('permutation_p_balanced_accuracy'))} | {fmt(rec_no.get('permutation_p_balanced_accuracy'))} |

Classification above chance would indicate some disease information, but not necessarily biological specificity: journals, databases, formula traditions, or repeated pipelines can create a learnable signature. Failure to exceed chance is stronger evidence that the reported lists are label-invariant.

## Recurrent generic targets

| Rank | Gene | Articles | Article share | Disease labels | Normalized disease entropy |
|---:|---|---:|---:|---:|---:|
{top_table}

Summary concentration measures:

- Total article–gene edges: **{hub.get('total_article_gene_edges', 'NA')}**
- Unique genes: **{hub.get('unique_genes', 'NA')}**
- Top-10 edge share: **{100*hub.get('top_10_edge_share', float('nan')):.1f}%**
- Top-20 edge share: **{100*hub.get('top_20_edge_share', float('nan')):.1f}%**
- Gene-frequency Gini: **{fmt(hub.get('gene_frequency_gini'), 3)}**
- Genes reported in at least 10% of articles: **{hub.get('genes_in_at_least_10_percent', 'NA')}**

High frequency plus high disease entropy identifies a gene that is repeatedly reported across many unrelated disease labels—the pattern expected for a generic network hub rather than a disease fingerprint.

## Disease representation

| Disease label | Articles |
|---|---:|
{class_table}

## What was solved, and what was not

### Solved by this audit

The analysis provides a reproducible answer to the following bounded question:

> Do explicit author-reported hub/core/key gene lists from open-access herbal network-pharmacology papers contain a statistically recoverable disease-specific signature under label permutation, list-size/year stratification, recurrent-hub removal, and an article–gene frequency-preserving null?

Every included article, target symbol, source span, exclusion reason, code path, random seed, and file checksum is in the accompanying package.

### Not established

- clinical efficacy or inefficacy of any herb or formula;
- causal involvement of a target;
- biochemical binding;
- product identity, dose, bioavailability, or safety;
- validity of network pharmacology outside the frozen corpus and extraction contract.

## Limitations

1. Target lists were machine-extracted from XML. The primary set used strict, predeclared high-confidence rules, but machine extraction can still miss lists embedded only in images or supplementary files.
2. Disease labels came from a frozen title lexicon. Unmapped or multiply framed diseases were omitted from the specific-label analysis.
3. The study evaluates what articles explicitly reported, not every target generated internally by their pipelines.
4. Open-access indexing and title wording create coverage selection.
5. Repeated use of the same public databases can be a property of the field rather than misconduct by individual authors.
6. A non-specific target list does not imply that the intervention has no biological or clinical effect.

## Reproducibility files

- `tables/articles.csv`: one row per extracted article.
- `tables/article_gene_edges.csv`: one row per article–gene edge.
- `tables/exclusions.csv`: all recorded exclusions.
- `tables/extraction_candidates.csv`: candidate lists and exact source spans.
- `results/results.json`: complete machine-readable results.
- `results/flow.json`: corpus flow.
- `figures/`: publication-ready plots.
- `MANIFEST.json`: SHA-256 hashes and sizes.
- `audit.py`: complete executable pipeline.

## Final bounded conclusion

{conclusion}
"""
    path = out_dir / "FINAL_REPORT.md"
    path.write_text(report, encoding="utf-8")
    return report


def generate_manuscript(report: str, out_dir: Path, specific: Mapping[str, Any], hub: Mapping[str, Any]) -> None:
    # A compact journal-style draft; numerical values are inherited from the locked report.
    sp = specific.get("primary_unrestricted", {}) if specific.get("status") == "completed" else {}
    title = "Does the Disease Label Matter? A Corpus-Scale Specificity Audit of Herbal Network Pharmacology"
    manuscript = f"""# {title}

**Byungwoong Yoo**  
Independent Researcher, Republic of Korea

## Abstract

**Background:** Herbal network-pharmacology studies commonly interpret hub or core targets as disease mechanisms, but the disease specificity of those target sets has not been tested at corpus scale.  
**Methods:** We froze a Europe PMC Open Access query through {SNAPSHOT_DATE}, extracted explicit author-reported hub/core/key target lists from article XML using an NCBI human-gene dictionary, and mapped diseases from titles with a predeclared lexicon. The primary statistic compared mean within-disease and between-disease Jaccard similarity. Significance was assessed by disease-label permutation, a list-size/publication-era-stratified permutation, and a bipartite incidence null preserving article list lengths and gene frequencies. A cross-validated linear classifier tested whether disease labels could be recovered from gene sets alone.  
**Results:** The primary analysis included {specific.get('n_articles', 'NA')} articles across {specific.get('n_classes', 'NA')} repeatedly represented disease labels. Mean Jaccard similarity was {fmt(sp.get('within_mean_jaccard'))} within diseases and {fmt(sp.get('between_mean_jaccard'))} between diseases, for a difference of {fmt(sp.get('delta'))} (one-sided permutation p={fmt(sp.get('p_one_sided_greater'))}). The top 10 genes accounted for {100*hub.get('top_10_edge_share', float('nan')):.1f}% of all article–gene edges. Full null-control and classification results are reported in the accompanying tables.  
**Conclusions:** The result is restricted to the specificity of reported computational target lists and does not address clinical efficacy or causal mechanism. The corpus-level tests quantify whether disease labels leave a recoverable signature after controlling for generic hub reuse.

## Introduction

Network pharmacology is widely used to propose multi-component, multi-target mechanisms for herbal medicines. A recurring analytic sequence intersects predicted compound targets with disease-associated genes, constructs a protein–protein interaction network, and selects high-degree or central nodes as hub targets. The resulting targets are often discussed as disease-specific mechanistic explanations. However, highly annotated and highly connected genes can recur across many diseases, and the same databases and centrality algorithms may reproduce a generic set of inflammatory, apoptotic, and proliferative hubs.

This study asks a deliberately falsifiable question: if disease labels are removed, do the reported target lists retain enough structure to identify the disease? A disease-specific target set should be more similar to other studies of the same disease than to studies of different diseases, should outperform label-permuted controls, and should retain signal after article list length and global gene frequency are controlled.

## Methods

### Protocol and corpus

The query, date window, inclusion rules, disease lexicon, target-list anchors, random seed, and primary estimand were fixed in code before full-corpus outcome inspection. Europe PMC Open Access full-text XML was searched with the query reproduced in the final report. Reviews, meta-analyses, bibliometric studies, protocols, and papers without a finite explicit target list were excluded.

### Target extraction

Paragraphs, table cells, and captions containing predeclared phrases such as “hub genes,” “core targets,” and “key genes” were searched. Candidate gene symbols were normalized against the NCBI Homo sapiens gene information file. The primary analysis retained high-confidence candidates defined by an explicit count match or a dense gene list in an abstract, Results, Conclusion, or table context. The highest-scoring high-confidence set was selected per article. Exact source spans were retained for audit.

### Disease labels

Disease labels were assigned from article titles by a frozen, longest-match lexicon. Labels with fewer than the predeclared minimum number of articles were excluded from disease-level analysis but remained in corpus-wide hub-frequency summaries.

### Statistical analysis

Jaccard similarity was calculated for every article pair. The primary effect was the mean within-disease similarity minus the mean between-disease similarity. Label permutation preserved target lists and disease-group sizes. A stratified permutation operated within list-size quartiles and publication-era bins. A bipartite double-edge-swap null preserved each article's target-list length and each gene's total frequency. Sensitivity analyses removed the 5, 10, and 20 most frequently reported genes. Disease recoverability was evaluated with repeated stratified cross-validation using TF–IDF-weighted binary gene features and a class-weighted linear support-vector classifier. All tests used seed {BASE_SEED}.

## Results

See `FINAL_REPORT.md`, `results/results.json`, and the accompanying figures and tables. The numerical outputs in the Abstract were generated directly from the frozen analysis object, not entered manually.

## Discussion

The central interpretation depends on the direction and robustness of the corpus-level signal. A null result indicates that the reported hub/core/key lists do not constitute a reproducible disease fingerprint under the tested controls. A positive result indicates some disease-associated structure, but it remains necessary to distinguish biological specificity from shared databases, repeated pipelines, and publication subcultures. Generic hub concentration is therefore interpreted jointly with label-permutation, frequency-preserving randomization, hub-removal sensitivity, and out-of-sample disease recoverability.

The study does not evaluate whether any herbal intervention is effective. It evaluates whether a common computational evidence object—the reported target list—supports the disease-specific interpretation often attached to it.

## Data and code availability

All derived article-level data, source locators, exclusions, code, environment, random seeds, and checksums are included in the reproducibility package. Europe PMC full texts are not redistributed; the package contains identifiers and short extraction spans sufficient to locate each claim in the source article.
"""
    (out_dir / "MANUSCRIPT_DRAFT.md").write_text(manuscript, encoding="utf-8")


def build_manifest(out_dir: Path, script_path: Path) -> dict[str, Any]:
    manifest = {
        "project": "Herbal network-pharmacology disease specificity audit",
        "snapshot_date": SNAPSHOT_DATE,
        "base_seed": BASE_SEED,
        "generated_at_utc": utc_now(),
        "python": sys.version,
        "files": {},
    }
    for path in sorted(out_dir.rglob("*")):
        if path.is_file() and path.name != "MANIFEST.json":
            data = path.read_bytes()
            manifest["files"][str(path.relative_to(out_dir))] = {
                "bytes": len(data),
                "sha256": sha256_bytes(data),
            }
    manifest["source_script"] = {
        "path": str(script_path),
        "sha256": sha256_bytes(script_path.read_bytes()),
    }
    write_json(out_dir / "MANIFEST.json", manifest)
    return manifest


def package_outputs(out_dir: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for path in sorted(out_dir.rglob("*")):
            if path.is_file():
                zf.write(path, arcname=str(Path(out_dir.name) / path.relative_to(out_dir)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="artifact/herbal_specificity_audit")
    p.add_argument("--max-xml", type=int, default=None, help="Deterministic pilot sample; omit for full corpus")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--full", action="store_true", help="Use full permutation/bootstrap counts")
    p.add_argument("--min-specific-n", type=int, default=5)
    p.add_argument("--min-broad-n", type=int, default=10)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    for sub in ["raw", "tables", "results", "figures"]:
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    write_json(out_dir / "protocol.json", {
        "snapshot_date": SNAPSHOT_DATE,
        "query": EUROPE_PMC_QUERY,
        "base_seed": BASE_SEED,
        "max_xml": args.max_xml,
        "full_analysis": args.full,
        "min_specific_class_n": args.min_specific_n,
        "min_broad_class_n": args.min_broad_n,
        "title_exclusion_patterns": TITLE_EXCLUDE_PATTERNS,
        "anchors": [{"name": n, "pattern": p, "weight": w} for n, p, w in ANCHOR_PATTERNS],
        "disease_lexicon": [{"label": l, "broad": b, "patterns": p} for l, b, p in DISEASE_LEXICON],
    })

    metadata, meta_manifest = harvest_metadata(out_dir)
    selected, metadata_exclusions = select_metadata_for_xml(metadata, max_xml=args.max_xml)
    xml_manifest = download_xmls(selected, out_dir, workers=args.workers)
    downloaded_pmcids = {r["pmcid"] for r in xml_manifest if r.get("status") != "failed"}
    selected_downloaded = [r for r in selected if r.get("pmcid") in downloaded_pmcids]

    gene_dict = build_gene_dictionary(out_dir)
    extractions, extraction_exclusions, candidate_rows = extract_all(selected_downloaded, out_dir, gene_dict)
    exclusions = metadata_exclusions + extraction_exclusions

    article_rows = []
    edge_rows = []
    for a in extractions:
        row = asdict(a)
        row["genes"] = "|".join(a.genes)
        row["gene_count"] = len(a.genes)
        article_rows.append(row)
        for gene in a.genes:
            edge_rows.append({
                "pmcid": a.pmcid,
                "disease_label": a.disease_label,
                "broad_disease": a.broad_disease,
                "gene": gene,
                "confidence": a.extraction_confidence,
            })
    write_csv(out_dir / "tables" / "articles.csv", article_rows)
    write_csv(out_dir / "tables" / "article_gene_edges.csv", edge_rows)
    write_csv(out_dir / "tables" / "exclusions.csv", exclusions)
    write_csv(out_dir / "tables" / "extraction_candidates.csv", candidate_rows)

    flow = {
        "query_hits": meta_manifest.get("reported_hit_count"),
        "metadata_deduplicated": len(metadata),
        "metadata_excluded": len(metadata_exclusions),
        "selected_for_xml": len(selected),
        "xml_downloaded": len(selected_downloaded),
        "xml_failed": len(selected) - len(selected_downloaded),
        "target_lists_extracted": len(extractions),
        "high_confidence": sum(a.extraction_confidence == "high" for a in extractions),
        "medium_confidence": sum(a.extraction_confidence == "medium" for a in extractions),
        "specific_disease_mapped": sum(a.disease_label is not None for a in extractions),
        "broad_disease_mapped": sum(a.broad_disease is not None for a in extractions),
        "pilot_max_xml": args.max_xml,
    }
    write_json(out_dir / "results" / "flow.json", flow)

    high_articles = [a for a in extractions if a.extraction_confidence == "high"]
    hub = hub_statistics(high_articles, label_attr="disease_label")
    specific = analyze_label_level(
        extractions,
        label_attr="disease_label",
        min_class_n=args.min_specific_n,
        out_dir=out_dir,
        full=args.full,
    )
    broad = analyze_label_level(
        extractions,
        label_attr="broad_disease",
        min_class_n=args.min_broad_n,
        out_dir=out_dir,
        full=args.full,
    )
    results = {
        "flow": flow,
        "hub_statistics": hub,
        "specific_disease_analysis": specific,
        "broad_disease_analysis": broad,
    }
    write_json(out_dir / "results" / "results.json", results)
    write_csv(out_dir / "tables" / "top_genes.csv", hub.get("top_genes", []))

    save_figures(extractions, specific, broad, hub, out_dir)
    report = generate_report(
        out_dir=out_dir,
        flow=flow,
        specific=specific,
        broad=broad,
        hub=hub,
        all_articles=extractions,
        full=args.full,
    )
    generate_manuscript(report, out_dir, specific, hub)

    # Copy executable source and environment lock hints into package.
    script_path = Path(__file__).resolve()
    shutil.copy2(script_path, out_dir / "audit.py")
    (out_dir / "requirements.txt").write_text(
        "numpy\npandas\nscipy\nscikit-learn\nlxml\nmatplotlib\n",
        encoding="utf-8",
    )
    build_manifest(out_dir, script_path)
    package_outputs(out_dir, out_dir.parent / "herbal_specificity_audit_package.zip")

    print(json.dumps({
        "flow": flow,
        "specific_primary": specific.get("primary_unrestricted"),
        "recoverability": specific.get("recoverability"),
        "top_10_edge_share": hub.get("top_10_edge_share"),
        "package": str(out_dir.parent / "herbal_specificity_audit_package.zip"),
    }, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
