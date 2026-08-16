#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import io
import json
import math
import os
import re
import shutil
import sys
import time
import traceback
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests
from lxml import etree
from scipy import stats

OUT = Path(os.environ.get("OUT_DIR", "smoke_output"))
RAW = OUT / "raw"
TABLES = OUT / "tables"
LOGS = OUT / "logs"
for p in (OUT, RAW, TABLES, LOGS):
    p.mkdir(parents=True, exist_ok=True)

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "OpenProblemTournament/0.1 (independent reproducibility research; contact yoonge3@gmail.com)"})

@dataclass
class TrackResult:
    track: str
    question: str
    status: str
    gate: str
    headline: str
    metrics: dict[str, Any]
    limitations: list[str]
    next_decision: str
    files: list[str]


def log(msg: str) -> None:
    print(msg, flush=True)


def get(url: str, *, timeout: int = 120, retries: int = 4, stream: bool = False) -> requests.Response:
    last = None
    for k in range(retries):
        try:
            r = SESSION.get(url, timeout=timeout, stream=stream)
            if r.status_code == 429:
                time.sleep(2 ** k)
                continue
            r.raise_for_status()
            return r
        except Exception as e:
            last = e
            time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"GET failed after {retries}: {url}: {last}")


def download(url: str, path: Path, *, min_bytes: int = 1) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size >= min_bytes:
        return path
    r = get(url, stream=True, timeout=240)
    with path.open("wb") as f:
        for chunk in r.iter_content(1024 * 1024):
            if chunk:
                f.write(chunk)
    if path.stat().st_size < min_bytes:
        raise RuntimeError(f"Downloaded file too small: {url}, {path.stat().st_size}")
    return path


def sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def clean_text(x: Any) -> str:
    return re.sub(r"\s+", " ", str(x or "")).strip()


def safe_json(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): safe_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [safe_json(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.bool_,)):
 return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def parse_soft_samples(path: Path) -> pd.DataFrame:
    rows = []
    current: dict[str, Any] | None = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    rows.append(current)
                current = {"gsm": line.split("=", 1)[1].strip(), "characteristics": []}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_source_name_ch1 = "):
                current["source"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                current["characteristics"].append(line.split("=", 1)[1].strip())
            elif current is not None and line.startswith("!Sample_description = "):
                current.setdefault("description", []).append(line.split("=", 1)[1].strip())
    if current:
        rows.append(current)
    for r in rows:
        r["characteristics_text"] = " | ".join(r.get("characteristics", []))
        r["description_text"] = " | ".join(r.get("description", []))
    return pd.DataFrame(rows)


def syndrome_labels(text: str) -> list[str]:
    t = clean_text(text).lower().replace("-", " ")
    patterns = [
        ("qi_stagnation_blood_stasis", ["qi stagnation and blood stasis", "qi stagnation blood stasis", "qsbss", "qsbs"]),
        ("qi_deficiency_blood_stasis", ["qi deficiency and blood stasis", "qi deficiency blood stasis", "qdbs"]),
        ("yang_deficiency_blood_stasis", ["yang deficiency and blood stasis", "yang deficiency blood stasis", "ydbs"]),
        ("blood_stasis", ["blood stasis syndrome", "blood stasis", "bss"]),
        ("blood_heat", ["blood heat", "bhs"]),
        ("blood_dryness", ["blood dryness", "bds"]),
        ("cold_dampness", ["cold dampness", "cold-dampness"]),
        ("liver_kidney_deficiency", ["liver kidney deficiency", "liver-kidney deficiency"]),
        ("healthy", ["healthy", "normal control", "control"]),
    ]
    found = []
    for label, aliases in patterns:
        if any(a in t for a in aliases):
            found.append(label)
    if "blood_stasis" in found and any(x.endswith("_blood_stasis") for x in found if x != "blood_stasis"):
        found.remove("blood_stasis")
    return found


def extract_geo_dir_file(gse: str, filename: str, local_name: str | None = None) -> Path:
    prefix = re.sub(r"\d{3}$", "nnn", gse)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/{quote(filename)}"
    return download(url, RAW / (local_name or filename), min_bytes=100)


def extract_geo_soft(gse: str) -> Path:
    prefix = re.sub(r"\d{3}$", "nnn", gse)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/soft/{gse}_family.soft.gz"
    return download(url, RAW / f"{gse}_family.soft.gz", min_bytes=100)


def inspect_dataframe_file(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
    try:
        if path.suffix.lower() in {".xlsx", ".xls"}:
            xls = pd.ExcelFile(path)
            out["sheets"] = []
            for s in xls.sheet_names:
                df = pd.read_excel(path, sheet_name=s, nrows=10)
                out["sheets"].append({"name": s, "preview_rows": len(df), "columns": [str(c) for c in df.columns[:100]]})
        else:
            compression = "gzip" if path.suffix == ".gz" else "infer"
            df = pd.read_csv(path, sep=None, engine="python", compression=compression, nrows=10)
            out["columns"] = [str(c) for c in df.columns[:100]]
            out["preview_rows"] = len(df)
    except Exception as e:
        out["inspection_error"] = repr(e)
    return out


def track_pattern_biology() -> TrackResult:
    track_dir = OUT / "01_pattern_biology"
    track_dir.mkdir(exist_ok=True)
    inventory = []
    errors = []
    datasets = [
        {"gse": "GSE109265", "disease": "diabetes", "processed": "GSE109265_Processed_data.xlsx.gz", "note": "Serum-treated HUVECs; six samples, one per phenotype group."},
        {"gse": "GSE192867", "disease": "psoriasis", "processed": "GSE192867_Processed_data-file_3.xlsx", "note": "Human PBMC/monocyte transcriptome; blood stasis/heat/dryness and controls."},
        {"gse": "GSE303117", "disease": "ischemic_heart_failure", "processed": "GSE303117_gene_fpkm.txt.gz", "note": "Human transcriptome; QDBS/YDBS/YDBSFR/healthy."},
    ]
    for d in datasets:
        try:
            soft = extract_geo_soft(d["gse"])
            meta = parse_soft_samples(soft)
            meta["all_text"] = meta.fillna("").astype(str).agg(" | ".join, axis=1)
            meta["syndrome_labels"] = meta["all_text"].map(lambda x: ";".join(syndrome_labels(x)))
            meta.to_csv(track_dir / f"{d['gse']}_sample_metadata.csv", index=False)
            counts = Counter()
            for x in meta["syndrome_labels"]:
                for y in str(x).split(";"):
                    if y:
                        counts[y] += 1
            proc_inspect = None
            try:
                proc_path = extract_geo_dir_file(d["gse"], d["processed"])
                proc_inspect = inspect_dataframe_file(proc_path)
            except Exception as pe:
                errors.append(f"{d['gse']} processed: {pe}")
            inventory.append({"gse": d["gse"], "disease": d["disease"], "n_samples": len(meta), "syndrome_counts": dict(counts), "processed": proc_inspect, "note": d["note"]})
        except Exception as e:
            errors.append(f"{d['gse']}: {e}")
            inventory.append({"gse": d["gse"], "disease": d["disease"], "error": repr(e), "note": d["note"]})
    label_disease_counts: dict[str, dict[str, int]] = defaultdict(dict)
    for row in inventory:
        if "syndrome_counts" not in row:
            continue
        for label, n in row["syndrome_counts"].items():
            label_disease_counts[label][row["disease"]] = n
    eligible_exact = {label: {d: n for d, n in dc.items() if n >= 10} for label, dc in label_disease_counts.items() if len({d for d, n in dc.items() if n >= 10}) >= 2 and label != "healthy"}
    (track_dir / "inventory.json").write_text(json.dumps(safe_json(inventory), indent=2), encoding="utf-8")
    metrics = {"datasets_attempted": len(datasets), "datasets_parsed": sum("n_samples" in x for x in inventory), "inventory": inventory, "exact_syndrome_labels_with_two_diseases_n_ge_10": eligible_exact, "errors": errors}
    if eligible_exact:
        status, gate, headline, next_decision = "PASS_DATA_GATE", "Exact cross-disease label with usable replication exists.", "At least one exact syndrome label can be tested across two diseases.", "Proceed to harmonized gene/pathway leave-one-disease-out testing."
    else:
        status, gate = "FAIL_EXACT_CROSS_DISEASE_DATA_GATE", "No exact syndrome label has at least 10 samples in two independent public human disease datasets."
        headline = "The scientifically strongest question is not yet solvable from the verified public datasets without relaxing syndrome identity."
        next_decision = "Do not force a cross-disease biomarker result. Preserve as a data-gap result; seek author-shared QSBSS multi-disease data or newly released cohorts."
    return TrackResult("1_cross_disease_pattern_biology", "Does the same traditional-medicine pattern represent a shared biological state across diseases?", status, gate, headline, metrics, ["Related labels such as generic blood stasis and qi-deficiency blood stasis are not interchangeable.", "GSE109265 has one sample per phenotype and uses serum-treated HUVECs rather than replicated patient transcriptomes.", "Within-disease separability cannot substitute for cross-disease transfer."], next_decision, [str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])


def xml_tables(pmcid: str, track_dir: Path) -> tuple[list[dict[str, Any]], etree._Element]:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    xml_path = download(url, RAW / f"{pmcid}.xml", min_bytes=1000)
    root = etree.parse(str(xml_path)).getroot()
    tables = []
    for idx, tw in enumerate(root.xpath(".//table-wrap"), 1):
        label = clean_text(" ".join(tw.xpath("./label//text()"))) or f"table_{idx}"
        caption = clean_text(" ".join(tw.xpath("./caption//text()")))
        rows = []
        for tr in tw.xpath(".//tr"):
            cells = [clean_text(" ".join(c.xpath(".//text()"))) for c in tr.xpath("./th|./td")]
            if cells:
                rows.append(cells)
        tables.append({"label": label, "caption": caption, "rows": rows})
    (track_dir / f"{pmcid}_tables.json").write_text(json.dumps(tables, ensure_ascii=False, indent=2), encoding="utf-8")
    return tables, root


def find_supplement_links(root: etree._Element) -> list[str]:
    hrefs = []
    for el in root.xpath(".//supplementary-material|.//media|.//inline-supplementary-material"):
        href = el.get("{http://www.w3.org/1999/xlink}href") or el.get("href")
        if href:
            hrefs.append(href)
    return sorted(set(hrefs))


def normalize_herb(text: str) -> str:
    t = re.sub(r"\([^)]*\)", "", text)
    t = re.sub(r"\b\d+(?:\.\d+)?\s*(?:g|mg|ml)\b", "", t, flags=re.I)
    t = re.sub(r"[^A-Za-z ]", " ", t)
    return re.sub(r"\s+", " ", t).strip().lower()


def track_formula_synergy() -> TrackResult:
    track_dir = OUT / "02_formula_synergy"; track_dir.mkdir(exist_ok=True)
    pmcid = "PMC8345858"
    try:
        tables, root = xml_tables(pmcid, track_dir)
    except Exception as e:
        return TrackResult("2_formula_synergy", "Do multi-herb formulas show effects beyond additive component effects?", "DOWNLOAD_OR_PARSE_FAILED", "Article XML unavailable", str(e), {}, ["No analysis was forced without source tables."], "Retry source acquisition.", [])
    ingredient_table = next((t for t in tables if "ingredient" in (t["caption"] + " " + t["label"]).lower()), None)
    formulas = []
    if ingredient_table:
        for r in ingredient_table["rows"][1:]:
            if not r: continue
            study, herbs = r[0], []
            for c in r[1:]:
                for p in re.split(r"\s{2,}|;|\|", clean_text(c)):
                    h = normalize_herb(p)
                    if len(h) >= 4 and not h.startswith(("quality assessment", "yes", "no")):
                        herbs.append(h)
            herbs = sorted(set(herbs))
            if herbs: formulas.append({"study": study, "herbs": herbs, "n_herbs": len(herbs)})
    unique_herbs = sorted({h for f in formulas for h in f["herbs"]})
    X = np.array([[int(h in f["herbs"]) for h in unique_herbs] for f in formulas], dtype=float) if formulas else np.empty((0,0))
    rank = int(np.linalg.matrix_rank(X)) if X.size else 0
    herb_counts = Counter(h for f in formulas for h in f["herbs"])
    pair_counts = Counter()
    for f in formulas:
        hs = sorted(set(f["herbs"]))
        for i in range(len(hs)):
            for j in range(i+1, len(hs)):
                pair_counts[(hs[i], hs[j])] += 1
    supported_pairs = [{"a": a, "b": b, "n_formulas": n} for (a,b),n in pair_counts.items() if n >= 3]
    outcome_keywords = re.compile(r"respond|effective|event|mean|sd|sample|risk ratio|odds ratio|total", re.I)
    outcome_tables = []
    has_outcomes = False
    for t in tables:
        flat = " ".join(" ".join(r) for r in t["rows"][:5])
        if outcome_keywords.search(flat): outcome_tables.append({"label": t["label"], "caption": t["caption"], "n_rows": len(t["rows"]), "preview": t["rows"][:4]})
        if len(t["rows"]) >= 5:
            header = " | ".join(t["rows"][0]).lower()
            if ("sample" in header or " n " in f" {header} ") and any(k in header for k in ["effective","response","mean","sd","events"]): has_outcomes = True
    (track_dir / "formula_components.json").write_text(json.dumps(formulas, ensure_ascii=False, indent=2), encoding="utf-8")
    (track_dir / "supported_pairs.json").write_text(json.dumps(supported_pairs, ensure_ascii=False, indent=2), encoding="utf-8")
    has_component_variation = len(formulas)>=10 and len(unique_herbs)>=10 and rank>=5
    metrics = {"pmcid": pmcid, "tables_found": len(tables), "ingredient_table_found": bool(ingredient_table), "formulas_parsed": len(formulas), "unique_herbs": len(unique_herbs), "component_matrix_rank": rank, "top_herbs": herb_counts.most_common(20), "pairs_supported_in_3plus_formulas": len(supported_pairs), "machine_readable_arm_outcomes_found": has_outcomes, "outcome_like_tables": outcome_tables, "supplement_links": find_supplement_links(root)}
    if has_component_variation and has_outcomes and len(supported_pairs)>=3:
        status, gate, headline, next_decision = "PASS_SMOKE", "Component variation and arm-level outcome data are sufficient for a pilot CNMA.", "A component-additivity versus interaction model can be fitted from this corpus.", "Fit preregistered additive and sparse-interaction CNMA with held-out trials."
    elif has_component_variation:
        status, gate = "PASS_STRUCTURE_FAIL_OUTCOME_DATA", "Formula composition is machine-readable, but article-level arm outcomes are not supplied in an analysis-ready table."
        headline, next_decision = "The synergy question is promising, but the smoke corpus requires trial-level outcome reconstruction before modeling.", "Reconstruct one disease corpus from primary RCTs or obtain the review extraction workbook; do not infer synergy from co-occurrence alone."
    else:
        status, gate, headline, next_decision = "FAIL_IDENTIFIABILITY_GATE", "Insufficient independent component variation.", "The selected corpus cannot identify additive and interaction effects separately.", "Try another disease with more formula variation and reusable trial data."
    return TrackResult("2_formula_synergy", "Do multi-herb formulas outperform the additive effects of their components?", status, gate, headline, metrics, ["Component co-occurrence is not synergy.", "Aggregate RCT data cannot resolve patient-level effect modification.", "Herb combinations are highly collinear and may be confounded by dose, formulation, and trial quality."], next_decision, [str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])


def track_acupoint_specificity() -> TrackResult:
    track_dir = OUT / "03_acupoint_specificity"; track_dir.mkdir(exist_ok=True)
    corpus, errors = [], []
    for pmcid in ["PMC4965798", "PMC12403578"]:
        try:
            tables, root = xml_tables(pmcid, track_dir); corpus.append({"pmcid":pmcid,"n_tables":len(tables),"supplement_links":find_supplement_links(root),"tables":tables})
        except Exception as e: errors.append(f"{pmcid}: {e}")
    complete_trial_rows=[]; feature_rows=0; effect_rows=0
    for art in corpus:
        for t in art["tables"]:
            for r in t["rows"][1:]:
                txt=" | ".join(r).lower(); has_feature=any(k in txt for k in ["nonpenetr","non-penetr","superficial","deep need","nonacupoint","non-acupoint","acupoint"]); has_effect=bool(re.search(r"(?:smd|md|mean difference|effect size|[-−]?\d+\.\d+\s*\(?[-−]?\d+\.\d+)",txt))
                feature_rows += int(has_feature); effect_rows += int(has_effect)
                if has_feature and has_effect: complete_trial_rows.append({"pmcid":art["pmcid"],"table":t["label"],"row":r})
    (track_dir / "complete_trial_rows.json").write_text(json.dumps(complete_trial_rows,ensure_ascii=False,indent=2),encoding="utf-8")
    metrics={"articles":[{"pmcid":x["pmcid"],"n_tables":x["n_tables"],"n_supplement_links":len(x["supplement_links"])} for x in corpus],"feature_rows":feature_rows,"effect_rows":effect_rows,"rows_with_joint_effect_and_sham_location_depth":len(complete_trial_rows),"errors":errors}
    if len(complete_trial_rows)>=20:
        status,gate,headline,next_decision="PASS_SMOKE","At least 20 machine-readable trial rows jointly contain effect and sham/acupoint features.","A trial-level meta-regression of location versus penetration can be executed.","Normalize STRICTA fields and fit multivariable random-effects meta-regression."
    else:
        status,gate="FAIL_MACHINE_READABLE_JOINT_DATA_GATE","Published open tables do not jointly expose enough trial effect sizes and acupoint/sham features."
        headline,next_decision="The scientific question is strong, but current open review outputs do not support a clean automated causal component analysis without reconstructing primary trials.","Build a primary-trial extraction set for one condition or obtain IPD; do not claim location specificity from aggregate subgroup summaries."
    return TrackResult("3_acupoint_specificity","Does correct acupoint location add benefit beyond penetration, stimulation, treatment dose, and sham design?",status,gate,headline,metrics,["Sham location and penetration are often bundled rather than independently randomized.","Published reviews frequently report subgroup summaries rather than reusable trial-level component data.","Acupoint selection is condition-specific and correlated with session dose and practitioner style."],next_decision,[str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])


def find_patient_level_sheet(path: Path):
    xls=pd.ExcelFile(path); info={"sheets":[]}; best=None; best_sheet=None; best_df=None
    for s in xls.sheet_names:
        df=pd.read_excel(path,sheet_name=s); info["sheets"].append({"sheet":s,"rows":len(df),"columns":[str(c) for c in df.columns]})
        if len(df)>=80 and df.shape[1]>=8:
            score=len(df)*min(df.shape[1],50)
            if best is None or score>best: best=score; best_sheet=s; best_df=df
    return best_sheet,best_df,info


def fuzzy_col(columns: Iterable[str], patterns:list[str]):
    cols=[str(c) for c in columns]
    for p in patterns:
        for c in cols:
            if p.lower()==c.lower().strip(): return c
    for p in patterns:
        for c in cols:
            if p.lower() in c.lower(): return c
    return None


def track_pattern_incremental_value() -> TrackResult:
    track_dir=OUT/"04_pattern_incremental_value"; track_dir.mkdir(exist_ok=True)
    url="https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12913-026-15077-x/MediaObjects/12913_2026_15077_MOESM1_ESM.xlsx"; path=RAW/"12913_2026_15077_MOESM1_ESM.xlsx"
    try:
        download(url,path,min_bytes=1000); best_sheet,df,info=find_patient_level_sheet(path); (track_dir/"supplement_inspection.json").write_text(json.dumps(safe_json(info),ensure_ascii=False,indent=2),encoding="utf-8")
    except Exception as e:
        return TrackResult("4_pattern_incremental_value","Does pattern diagnosis add out-of-sample predictive value beyond ordinary clinical predictors?","SUPPLEMENT_DOWNLOAD_OR_PARSE_FAILED","Could not inspect the official supplement.",str(e),{"url":url},["No incremental analysis was reconstructed from published AUCs alone."],"Retry official supplement retrieval or contact authors for anonymized row-level data.",[])
    metrics={"supplement_bytes":path.stat().st_size,"patient_level_sheet":best_sheet,"inspection":info}
    if df is None:
        return TrackResult("4_pattern_incremental_value","Does pattern diagnosis add out-of-sample predictive value beyond ordinary clinical predictors?","FAIL_PATIENT_LEVEL_DATA_GATE","Supplement contains no identifiable patient-level table.","The published supplement does not permit a direct clinical-only versus clinical-plus-pattern reanalysis.",metrics,["Published combined-model AUC cannot identify incremental value of pattern diagnosis."],"Request row-level data or find another open cohort.",[str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])
    df.to_csv(track_dir/"candidate_patient_level_sheet.csv",index=False)
    outcome=fuzzy_col(df.columns,["responder","response","non-responder","nonresponse","outcome","efficacy"]); pattern=fuzzy_col(df.columns,["TCM syndrome","syndrome pattern","pattern","syndrome"])
    clinical_candidates={"age":fuzzy_col(df.columns,["age"]),"sex":fuzzy_col(df.columns,["sex","gender"]),"bmi":fuzzy_col(df.columns,["BMI","body mass index"]),"duration":fuzzy_col(df.columns,["disease duration","duration"]),"baseline_vas":fuzzy_col(df.columns,["pre-treatment VAS","baseline VAS","VAS before","pre VAS"]),"radiating_pain":fuzzy_col(df.columns,["radiating leg pain","radiating pain"]),"previous_treatment":fuzzy_col(df.columns,["previous treatment","treatment history"])}
    clinical_cols=[c for c in clinical_candidates.values() if c]; metrics.update({"outcome_column":outcome,"pattern_column":pattern,"clinical_columns":clinical_cols,"n_rows":len(df),"n_columns":df.shape[1]})
    if not outcome or not pattern or len(clinical_cols)<2:
        status,gate,headline,next_decision="FAIL_VARIABLE_MAPPING_GATE","Patient-level-looking sheet lacks mappable outcome/pattern/clinical fields.","The supplement is not analysis-ready for the incremental-value question.","Manually inspect dictionary or obtain raw data."
    else:
        from sklearn.compose import ColumnTransformer
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import OneHotEncoder,StandardScaler,LabelEncoder
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import RepeatedStratifiedKFold
        from sklearn.metrics import roc_auc_score,brier_score_loss,log_loss
        sub=df[[outcome,pattern]+clinical_cols].copy().dropna(subset=[outcome]); yraw=sub[outcome]
        if yraw.nunique(dropna=True)!=2:
            status,gate,headline,next_decision="FAIL_BINARY_OUTCOME_GATE",f"Outcome has {yraw.nunique(dropna=True)} levels, not a clean binary endpoint.","Incremental prediction was not forced.","Define endpoint from raw VAS if available."
        else:
            le=LabelEncoder(); y=le.fit_transform(yraw.astype(str))
            def build(cols):
                num=[c for c in cols if pd.api.types.is_numeric_dtype(sub[c])]; cat=[c for c in cols if c not in num]
                pre=ColumnTransformer([("num",Pipeline([("imp",SimpleImputer(strategy="median")),("sc",StandardScaler())]),num),("cat",Pipeline([("imp",SimpleImputer(strategy="most_frequent")),("oh",OneHotEncoder(handle_unknown="ignore"))]),cat)])
                return Pipeline([("pre",pre),("model",LogisticRegression(max_iter=5000,class_weight="balanced"))])
            nmin=int(np.bincount(y).min()); n_splits=min(5,max(2,nmin)); cv=RepeatedStratifiedKFold(n_splits=n_splits,n_repeats=10,random_state=20260816); fold_rows=[]
            for fold,(tr,te) in enumerate(cv.split(sub,y),1):
                for name,cols in [("clinical",clinical_cols),("clinical_plus_pattern",clinical_cols+[pattern])]:
                    model=build(cols); model.fit(sub.iloc[tr][cols],y[tr]); p=model.predict_proba(sub.iloc[te][cols])[:,1]; fold_rows.append({"fold":fold,"model":name,"auc":roc_auc_score(y[te],p),"brier":brier_score_loss(y[te],p),"logloss":log_loss(y[te],p,labels=[0,1])})
            fr=pd.DataFrame(fold_rows); fr.to_csv(track_dir/"incremental_cv_folds.csv",index=False); piv=fr.pivot(index="fold",columns="model",values=["auc","brier","logloss"])
            arrays={"delta_auc":(piv["auc"]["clinical_plus_pattern"]-piv["auc"]["clinical"]).to_numpy(),"delta_brier_improvement":(piv["brier"]["clinical"]-piv["brier"]["clinical_plus_pattern"]).to_numpy(),"delta_logloss_improvement":(piv["logloss"]["clinical"]-piv["logloss"]["clinical_plus_pattern"]).to_numpy()}
            diffs={k:float(v.mean()) for k,v in arrays.items()}; rng=np.random.default_rng(20260816); boot={}
            for k,arr in arrays.items():
                vals=[np.mean(rng.choice(arr,size=len(arr),replace=True)) for _ in range(5000)]; boot[k+"_95ci"]=[float(np.quantile(vals,.025)),float(np.quantile(vals,.975))]
            metrics.update({"cv":diffs,"cv_intervals":boot,"outcome_levels":list(le.classes_),"n_analyzed":len(sub),"n_splits":n_splits,"n_repeats":10}); pass_all=boot["delta_auc_95ci"][0]>0 and boot["delta_brier_improvement_95ci"][0]>0
            if pass_all:
                status,gate,headline,next_decision="PROMISING_INCREMENTAL_SIGNAL","Pattern improved both AUC and Brier score in repeated leakage-safe cross-validation.","This is the first candidate that directly addresses whether pattern information adds predictive value beyond ordinary clinical variables.","Run external-site analysis, calibration, decision curves, and prespecified sensitivity models; test pattern-only versus conventional-only and treatment interactions."
            else:
                status,gate,headline,next_decision="NO_ROBUST_INCREMENTAL_SIGNAL","Pattern did not robustly improve both discrimination and proper scoring.","The open cohort does not support a strong incremental-value claim under this smoke analysis.","Audit outcome coding and external cohort; do not promote the combined model AUC as evidence for pattern value."
    return TrackResult("4_pattern_incremental_value","Does pattern diagnosis add out-of-sample predictive value beyond ordinary clinical predictors?",status,gate,headline,metrics,["Retrospective treatment data are vulnerable to confounding and treatment-selection leakage.","The published model mixes pretreatment patient features with acupoint and combined-treatment choices.","A single center cannot establish generalizable incremental value."],next_decision,[str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])


def openfda_count(search:str):
    try:
        r=SESSION.get("https://api.fda.gov/drug/event.json",params={"search":search,"limit":1},timeout=120)
        if r.status_code==404:return 0
        r.raise_for_status(); return int(r.json().get("meta",{}).get("results",{}).get("total",0))
    except Exception:return None

def qterm(field,text):return f'{field}:"{text}"'

def track_herb_drug_safety() -> TrackResult:
    track_dir=OUT/"05_herb_drug_safety"; track_dir.mkdir(exist_ok=True)
    pairs=[{"name":"ginkgo_warfarin_bleeding","herb_aliases":["GINKGO","GINKGO BILOBA"],"drug_aliases":["WARFARIN","COUMADIN"],"event_terms":["HAEMORRHAGE","GASTROINTESTINAL HAEMORRHAGE","INTERNATIONAL NORMALISED RATIO INCREASED"]},{"name":"st_johns_wort_ssri_serotonin","herb_aliases":["ST JOHNS WORT","HYPERICUM"],"drug_aliases":["SERTRALINE","FLUOXETINE","PAROXETINE","CITALOPRAM","ESCITALOPRAM"],"event_terms":["SEROTONIN SYNDROME"]},{"name":"licorice_diuretic_hypokalaemia","herb_aliases":["LICORICE","GLYCYRRHIZA"],"drug_aliases":["FUROSEMIDE","HYDROCHLOROTHIAZIDE"],"event_terms":["HYPOKALAEMIA","BLOOD POTASSIUM DECREASED"]}]
    date_windows={"early":"[20040101+TO+20191231]","late":"[20200101+TO+20260428]"}; results=[]
    for pair in pairs:
        for period,drange in date_windows.items():
            herb_q="("+"+OR+".join(qterm("patient.drug.medicinalproduct",x) for x in pair["herb_aliases"])+")"; drug_q="("+"+OR+".join(qterm("patient.drug.medicinalproduct",x) for x in pair["drug_aliases"])+")"; event_q="("+"+OR+".join(qterm("patient.reaction.reactionmeddrapt",x) for x in pair["event_terms"])+")"; date_q=f"receivedate:{drange}"
            queries={"a_exposed_event":f"{drug_q}+AND+{herb_q}+AND+{event_q}+AND+{date_q}","exposed_total":f"{drug_q}+AND+{herb_q}+AND+{date_q}","drug_event_total":f"{drug_q}+AND+{event_q}+AND+{date_q}","drug_total":f"{drug_q}+AND+{date_q}"}; counts={k:openfda_count(v) for k,v in queries.items()}
            if all(v is not None for v in counts.values()):
                a=counts["a_exposed_event"]; exposed_total=counts["exposed_total"]; drug_event_total=counts["drug_event_total"]; drug_total=counts["drug_total"]; b=max(exposed_total-a,0); c=max(drug_event_total-a,0); d=max(drug_total-exposed_total-c,0); aa,bb,cc,dd=a+.5,b+.5,c+.5,d+.5; ror=(aa/bb)/(cc/dd) if bb>0 and cc>0 and dd>0 else None; se=math.sqrt(1/aa+1/bb+1/cc+1/dd); ci=[math.exp(math.log(ror)-1.96*se),math.exp(math.log(ror)+1.96*se)] if ror and ror>0 else [None,None]
            else:a=b=c=d=None;ror=None;ci=[None,None]
            results.append({"pair":pair["name"],"period":period,"counts":counts,"a":a,"b":b,"c":c,"d":d,"ror":ror,"ror95":ci}); time.sleep(.25)
    pd.DataFrame([{**{"pair":r["pair"],"period":r["period"]},**r["counts"],"ror":r["ror"],"ror_lo":r["ror95"][0],"ror_hi":r["ror95"][1]} for r in results]).to_csv(track_dir/"openfda_smoke.csv",index=False)
    replicated=[]
    for pair in {r["pair"] for r in results}:
        rr={r["period"]:r for r in results if r["pair"]==pair}
        if all(p in rr for p in ["early","late"]):
            e,l=rr["early"],rr["late"]
            if e["a"] is not None and l["a"] is not None and e["a"]>=5 and l["a"]>=5 and e["ror"] and l["ror"] and e["ror"]>1 and l["ror"]>1:replicated.append(pair)
    metrics={"pairs":results,"temporally_replicated_positive_control_signals":replicated}
    if replicated:
        status,gate,headline,next_decision="PASS_SIGNAL_DETECTION_SMOKE","At least one prespecified known interaction signal was directionally reproduced in two time windows.","openFDA can support an automated herb–drug signal-discovery pipeline, but not causal risk estimation.","Expand alias normalization, deduplicate cases, use negative controls, and require independent confirmation in a second reporting system."
    else:
        status,gate,headline,next_decision="FAIL_SPARSE_OR_UNSTABLE_SIGNAL_GATE","No prespecified interaction had >=5 exposed-event reports and ROR>1 in both time windows.","The public reporting data are too sparse or inconsistently coded for a robust automated signal in this smoke test.","Improve product normalization or move to a different pharmacovigilance source."
    return TrackResult("5_herb_drug_safety","Can reproducible hidden herb–drug safety signals be detected from public adverse-event reports?",status,gate,headline,metrics,["Spontaneous reports do not establish causality or incidence.","A report may contain many products and reactions; indication and reporting biases remain.","Product-name normalization is incomplete in this smoke test."],next_decision,[str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])


def track_icu_heterogeneity() -> TrackResult:
    track_dir=OUT/"06_icu_heterogeneity"; track_dir.mkdir(exist_ok=True); base="https://physionet.org/files/eicu-crd-demo/2.0.1"; names=["patient.csv.gz","diagnosis.csv.gz","infusiondrug.csv.gz","apachePatientResult.csv.gz"]
    try:
        files={name:download(f"{base}/{name}",RAW/f"eicu_{name}",min_bytes=100) for name in names}; patient=pd.read_csv(files["patient.csv.gz"],compression="gzip",low_memory=False); dx=pd.read_csv(files["diagnosis.csv.gz"],compression="gzip",low_memory=False); inf=pd.read_csv(files["infusiondrug.csv.gz"],compression="gzip",low_memory=False); apr=pd.read_csv(files["apachePatientResult.csv.gz"],compression="gzip",low_memory=False)
    except Exception as e:
        return TrackResult("6_icu_treatment_heterogeneity","Which ICU patients benefit or are harmed by early treatment strategies?","DOWNLOAD_OR_PARSE_FAILED","eICU demo unavailable",str(e),{},["No causal estimate was attempted."],"Retry download.",[])
    dx_text_col=next((c for c in dx.columns if c.lower() in {"diagnosisstring","diagnosis"}),None) or [c for c in dx.columns if "diagnosis" in c.lower()][0]; sepsis_ids=set(dx.loc[dx[dx_text_col].astype(str).str.contains(r"sepsis|septic",case=False,regex=True,na=False),"patientunitstayid"].astype(int)); cohort=patient[patient["patientunitstayid"].isin(sepsis_ids)].copy(); drug_col=next((c for c in inf.columns if "drugname" in c.lower()),None); offset_col=next((c for c in inf.columns if c.lower() in {"infusionoffset","drugstartoffset"}),None)
    if drug_col and offset_col:
        vaso=inf[inf[drug_col].astype(str).str.contains(r"norepinephrine|levophed|noradrenaline",case=False,regex=True,na=False)].copy(); vaso["early"]=pd.to_numeric(vaso[offset_col],errors="coerce").between(0,360); early_ids=set(vaso.loc[vaso["early"],"patientunitstayid"].astype(int))
    else:early_ids=set()
    cohort["early_norepi"]=cohort["patientunitstayid"].astype(int).isin(early_ids).astype(int); cohort["mortality"]=cohort["hospitaldischargestatus"].astype(str).str.lower().str.contains("expired").astype(int); score_col=next((c for c in apr.columns if c.lower() in {"apachescore","acutephysiologyscore"}),None)
    if score_col:cohort=cohort.merge(apr[["patientunitstayid",score_col]].drop_duplicates("patientunitstayid"),on="patientunitstayid",how="left")
    else:cohort["apache_score"]=np.nan;score_col="apache_score"
    cohort["age_num"]=pd.to_numeric(cohort["age"].astype(str).replace({"> 89":"90",">89":"90"}),errors="coerce"); cohort.to_csv(track_dir/"eicu_sepsis_smoke_cohort.csv",index=False); n=len(cohort);n_t=int(cohort["early_norepi"].sum());deaths=int(cohort["mortality"].sum());metrics={"n_sepsis_stays":n,"n_early_norepi":n_t,"n_no_early_norepi":n-n_t,"deaths":deaths}
    if n>=100 and n_t>=20 and (n-n_t)>=20 and deaths>=20:
        import statsmodels.api as sm
        dat=cohort[["early_norepi","mortality","age_num",score_col,"gender"]].copy(); dat=pd.get_dummies(dat,columns=["gender"],drop_first=True,dtype=float).replace([np.inf,-np.inf],np.nan).dropna(); Xp=sm.add_constant(dat.drop(columns=["early_norepi","mortality"]),has_constant="add")
        try:
            psmod=sm.Logit(dat["early_norepi"],Xp).fit(disp=0);ps=np.clip(psmod.predict(Xp),.05,.95);w=np.where(dat["early_norepi"].to_numpy()==1,1/ps,1/(1-ps));sev=pd.to_numeric(dat[score_col],errors="coerce");dat["high_severity"]=(sev>=sev.median()).astype(int);dat["interaction"]=dat["early_norepi"]*dat["high_severity"];Xo=sm.add_constant(dat[["early_norepi","high_severity","interaction","age_num"]],has_constant="add");fit=sm.GLM(dat["mortality"],Xo,family=sm.families.Binomial(),freq_weights=w).fit(cov_type="HC0");metrics["technical_weighted_model"]={k:{"coef":float(fit.params[k]),"or":float(np.exp(fit.params[k])),"p":float(fit.pvalues[k])} for k in fit.params.index};metrics["n_model_complete"]=len(dat);status,gate,headline,next_decision="PASS_TECHNICAL_SMOKE_ONLY","Demo cohort supports end-to-end cohort construction and heterogeneous-effect code.","The pipeline runs, but the demo cannot establish a causal treatment-effect subgroup.","Obtain credentialed full MIMIC-IV/eICU data, emulate a prespecified target trial, and externally validate subgroup effects."
        except Exception as e:metrics["model_error"]=repr(e);status,gate,headline,next_decision="FAIL_MODEL_STABILITY_GATE","Cohort exists but the demo model was unstable.","Technical feasibility is incomplete.","Use the full dataset or simplify exposure definition."
    else:status,gate,headline,next_decision="FAIL_SAMPLE_GATE","Demo cohort lacks the sample/event/exposure support needed even for a technical heterogeneity smoke test.","The open demo is useful for schema testing, not treatment-effect discovery.","Credentialed full data are required."
    return TrackResult("6_icu_treatment_heterogeneity","Which ICU patients benefit or are harmed by early treatment strategies?",status,gate,headline,metrics,["This is observational and confounded by indication.","The eICU demo is a nonrepresentative subset designed for code testing.","A treatment-effect claim requires a target-trial protocol, stronger covariate history, and external validation."],next_decision,[str(p.relative_to(OUT)) for p in track_dir.rglob("*") if p.is_file()])


def make_report(results:list[TrackResult]):
    pd.DataFrame([{"track":r.track,"status":r.status,"gate":r.gate,"headline":r.headline,"next":r.next_decision} for r in results]).to_csv(TABLES/"tournament_summary.csv",index=False);(OUT/"smoke_results.json").write_text(json.dumps([safe_json(asdict(r)) for r in results],ensure_ascii=False,indent=2),encoding="utf-8");rank_map={"PROMISING_INCREMENTAL_SIGNAL":100,"PASS_SMOKE":90,"PASS_SIGNAL_DETECTION_SMOKE":80,"PASS_DATA_GATE":75,"PASS_TECHNICAL_SMOKE_ONLY":40,"PASS_STRUCTURE_FAIL_OUTCOME_DATA":35};ranked=sorted(results,key=lambda r:rank_map.get(r.status,0),reverse=True)
    lines=["# Open-Problem Tournament — Smoke-Test Report","","**Run date:** 2026-08-16  ","**Rule:** no track advances because it sounds important; it advances only if public data support a falsifiable, independently rerunnable test.","","## Tournament table","","| Rank | Track | Status | Gate conclusion |","|---:|---|---|---|"]
    for i,r in enumerate(ranked,1):lines.append(f"| {i} | {r.track} | `{r.status}` | {r.headline.replace('|','/')} |")
    for r in ranked:
        lines += ["",f"## {r.track}","",f"**Question:** {r.question}","",f"**Status:** `{r.status}`","",f"**Gate:** {r.gate}","",f"**Finding:** {r.headline}","","### Key metrics","","```json",json.dumps(safe_json(r.metrics),ensure_ascii=False,indent=2)[:25000],"```","","### Limitations",""]+[f"- {x}" for x in r.limitations]+["",f"**Decision:** {r.next_decision}"]
    (OUT/"SMOKE_TEST_REPORT.md").write_text("\n".join(lines),encoding="utf-8")


def main():
    tracks=[("pattern_biology",track_pattern_biology),("formula_synergy",track_formula_synergy),("acupoint_specificity",track_acupoint_specificity),("pattern_incremental_value",track_pattern_incremental_value),("herb_drug_safety",track_herb_drug_safety),("icu_heterogeneity",track_icu_heterogeneity)];results=[]
    for name,fn in tracks:
        log(f"\n===== {name} =====")
        try:r=fn()
        except Exception as e:traceback.print_exc();r=TrackResult(name,"","UNEXPECTED_ERROR","Unhandled error",repr(e),{"traceback":traceback.format_exc()},[],"Debug before interpretation.",[])
        results.append(r);log(f"{r.track}: {r.status} — {r.headline}");(OUT/f"checkpoint_{name}.json").write_text(json.dumps(safe_json(asdict(r)),ensure_ascii=False,indent=2),encoding="utf-8")
    make_report(results);manifest={}
    for p in sorted(OUT.rglob("*")):
        if p.is_file():manifest[str(p.relative_to(OUT))]={"bytes":p.stat().st_size,"sha256":sha256(p)}
    (OUT/"MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8");archive=OUT.parent/"open_problem_tournament_smoke.zip"
    if archive.exists():archive.unlink()
    with zipfile.ZipFile(archive,"w",zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(OUT.rglob("*")):
            if p.is_file() and p.stat().st_size<50_000_000:zf.write(p,arcname=str(Path(OUT.name)/p.relative_to(OUT)))
    log(f"DONE {archive}")

if __name__=="__main__":main()
