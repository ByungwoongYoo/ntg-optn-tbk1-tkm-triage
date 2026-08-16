#!/usr/bin/env python3
"""Family-cluster-held-out ProteinGym clinical missense stacking and dynamic-gating audit.

Primary scope
-------------
* ProteinGym v1.3 clinical substitutions.
* Five public zero-shot component scores: PoET, TranceptEVE_L, GEMME, EVE, ESM1b.
* Complete-case mutation set: every method is evaluated on exactly the same variants.
* Homology components are created outside this script with MMseqs2 and supplied as an edge list.
* Every directly connected sequence pair is kept in the same outer/inner fold.
* Primary unit is the protein; uncertainty resamples homology clusters.

The dynamic gate predicts, for each variant, which component score is likely to assign more
probability to the true class. It is trained only within outer-training homology clusters.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.special import softmax
from scipy.stats import wilcoxon
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20260817
ALLOWED = ["PoET", "TranceptEVE_L", "GEMME", "EVE", "ESM1b"]
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}
HYDRO = dict(zip(AA, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                      1.9, -3.5, -1.6, -0.8, -0.7, -0.9, -0.8, 4.2, -1.3, -1.3]))
VOLUME = dict(zip(AA, [88.6,108.5,111.1,138.4,189.9,60.1,153.2,166.7,168.6,166.7,
                       162.9,114.1,112.7,143.9,173.4,89.0,116.1,140.0,227.8,193.6]))
POLAR = dict(zip(AA, [8.1,5.5,13.0,12.3,5.2,9.0,10.4,5.2,11.3,4.9,5.7,11.6,8.0,10.5,10.5,9.2,8.6,5.9,5.4,6.2]))
CHARGE = {a: 0.0 for a in AA}
for a in "DE": CHARGE[a] = -1.0
for a in "KR": CHARGE[a] = 1.0
CHARGE["H"] = 0.1
AROMATIC = {a: float(a in "FWY") for a in AA}
SPECIAL = {a: float(a in "GPC") for a in AA}


def label01(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().map({
        "pathogenic": 1, "likely pathogenic": 1,
        "benign": 0, "likely benign": 0,
    })


def rank01(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").rank(method="average", pct=True)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def safe_auc(y: Iterable[float], p: Iterable[float]) -> float:
    y = np.asarray(y, dtype=float); p = np.asarray(p, dtype=float)
    m = np.isfinite(y) & np.isfinite(p)
    if m.sum() < 4 or np.unique(y[m]).size != 2:
        return float("nan")
    return float(roc_auc_score(y[m].astype(int), p[m]))


def safe_ap(y: Iterable[float], p: Iterable[float]) -> float:
    y = np.asarray(y, dtype=float); p = np.asarray(p, dtype=float)
    m = np.isfinite(y) & np.isfinite(p)
    if m.sum() < 4 or np.unique(y[m]).size != 2:
        return float("nan")
    return float(average_precision_score(y[m].astype(int), p[m]))


def clip_prob(p: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)


class UnionFind:
    def __init__(self, items: Iterable[str]):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x: str) -> str:
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb: return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]: self.rank[ra] += 1


def parse_cluster_edges(headers: list[str], edge_path: Path) -> tuple[dict[str, str], dict[str, Any]]:
    uf = UnionFind(headers)
    edges = 0
    max_fident = 0.0
    if edge_path.exists() and edge_path.stat().st_size > 0:
        ed = pd.read_csv(edge_path, sep="\t", header=None,
                         names=["query", "target", "fident", "qcov", "tcov", "evalue", "bits"])
        for row in ed.itertuples(index=False):
            q, t = str(row.query), str(row.target)
            if q in uf.parent and t in uf.parent and q != t:
                uf.union(q, t); edges += 1
                try: max_fident = max(max_fident, float(row.fident))
                except Exception: pass
    components: dict[str, list[str]] = {}
    for h in headers:
        components.setdefault(uf.find(h), []).append(h)
    mapping = {}
    for members in components.values():
        cid = "C_" + sha256_text("|".join(sorted(members)))[:16]
        for h in members: mapping[h] = cid
    sizes = sorted((len(v) for v in components.values()), reverse=True)
    diag = {
        "edge_rows_used": edges,
        "components": len(components),
        "singletons": int(sum(s == 1 for s in sizes)),
        "largest_component": int(sizes[0] if sizes else 0),
        "components_ge2": int(sum(s >= 2 for s in sizes)),
        "max_reported_fident": max_fident,
    }
    return mapping, diag


def assign_cluster_folds(df: pd.DataFrame, cluster_col: str, n_folds: int) -> dict[str, int]:
    stats = df[["protein_file", cluster_col]].drop_duplicates().groupby(cluster_col).agg(
        proteins=("protein_file", "nunique")
    )
    var_counts = df.groupby(cluster_col).size().rename("variants")
    stats = stats.join(var_counts)
    stats["tie"] = [int(sha256_text(str(x))[:12], 16) for x in stats.index]
    stats = stats.sort_values(["proteins", "variants", "tie"], ascending=[False, False, True])
    load_p = np.zeros(n_folds, dtype=int)
    load_v = np.zeros(n_folds, dtype=int)
    out = {}
    for cid, row in stats.iterrows():
        candidates = sorted(range(n_folds), key=lambda f: (load_p[f], load_v[f], f))
        f = candidates[0]
        out[str(cid)] = f
        load_p[f] += int(row.proteins); load_v[f] += int(row.variants)
    return out


def parse_mutant(mut: str) -> tuple[str, int, str] | None:
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", str(mut).strip().upper())
    if not m: return None
    return m.group(1), int(m.group(2)), m.group(3)


def load_common_scores(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    frames = []
    inventory = []
    failures = []
    original_labeled = 0
    with zipfile.ZipFile(path) as z:
        members = sorted(n for n in z.namelist() if n.lower().endswith(".csv"))
        for i, member in enumerate(members, 1):
            try:
                with z.open(member) as f:
                    d = pd.read_csv(f, low_memory=False)
                required = {"mutant", "DMS_bin_score", "protein_sequence"}.union(ALLOWED)
                if not required.issubset(d.columns):
                    failures.append({"file": member, "reason": "missing required or component score fields"})
                    continue
                seqs = d["protein_sequence"].dropna().astype(str).str.strip().str.upper().unique()
                if len(seqs) != 1:
                    failures.append({"file": member, "reason": f"sequence_count={len(seqs)}"}); continue
                seq = seqs[0]
                protein_file = Path(member).name
                header = protein_file.removesuffix(".csv")
                d["label"] = label01(d["DMS_bin_score"])
                d = d[d["label"].isin([0,1])].drop_duplicates("mutant").copy()
                original_labeled += len(d)
                missing_by_model = {m: int(pd.to_numeric(d[m], errors="coerce").isna().sum()) for m in ALLOWED}
                for m in ALLOWED: d[m] = pd.to_numeric(d[m], errors="coerce")
                before = len(d)
                d = d.dropna(subset=ALLOWED).copy()
                if len(d) < 4 or d["label"].nunique() != 2:
                    failures.append({"file": member, "reason": "insufficient complete-case two-class variants"}); continue
                for m in ALLOWED: d[m] = rank01(d[m])
                d["protein_file"] = protein_file
                d["header"] = header
                d["sequence"] = seq
                d["sequence_hash"] = sha256_text(seq)
                d["length"] = len(seq)
                parsed = d["mutant"].map(parse_mutant)
                d["wt"] = parsed.map(lambda x: x[0] if x else "X")
                d["position"] = parsed.map(lambda x: x[1] if x else np.nan)
                d["mut"] = parsed.map(lambda x: x[2] if x else "X")
                d["sequence_match"] = [
                    bool(np.isfinite(pos) and 1 <= int(pos) <= len(seq) and seq[int(pos)-1] == wt)
                    for pos, wt in zip(d["position"], d["wt"])
                ]
                frames.append(d[["protein_file","header","sequence","sequence_hash","length","mutant","wt","position","mut","sequence_match","label"] + ALLOWED])
                inventory.append({
                    "protein_file": protein_file, "header": header, "sequence_hash": sha256_text(seq),
                    "length": len(seq), "labeled_before_complete_case": before,
                    "complete_case_variants": len(d), "missing_by_model": missing_by_model,
                })
            except Exception as exc:
                failures.append({"file": member, "reason": repr(exc)})
            if i % 500 == 0:
                print(f"load {i}/{len(members)} accepted={len(frames)}", flush=True)
    if not frames: raise RuntimeError("No complete-case clinical proteins loaded")
    data = pd.concat(frames, ignore_index=True)
    inv = pd.DataFrame(inventory)
    diagnostics = {
        "members_scanned": len(members), "accepted_proteins": int(data["protein_file"].nunique()),
        "complete_case_variants": int(len(data)), "original_labeled_variants_in_accepted_or_failed_files": int(original_labeled),
        "variant_retention_fraction": float(len(data) / original_labeled) if original_labeled else None,
        "sequence_mismatch_rows": int((~data["sequence_match"]).sum()),
        "failures": failures,
    }
    return data, inv, diagnostics


def protein_class_weights(df: pd.DataFrame) -> np.ndarray:
    counts = df.groupby(["protein_file","label"]).size().to_dict()
    w = np.array([0.5 / counts[(p, int(y))] for p, y in zip(df["protein_file"], df["label"])], dtype=float)
    return w / np.mean(w)


def orientation_performance(df: pd.DataFrame) -> tuple[dict[str,int], dict[str,float]]:
    signs = {}; perf = {}
    for m in ALLOWED:
        vals = [safe_auc(g["label"], g[m]) for _, g in df.groupby("protein_file", sort=False)]
        mean = float(np.nanmean(vals))
        signs[m] = 1 if not np.isfinite(mean) or mean >= 0.5 else -1
        perf[m] = mean if signs[m] > 0 else 1.0 - mean
    return signs, perf


def oriented_scores(df: pd.DataFrame, signs: dict[str,int]) -> np.ndarray:
    x = df[ALLOWED].to_numpy(float, copy=True)
    for j,m in enumerate(ALLOWED):
        if signs[m] < 0: x[:,j] = 1.0 - x[:,j]
    return x


def aa_one_hot(values: Iterable[str]) -> np.ndarray:
    vals = list(values); out = np.zeros((len(vals), len(AA)), dtype=float)
    for i,a in enumerate(vals):
        if a in AA_INDEX: out[i,AA_INDEX[a]] = 1.0
    return out


def local_context(df: pd.DataFrame, offsets=(-2,-1,1,2)) -> np.ndarray:
    out = np.zeros((len(df), len(offsets)*(len(AA)+1)), dtype=float)
    for i,(seq,pos) in enumerate(zip(df["sequence"],df["position"])):
        if not np.isfinite(pos): continue
        pos0 = int(pos)-1
        for oi,off in enumerate(offsets):
            j = pos0+off
            aa = seq[j] if 0 <= j < len(seq) else "X"
            idx = AA_INDEX.get(aa, len(AA))
            out[i,oi*(len(AA)+1)+idx] = 1.0
    return out


def feature_matrix(df: pd.DataFrame, signs: dict[str,int]) -> np.ndarray:
    s = oriented_scores(df, signs)
    summaries = np.column_stack([
        np.mean(s,axis=1), np.median(s,axis=1), np.std(s,axis=1),
        np.min(s,axis=1), np.max(s,axis=1), np.ptp(s,axis=1),
    ])
    pairs = []
    for i in range(s.shape[1]):
        for j in range(i+1,s.shape[1]): pairs.append(np.abs(s[:,i]-s[:,j]))
    pairmat = np.column_stack(pairs)
    wt = df["wt"].astype(str).tolist(); mut = df["mut"].astype(str).tolist()
    wt_oh = aa_one_hot(wt); mut_oh = aa_one_hot(mut)
    pos = pd.to_numeric(df["position"], errors="coerce").to_numpy(float)
    length = df["length"].to_numpy(float)
    pos_frac = np.divide(pos, length, out=np.zeros_like(pos), where=length>0)
    scalars = np.column_stack([pos_frac, np.log1p(length), df["sequence_match"].astype(float)])
    phys = np.column_stack([
        [HYDRO.get(b,0)-HYDRO.get(a,0) for a,b in zip(wt,mut)],
        [abs(HYDRO.get(b,0)-HYDRO.get(a,0)) for a,b in zip(wt,mut)],
        [VOLUME.get(b,0)-VOLUME.get(a,0) for a,b in zip(wt,mut)],
        [abs(VOLUME.get(b,0)-VOLUME.get(a,0)) for a,b in zip(wt,mut)],
        [POLAR.get(b,0)-POLAR.get(a,0) for a,b in zip(wt,mut)],
        [CHARGE.get(b,0)-CHARGE.get(a,0) for a,b in zip(wt,mut)],
        [AROMATIC.get(b,0)-AROMATIC.get(a,0) for a,b in zip(wt,mut)],
        [SPECIAL.get(b,0)-SPECIAL.get(a,0) for a,b in zip(wt,mut)],
    ])
    return np.column_stack([s,summaries,pairmat,scalars,phys,wt_oh,mut_oh,local_context(df)])


@dataclass(frozen=True)
class FixedSpec:
    name: str; kind: str; k: int; power: float=1.0

FIXED_SPECS = [
    FixedSpec("mean3","mean",3), FixedSpec("mean5","mean",5),
    FixedSpec("weighted3_p1","weighted",3,1), FixedSpec("weighted3_p2","weighted",3,2),
    FixedSpec("weighted5_p1","weighted",5,1), FixedSpec("weighted5_p2","weighted",5,2),
]

@dataclass(frozen=True)
class GateSpec:
    name: str; learner: str; temperature: float

GATE_SPECS = [
    GateSpec("ridge_gate_t010","ridge",0.10),
    GateSpec("hgb_gate_t010","hgb",0.10),
]

@dataclass(frozen=True)
class DirectSpec:
    name: str; C: float

DIRECT_SPECS = [DirectSpec("direct_logit_c001",0.01), DirectSpec("direct_logit_c01",0.1), DirectSpec("direct_logit_c1",1.0)]


class FixedModel:
    def __init__(self, spec: FixedSpec, signs: dict[str,int], perf: dict[str,float]):
        self.spec=spec; self.signs=signs; self.perf=perf
        self.models = sorted(perf, key=lambda m:(-np.nan_to_num(perf[m],nan=-9),m))[:spec.k]
    def predict(self, df: pd.DataFrame) -> np.ndarray:
        all_s = oriented_scores(df,self.signs)
        idx=[ALLOWED.index(m) for m in self.models]; x=all_s[:,idx]
        if self.spec.kind=="mean": return np.mean(x,axis=1)
        w=np.array([max(self.perf[m]-0.5,0.001)**self.spec.power for m in self.models]); w/=w.sum()
        return x@w


class GateModel:
    def __init__(self, spec: GateSpec, signs: dict[str,int], regressors: list[Any]):
        self.spec=spec; self.signs=signs; self.regressors=regressors
    def predict(self, df: pd.DataFrame, return_weights: bool=False):
        X=feature_matrix(df,self.signs); s=oriented_scores(df,self.signs)
        util=np.column_stack([r.predict(X) for r in self.regressors])
        w=softmax((util-util.mean(axis=1,keepdims=True))/self.spec.temperature,axis=1)
        pred=np.sum(w*s,axis=1)
        return (pred,w) if return_weights else pred


class DirectModel:
    def __init__(self,spec:DirectSpec,signs:dict[str,int],model:Any): self.spec=spec;self.signs=signs;self.model=model
    def predict(self,df:pd.DataFrame)->np.ndarray: return self.model.predict_proba(feature_matrix(df,self.signs))[:,1]


def fit_fixed(df: pd.DataFrame, spec: FixedSpec) -> FixedModel:
    signs,perf=orientation_performance(df); return FixedModel(spec,signs,perf)


def fit_gate(df: pd.DataFrame, spec: GateSpec) -> GateModel:
    signs,_=orientation_performance(df); X=feature_matrix(df,signs); s=oriented_scores(df,signs)
    y=df["label"].to_numpy(float); sw=protein_class_weights(df); regs=[]
    for j in range(len(ALLOWED)):
        target=1.0-(s[:,j]-y)**2
        if spec.learner=="ridge":
            model=make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),Ridge(alpha=20.0))
            fit_params={"ridge__sample_weight":sw}
        else:
            model=make_pipeline(SimpleImputer(strategy="median"),HistGradientBoostingRegressor(
                max_iter=120,learning_rate=0.05,max_leaf_nodes=15,min_samples_leaf=80,l2_regularization=3.0,random_state=SEED+j))
            fit_params={"histgradientboostingregressor__sample_weight":sw}
        model.fit(X,target,**fit_params)
        regs.append(model)
    return GateModel(spec,signs,regs)


def fit_direct(df:pd.DataFrame,spec:DirectSpec)->DirectModel:
    signs,_=orientation_performance(df); X=feature_matrix(df,signs); y=df["label"].to_numpy(int); sw=protein_class_weights(df)
    model=make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),LogisticRegression(C=spec.C,max_iter=3000,class_weight=None,solver="lbfgs"))
    model.fit(X,y,logisticregression__sample_weight=sw)
    return DirectModel(spec,signs,model)


def fit_best_individual(df:pd.DataFrame):
    signs,perf=orientation_performance(df); best=max(perf,key=lambda m:np.nan_to_num(perf[m],nan=-9)); return best,signs


def predict_best_individual(df:pd.DataFrame,best:str,signs:dict[str,int])->np.ndarray:
    p=df[best].to_numpy(float,copy=True)
    return p if signs[best]>0 else 1-p


def per_protein_auc(df:pd.DataFrame,p:np.ndarray,name:str)->pd.DataFrame:
    t=df[["protein_file","cluster_id","label"]].copy();t[name]=p;rows=[]
    for protein,g in t.groupby("protein_file",sort=False):
        rows.append({"protein_file":protein,"cluster_id":g["cluster_id"].iloc[0],"n_variants":len(g),
                     "n_pathogenic":int(g["label"].sum()),name:safe_auc(g["label"],g[name])})
    return pd.DataFrame(rows)


def mean_protein_auc(df:pd.DataFrame,p:np.ndarray)->float:
    return float(per_protein_auc(df,p,"auc")["auc"].mean())


def choose_spec(outer_train:pd.DataFrame,kind:str,n_inner:int=3):
    folds=assign_cluster_folds(outer_train,"cluster_id",n_inner)
    fseries=outer_train["cluster_id"].map(folds)
    specs=FIXED_SPECS if kind=="fixed" else GATE_SPECS if kind=="gate" else DIRECT_SPECS
    rows=[]
    for spec in specs:
        vals=[]
        for f in range(n_inner):
            tr=outer_train[fseries!=f].copy();va=outer_train[fseries==f].copy()
            if kind=="fixed": model=fit_fixed(tr,spec)
            elif kind=="gate": model=fit_gate(tr,spec)
            else: model=fit_direct(tr,spec)
            vals.append(mean_protein_auc(va,model.predict(va)))
        rows.append({"name":spec.name,"mean_inner_auc":float(np.mean(vals)),"fold_auc":vals})
        print(f"inner {kind} {spec.name}: {np.mean(vals):.6f}",flush=True)
    tab=pd.DataFrame(rows).sort_values(["mean_inner_auc","name"],ascending=[False,True])
    chosen_name=str(tab.iloc[0]["name"]);chosen=next(s for s in specs if s.name==chosen_name)
    return chosen,tab,folds


def crossfit_method(df:pd.DataFrame,kind:str,spec:Any,folds:dict[str,int],n_inner:int=3)->np.ndarray:
    out=np.full(len(df),np.nan);fs=df["cluster_id"].map(folds).to_numpy()
    for f in range(n_inner):
        tr=df[fs!=f].copy();va=df[fs==f].copy()
        model=fit_fixed(tr,spec) if kind=="fixed" else fit_gate(tr,spec) if kind=="gate" else fit_direct(tr,spec)
        out[np.flatnonzero(fs==f)]=model.predict(va)
    return out


def fit_platt(raw:np.ndarray,y:np.ndarray,weights:np.ndarray):
    raw=np.asarray(raw,float);y=np.asarray(y,int);m=np.isfinite(raw)&np.isfinite(y)
    if m.sum()<20 or np.unique(y[m]).size<2: return None
    model=LogisticRegression(C=1e4,solver="lbfgs",max_iter=2000)
    model.fit(raw[m].reshape(-1,1),y[m],sample_weight=weights[m]);return model


def apply_platt(model,raw:np.ndarray)->np.ndarray:
    if model is None:return clip_prob(raw)
    return clip_prob(model.predict_proba(np.asarray(raw).reshape(-1,1))[:,1])


def compute_variant_metrics(y,p)->dict[str,float]:
    y=np.asarray(y,int);p=clip_prob(p)
    return {"auc":safe_auc(y,p),"average_precision":safe_ap(y,p),"brier":float(brier_score_loss(y,p)),"log_loss":float(log_loss(y,p,labels=[0,1]))}


def cluster_bootstrap(diff_tab:pd.DataFrame,diff_col:str,n:int=20000)->list[float|None]:
    grouped={cid:g[diff_col].dropna().to_numpy(float) for cid,g in diff_tab.groupby("cluster_id")}
    ids=list(grouped)
    if len(ids)<3:return [None,None]
    rng=np.random.default_rng(SEED);vals=np.empty(n)
    for b in range(n):
        sampled=rng.choice(ids,size=len(ids),replace=True);parts=[grouped[x] for x in sampled]
        vals[b]=np.mean(np.concatenate(parts))
    return [float(np.quantile(vals,.025)),float(np.quantile(vals,.975))]


def paired_summary(tab:pd.DataFrame,a:str,b:str)->dict[str,Any]:
    t=tab[["cluster_id",a,b]].dropna().copy();t["diff"]=t[a]-t[b];d=t["diff"].to_numpy(float)
    cluster_means=t.groupby("cluster_id")["diff"].mean().to_numpy(float)
    try:p=float(wilcoxon(cluster_means,alternative="greater").pvalue) if len(cluster_means)>=5 and np.any(cluster_means!=0) else None
    except Exception:p=None
    ci=cluster_bootstrap(t,"diff")
    return {"n_proteins":int(len(d)),"n_clusters":int(t["cluster_id"].nunique()),
            "mean_a":float(t[a].mean()),"mean_b":float(t[b].mean()),"mean_difference":float(np.mean(d)),
            "cluster_bootstrap_95ci":ci,"fraction_proteins_improved":float(np.mean(d>1e-12)),
            "fraction_tied":float(np.mean(np.abs(d)<=1e-12)),"cluster_mean_wilcoxon_one_sided_p":p,
            "confirmed":bool(ci[0] is not None and ci[0]>0)}


def calibration_slope_intercept(y:np.ndarray,p:np.ndarray,w:np.ndarray)->dict[str,float|None]:
    p=clip_prob(p);x=np.log(p/(1-p)).reshape(-1,1)
    try:
        m=LogisticRegression(C=1e4,solver="lbfgs",max_iter=2000);m.fit(x,y,sample_weight=w)
        return {"intercept":float(m.intercept_[0]),"slope":float(m.coef_[0,0])}
    except Exception:return {"intercept":None,"slope":None}


def run_scheme(data:pd.DataFrame,cluster_map:dict[str,str],scheme:str,out:Path)->dict[str,Any]:
    d=data.copy();d["cluster_id"]=d["header"].map(cluster_map).fillna(d["header"].map(lambda x:"S_"+sha256_text(str(x))[:16]))
    outer_map=assign_cluster_folds(d,"cluster_id",5);d["outer_fold"]=d["cluster_id"].map(outer_map)
    variant_parts=[];protein_parts=[];fold_records=[];tuning_parts=[];weight_parts=[]
    for fold in range(5):
        tr=d[d["outer_fold"]!=fold].copy();te=d[d["outer_fold"]==fold].copy()
        fixed_spec,fixed_tune,inner_folds=choose_spec(tr,"fixed",3)
        gate_spec,gate_tune,_=choose_spec(tr,"gate",3)
        direct_spec,direct_tune,_=choose_spec(tr,"direct",3)
        tuning_parts.append(pd.concat([
            fixed_tune.assign(family="fixed",outer_fold=fold),gate_tune.assign(family="gate",outer_fold=fold),direct_tune.assign(family="direct",outer_fold=fold)
        ],ignore_index=True))
        raw_fixed_oof=crossfit_method(tr,"fixed",fixed_spec,inner_folds,3)
        raw_gate_oof=crossfit_method(tr,"gate",gate_spec,inner_folds,3)
        raw_direct_oof=crossfit_method(tr,"direct",direct_spec,inner_folds,3)
        raw_base_oof=np.full(len(tr),np.nan);fs=tr["cluster_id"].map(inner_folds).to_numpy()
        for f in range(3):
            itr=tr[fs!=f].copy();iva=tr[fs==f].copy();best_i,signs_i=fit_best_individual(itr)
            raw_base_oof[np.flatnonzero(fs==f)]=predict_best_individual(iva,best_i,signs_i)
        sw_tr=protein_class_weights(tr);y_tr=tr["label"].to_numpy(int)
        calibrators={
            "best_individual":fit_platt(raw_base_oof,y_tr,sw_tr),"fixed":fit_platt(raw_fixed_oof,y_tr,sw_tr),
            "dynamic_gate":fit_platt(raw_gate_oof,y_tr,sw_tr),"direct_stack":fit_platt(raw_direct_oof,y_tr,sw_tr),
        }
        fixed=fit_fixed(tr,fixed_spec);gate=fit_gate(tr,gate_spec);direct=fit_direct(tr,direct_spec);best,signs=fit_best_individual(tr)
        raw_base=predict_best_individual(te,best,signs);raw_fixed=fixed.predict(te);raw_gate,w_gate=gate.predict(te,return_weights=True);raw_direct=direct.predict(te)
        raws={"best_individual":raw_base,"fixed":raw_fixed,"dynamic_gate":raw_gate,"direct_stack":raw_direct}
        preds={k:apply_platt(calibrators[k],v) for k,v in raws.items()}
        vp=te[["protein_file","header","cluster_id","outer_fold","mutant","label"]].copy()
        for k,v in raws.items():vp[f"{k}_raw"]=v;vp[f"{k}_prob"]=preds[k]
        variant_parts.append(vp)
        for j,m in enumerate(ALLOWED):
            weight_parts.append(pd.DataFrame({"protein_file":te["protein_file"],"mutant":te["mutant"],"outer_fold":fold,"model":m,"weight":w_gate[:,j]}))
        ptab=None
        for k,v in raws.items():
            pt=per_protein_auc(te,v,k)
            ptab=pt if ptab is None else ptab.merge(pt[["protein_file",k]],on="protein_file",how="inner")
        ptab["outer_fold"]=fold;protein_parts.append(ptab)
        fold_records.append({
            "fold":fold,"train_clusters":int(tr["cluster_id"].nunique()),"test_clusters":int(te["cluster_id"].nunique()),
            "train_proteins":int(tr["protein_file"].nunique()),"test_proteins":int(te["protein_file"].nunique()),
            "fixed_spec":fixed_spec.name,"gate_spec":gate_spec.name,"direct_spec":direct_spec.name,
            "best_individual":best,"fixed_models":fixed.models,
            "mean_auc_best":mean_protein_auc(te,raw_base),"mean_auc_fixed":mean_protein_auc(te,raw_fixed),
            "mean_auc_gate":mean_protein_auc(te,raw_gate),"mean_auc_direct":mean_protein_auc(te,raw_direct),
        })
        print(json.dumps(fold_records[-1],indent=2),flush=True)
    variants=pd.concat(variant_parts,ignore_index=True);proteins=pd.concat(protein_parts,ignore_index=True)
    tuning=pd.concat(tuning_parts,ignore_index=True);weights=pd.concat(weight_parts,ignore_index=True)
    y=variants["label"].to_numpy(int);sw=protein_class_weights(variants)
    variant_metrics={k:compute_variant_metrics(y,variants[f"{k}_prob"].to_numpy(float)) for k in ["best_individual","fixed","dynamic_gate","direct_stack"]}
    calibration={k:calibration_slope_intercept(y,variants[f"{k}_prob"].to_numpy(float),sw) for k in ["best_individual","fixed","dynamic_gate","direct_stack"]}
    comparisons={
        "fixed_vs_best":paired_summary(proteins,"fixed","best_individual"),
        "dynamic_vs_fixed":paired_summary(proteins,"dynamic_gate","fixed"),
        "dynamic_vs_best":paired_summary(proteins,"dynamic_gate","best_individual"),
        "direct_vs_fixed":paired_summary(proteins,"direct_stack","fixed"),
    }
    cluster_sizes=d[["protein_file","cluster_id"]].drop_duplicates().groupby("cluster_id").size()
    result={
        "scheme":scheme,"proteins":int(d["protein_file"].nunique()),"variants":int(len(d)),"clusters":int(d["cluster_id"].nunique()),
        "largest_cluster":int(cluster_sizes.max()),"clusters_ge2":int((cluster_sizes>=2).sum()),
        "folds":fold_records,"comparisons":comparisons,"variant_metrics":variant_metrics,"calibration":calibration,
        "primary_success":bool(comparisons["fixed_vs_best"]["confirmed"] and comparisons["dynamic_vs_fixed"]["mean_difference"]>=0.005 and comparisons["dynamic_vs_fixed"]["confirmed"] and variant_metrics["dynamic_gate"]["brier"]<=variant_metrics["fixed"]["brier"]+1e-6),
    }
    out.mkdir(parents=True,exist_ok=True)
    variants.to_csv(out/"family_oof_variant_predictions.csv.gz",index=False,compression="gzip")
    proteins.to_csv(out/"family_oof_protein_auc.csv",index=False)
    tuning.to_csv(out/"inner_tuning_results.csv",index=False)
    weights.to_csv(out/"dynamic_gate_weights.csv.gz",index=False,compression="gzip")
    (out/"result.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    ci=comparisons["dynamic_vs_fixed"]["cluster_bootstrap_95ci"]
    (out/"REPORT.md").write_text(f"""# ProteinGym family-held-out clinical missense audit: {scheme}

- Complete-case proteins: {result['proteins']}
- Complete-case variants: {result['variants']:,}
- Homology components: {result['clusters']} (largest {result['largest_cluster']})
- Fixed ensemble vs training-selected best individual: {comparisons['fixed_vs_best']['mean_difference']:+.4f}; cluster-bootstrap 95% CI {comparisons['fixed_vs_best']['cluster_bootstrap_95ci']}
- Dynamic gate vs fixed ensemble: {comparisons['dynamic_vs_fixed']['mean_difference']:+.4f}; cluster-bootstrap 95% CI {ci}
- Dynamic gate primary success: {result['primary_success']}

All four methods were evaluated on the identical complete-case mutations. Outer and inner splits kept every MMseqs-connected component intact. The result is a retrospective benchmark test, not a prospective clinical-variant solution.
""",encoding="utf-8")
    return result


def main()->None:
    ap=argparse.ArgumentParser();ap.add_argument("--clinical-score-zip",required=True);ap.add_argument("--sequence-inventory",required=True)
    ap.add_argument("--cluster-edge",action="append",required=True,help="NAME=path to MMseqs edge TSV");ap.add_argument("--out-dir",required=True)
    args=ap.parse_args();out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    data,inventory,load_diag=load_common_scores(Path(args.clinical_score_zip))
    seqinv=pd.read_csv(args.sequence_inventory)
    all_results={"load_diagnostics":load_diag,"schemes":{},"grand_problem_fully_solved":False}
    inventory.to_json(out/"complete_case_inventory.json",orient="records",indent=2)
    for item in args.cluster_edge:
        if "=" not in item:raise ValueError("--cluster-edge must be NAME=path")
        name,path=item.split("=",1);mapping,diag=parse_cluster_edges(seqinv["header"].astype(str).tolist(),Path(path))
        sdir=out/name;result=run_scheme(data,mapping,name,sdir);result["cluster_diagnostics"]=diag
        (sdir/"result.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
        all_results["schemes"][name]=result
    primary=all_results["schemes"].get("global30_cov80",next(iter(all_results["schemes"].values())))
    sensitivity=[v for k,v in all_results["schemes"].items() if k!="global30_cov80"]
    all_results["family_heldout_fixed_signal_confirmed"]=bool(primary["comparisons"]["fixed_vs_best"]["confirmed"])
    all_results["dynamic_gate_advance_confirmed"]=bool(primary["primary_success"] and all(v["comparisons"]["dynamic_vs_fixed"]["mean_difference"]>0 for v in sensitivity))
    all_results["open_problem_fully_solved"]=False
    (out/"FINAL_RESULT.json").write_text(json.dumps(all_results,ensure_ascii=False,indent=2),encoding="utf-8")
    (out/"FINAL_REPORT.md").write_text(f"""# Family-held-out missense experiment v12

## Decisions

- Family-held-out fixed ensemble signal confirmed: **{all_results['family_heldout_fixed_signal_confirmed']}**
- Variant-conditioned dynamic-gating advance confirmed: **{all_results['dynamic_gate_advance_confirmed']}**
- All-human-missense open problem fully solved: **False**

The first statement tests whether the original clinical stacking signal survives direct homology-component separation and a common mutation set. The second requires a dynamic-gating gain of at least 0.005 with a cluster-bootstrap confidence interval above zero, no Brier-score deterioration, and a positive direction in every sensitivity clustering.
""",encoding="utf-8")
    print(json.dumps(all_results,ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()
