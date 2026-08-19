#!/usr/bin/env python3
"""Population-adjusted, untouched-validation tests for frozen unitig candidates."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import norm

ACC_RE = re.compile(r"GC[AF]_\d+(?:\.\d+)?")


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--selection", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--rtab", required=True)
    p.add_argument("--distance", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args()


def normalize_acc(v: str) -> str:
    m = ACC_RE.search(str(v))
    return m.group(0) if m else Path(str(v)).name.removesuffix(".fna")


def normalize_seq(v: str) -> str:
    s = re.sub(r"[^ACGTNacgtn]", "", str(v)).upper()
    return s


def load_rtab(path: str, sample_ids: set[str], selected_sequences: set[str]) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", index_col=0)
    x.index = x.index.astype(str)
    x.columns = x.columns.astype(str)
    col_overlap = len({normalize_acc(c) for c in x.columns} & sample_ids)
    idx_overlap = len({normalize_acc(i) for i in x.index} & sample_ids)
    if idx_overlap > col_overlap:
        x = x.T
    x.columns = [normalize_acc(c) for c in x.columns]
    x.index = [normalize_seq(i) for i in x.index]
    x = x.apply(pd.to_numeric, errors="coerce").fillna(0)
    x = (x > 0).astype(np.int8)
    if x.index.duplicated().any():
        x = x.groupby(level=0).max()
    # Some versions may preserve FASTA IDs rather than sequences. If so, rows are remapped later by order.
    return x


def load_distance(path: str, ids: list[str]) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", index_col=0)
    x.index = [normalize_acc(v) for v in x.index]
    x.columns = [normalize_acc(v) for v in x.columns]
    if x.index.duplicated().any() or pd.Index(x.columns).duplicated().any():
        raise RuntimeError("Duplicate sample IDs in Mash distance matrix")
    missing = sorted(set(ids) - set(x.index))
    if missing:
        raise RuntimeError(f"Mash distance matrix missing {len(missing)} IDs: {missing[:10]}")
    return x.loc[ids, ids].astype(float)


def pcoa(distance: np.ndarray, max_dim: int = 30) -> np.ndarray:
    d2 = distance ** 2
    b = -0.5 * (d2 - d2.mean(1, keepdims=True) - d2.mean(0, keepdims=True) + d2.mean())
    vals, vecs = np.linalg.eigh(b)
    o = np.argsort(vals)[::-1]
    vals, vecs = vals[o], vecs[:, o]
    keep = vals > max(1e-12, vals[0] * 1e-10)
    vals, vecs = vals[keep][:max_dim], vecs[:, keep][:, :max_dim]
    return vecs * np.sqrt(vals)


def fit_null(y: np.ndarray, cov: np.ndarray):
    X = np.column_stack([np.ones(len(y)), cov])
    beta = np.zeros(X.shape[1])
    for _ in range(100):
        eta = np.clip(X @ beta, -30, 30)
        mu = 1 / (1 + np.exp(-eta))
        w = np.clip(mu * (1 - mu), 1e-8, None)
        z = eta + (y - mu) / w
        lhs = X.T @ (w[:, None] * X) + 1e-8 * np.eye(X.shape[1])
        new = np.linalg.solve(lhs, X.T @ (w * z))
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    eta = np.clip(X @ beta, -30, 30)
    mu = 1 / (1 + np.exp(-eta))
    w = np.clip(mu * (1 - mu), 1e-8, None)
    inv = np.linalg.pinv(X.T @ (w[:, None] * X))
    return X, mu, w, inv


def score(y: np.ndarray, cov: np.ndarray, G: np.ndarray) -> pd.DataFrame:
    X, mu, w, inv = fit_null(y, cov)
    gt = G - X @ (inv @ (X.T @ (w[:, None] * G)))
    u = gt.T @ (y - mu)
    v = np.sum(gt * gt * w[:, None], axis=0)
    valid = v > 1e-12
    z = np.full(G.shape[1], np.nan)
    beta = np.full(G.shape[1], np.nan)
    p = np.full(G.shape[1], np.nan)
    z[valid] = u[valid] / np.sqrt(v[valid])
    beta[valid] = u[valid] / v[valid]
    p[valid] = 2 * norm.sf(np.abs(z[valid]))
    return pd.DataFrame({"beta": beta, "z": z, "p": p})


def bh(s: pd.Series) -> pd.Series:
    p = pd.to_numeric(s, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if len(vals):
        order = np.argsort(vals)
        ranked = vals[order]
        q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out[np.flatnonzero(ok)[order]] = np.clip(q, 0, 1)
    return pd.Series(out, index=s.index)


def orci(a: int, b: int, c: int, d: int):
    aa, bb, cc, dd = map(float, (a, b, c, d))
    if min(aa, bb, cc, dd) == 0:
        aa += .5; bb += .5; cc += .5; dd += .5
    lor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    return math.exp(lor), math.exp(lor - 1.96*se), math.exp(lor + 1.96*se)


def counts(g: np.ndarray, y: np.ndarray):
    a = int(((g == 1) & (y == 1)).sum())
    b = int(((g == 1) & (y == 0)).sum())
    c = int(((g == 0) & (y == 1)).sum())
    d = int(((g == 0) & (y == 0)).sum())
    return a, b, c, d


def random_effects(g: np.ndarray, meta: pd.DataFrame, col: str):
    effects = []
    for _, sub in meta.assign(_g=g).groupby(col, dropna=False):
        a, b, c, d = counts(sub["_g"].to_numpy(int), sub["y"].to_numpy(int))
        if (a+c) == 0 or (b+d) == 0 or (a+b) == 0 or (c+d) == 0:
            continue
        aa, bb, cc, dd = map(float, (a,b,c,d))
        if min(aa,bb,cc,dd) == 0:
            aa += .5; bb += .5; cc += .5; dd += .5
        lor = math.log((aa*dd)/(bb*cc))
        var = 1/aa + 1/bb + 1/cc + 1/dd
        effects.append((lor,var))
    if len(effects) < 2:
        return {"n":len(effects),"or":np.nan,"lo":np.nan,"hi":np.nan,"p":np.nan,"positive":0,"I2":np.nan}
    yy = np.array([e[0] for e in effects]); vv = np.array([e[1] for e in effects]); w = 1/vv
    mu = np.sum(w*yy)/np.sum(w); q = np.sum(w*(yy-mu)**2); df = len(yy)-1
    cval = np.sum(w) - np.sum(w*w)/np.sum(w); tau = max(0,(q-df)/cval) if df>0 and cval>0 else 0
    wr = 1/(vv+tau); mur = np.sum(wr*yy)/np.sum(wr); se = math.sqrt(1/np.sum(wr))
    return {"n":len(effects),"or":math.exp(mur),"lo":math.exp(mur-1.96*se),"hi":math.exp(mur+1.96*se),"p":float(2*norm.sf(abs(mur/se))),"positive":int((yy>0).sum()),"I2":max(0,(q-df)/q*100) if q>0 else 0}


def main() -> None:
    a = args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    sel = pd.read_csv(a.selection, dtype=str)
    if sel.empty:
        summary={"status":"NO_DISCOVERY_CANDIDATE","n_strict_statistically_replicated":0,"boundary":"No marker claim is supported."}
        (out/"UNITIG_REPLICATION_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n")
        return
    sel["sequence"] = sel["sequence"].map(normalize_seq)
    manifest = pd.read_csv(a.manifest, dtype=str).drop_duplicates("assembly_ID")
    manifest["assembly_ID"] = manifest["assembly_ID"].map(normalize_acc)
    manifest["phenotype"] = manifest["phenotype"].str.upper(); manifest = manifest[manifest["phenotype"].isin(["R","S"])].copy()
    manifest["y"] = manifest["phenotype"].eq("R").astype(int)
    for col in ["source_group","ISO_country_code"]:
        if col not in manifest: manifest[col]="UNKNOWN"
        manifest[col]=manifest[col].fillna("UNKNOWN").astype(str)

    ids = manifest["assembly_ID"].tolist(); sample_set=set(ids); selected_set=set(sel["sequence"])
    rtab = load_rtab(a.rtab, sample_set, selected_set)
    # Prefer exact sequence rows. If the query tool preserved FASTA IDs, restore by row order only when dimensions match exactly.
    if not selected_set.issubset(set(rtab.index)):
        if len(rtab) == len(sel):
            rtab.index = sel["sequence"].tolist()
        else:
            missing = sorted(selected_set-set(rtab.index))
            raise RuntimeError(f"Rtab cannot be mapped to {len(missing)} selected sequences")
    rtab = rtab.loc[sel["sequence"].tolist()]
    base={c.split('.')[0]:c for c in rtab.columns}; rename={}
    for sid in ids:
        if sid not in rtab.columns and sid.split('.')[0] in base: rename[base[sid.split('.')[0]]]=sid
    rtab=rtab.rename(columns=rename)
    missing=sorted(set(ids)-set(rtab.columns))
    if missing: raise RuntimeError(f"Rtab missing {len(missing)} samples")
    rtab=rtab.loc[:,ids]

    dist=load_distance(a.distance,ids); pcs=pcoa(dist.to_numpy(),30); pcdf=pd.DataFrame(pcs,index=ids)
    # Mash clusters for an additional within-lineage replication gate.
    condensed=squareform(dist.to_numpy(),checks=False); link=linkage(condensed,method='average'); manifest['mash_cluster']=fcluster(link,t=0.005,criterion='distance').astype(str)

    results={}; long=[]
    for scope,sub in [('discovery',manifest[manifest['split'].eq('discovery')]),('validation',manifest[manifest['split'].eq('validation')]),('all',manifest)]:
        sids=sub['assembly_ID'].tolist(); y=sub.set_index('assembly_ID').loc[sids,'y'].to_numpy(int); G=rtab.loc[:,sids].T.to_numpy(float); pc=pcdf.loc[sids].to_numpy()
        frames=[]
        for dim in [10,20,30]:
            z=score(y,pc[:,:dim],G); z['sequence']=sel['sequence'].tolist(); z['dim']=dim; z['q']=bh(z['p']); frames.append(z); long.append(z.assign(scope=scope))
        z=pd.concat(frames,ignore_index=True)
        wide=[]
        for seq,ss in z.groupby('sequence'):
            wide.append({'sequence':seq,**{f'{k}_{int(r.dim)}':getattr(r,k) for _,r in ss.iterrows() for k in ['beta','p','q']}})
        results[scope]=pd.DataFrame(wide)

    merged=sel.merge(results['discovery'],on='sequence').merge(results['validation'],on='sequence',suffixes=('_disc','_val'))
    # Rename ambiguous columns created by merge, then add whole cohort explicitly.
    whole=results['all'].copy(); whole=whole.rename(columns={c:f'{c}_all' for c in whole.columns if c!='sequence'})
    merged=merged.merge(whole,on='sequence')
    # After the first merge, discovery columns carry _disc and validation _val.
    for scope,suffix in [('disc','_disc'),('val','_val'),('all','_all')]:
        betas=[f'beta_{d}{suffix}' for d in [10,20,30]]; ps=[f'p_{d}{suffix}' for d in [10,20,30]]; qs=[f'q_{d}{suffix}' for d in [10,20,30]]
        merged[f'{scope}_positive_all']=np.logical_and.reduce([pd.to_numeric(merged[c],errors='coerce')>0 for c in betas])
        merged[f'{scope}_p_max']=merged[ps].apply(pd.to_numeric,errors='coerce').max(axis=1)
        merged[f'{scope}_q_max']=merged[qs].apply(pd.to_numeric,errors='coerce').max(axis=1)

    nsel=max(1,len(merged))
    merged['discovery_adjusted_stable']=merged['disc_positive_all'] & (merged['disc_q_max']<=a.alpha)
    merged['validation_bonferroni']=merged['val_positive_all'] & (merged['val_p_max']<=a.alpha/nsel)
    merged['whole_adjusted_stable']=merged['all_positive_all'] & (merged['all_q_max']<=a.alpha)

    detail=[]; Gfull=rtab.T.to_numpy(int); seqs=sel['sequence'].tolist(); seqidx={s:i for i,s in enumerate(seqs)}
    for seq in seqs:
        g=Gfull[:,seqidx[seq]]
        vals=manifest['split'].eq('validation').to_numpy(); va,vb,vc,vd=counts(g[vals],manifest.loc[vals,'y'].to_numpy(int)); vor,vlo,vhi=orci(va,vb,vc,vd)
        src=random_effects(g,manifest,'source_group'); country=random_effects(g,manifest,'ISO_country_code'); lineage=random_effects(g,manifest,'mash_cluster')
        detail.append({'sequence':seq,'validation_R_present':va,'validation_S_present':vb,'validation_R_absent':vc,'validation_S_absent':vd,'validation_or':vor,'validation_ci_low':vlo,'validation_ci_high':vhi,
                       'source_n':src['n'],'source_or':src['or'],'source_ci_low':src['lo'],'source_ci_high':src['hi'],'source_p':src['p'],'source_I2':src['I2'],
                       'country_n':country['n'],'country_or':country['or'],'country_ci_low':country['lo'],'country_ci_high':country['hi'],'country_p':country['p'],'country_I2':country['I2'],
                       'lineage_n':lineage['n'],'lineage_or':lineage['or'],'lineage_ci_low':lineage['lo'],'lineage_ci_high':lineage['hi'],'lineage_p':lineage['p'],'lineage_I2':lineage['I2']})
    merged=merged.merge(pd.DataFrame(detail),on='sequence')
    merged['source_replication']=(merged['source_n']>=3)&(merged['source_ci_low']>1)&(merged['source_p']<=a.alpha)
    merged['country_replication']=(merged['country_n']>=3)&(merged['country_ci_low']>1)&(merged['country_p']<=a.alpha)
    merged['lineage_replication']=(merged['lineage_n']>=3)&(merged['lineage_ci_low']>1)&(merged['lineage_p']<=a.alpha)
    merged['strict_statistical_replication']=merged['discovery_adjusted_stable']&merged['validation_bonferroni']&merged['whole_adjusted_stable']&merged['source_replication']&merged['country_replication']&merged['lineage_replication']
    merged['evidence_score']=3*merged['discovery_adjusted_stable'].astype(int)+5*merged['validation_bonferroni'].astype(int)+3*merged['whole_adjusted_stable'].astype(int)+3*merged['source_replication'].astype(int)+2*merged['country_replication'].astype(int)+3*merged['lineage_replication'].astype(int)
    merged=merged.sort_values(['strict_statistical_replication','evidence_score','val_p_max','disc_p_max'],ascending=[False,False,True,True])
    merged.to_csv(out/'ALL_UNITIG_REPLICATION_EVIDENCE.csv',index=False); strict=merged[merged['strict_statistical_replication']].copy(); strict.to_csv(out/'STRICT_STATISTICALLY_REPLICATED_UNITIGS.csv',index=False); pd.concat(long,ignore_index=True).to_csv(out/'ALL_SCORE_TESTS_LONG.csv',index=False)
    with open(out/'strict_unitigs.fasta','w') as fh:
        for _,r in strict.iterrows(): fh.write(f">{r['candidate_id']}\n{r['sequence']}\n")
    summary={'status':'STRICT_STATISTICAL_SURVIVORS_REQUIRE_MECHANISM_AND_NOVELTY_AUDIT' if len(strict) else 'NO_UNITIG_SURVIVED_COMPLETE_STATISTICAL_GATE','n_selected':len(sel),'n_evaluated':len(merged),'n_discovery_adjusted':int(merged['discovery_adjusted_stable'].sum()),'n_validation_bonferroni':int(merged['validation_bonferroni'].sum()),'n_whole_adjusted':int(merged['whole_adjusted_stable'].sum()),'n_source_replicated':int(merged['source_replication'].sum()),'n_country_replicated':int(merged['country_replication'].sum()),'n_lineage_replicated':int(merged['lineage_replication'].sum()),'n_strict_statistically_replicated':len(strict),'strict_candidate_ids':strict['candidate_id'].astype(str).tolist(),'boundary':'Statistical replication is not novelty or causality. Known-mechanism, sequence-context, database/literature and laboratory audits remain mandatory.'}
    (out/'UNITIG_REPLICATION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    cols=[c for c in ['candidate_id','sequence_length','one_sided_fisher_p','bonferroni_p','disc_p_max','val_p_max','all_p_max','validation_or','validation_ci_low','validation_ci_high','source_n','source_or','source_ci_low','country_n','country_or','country_ci_low','lineage_n','lineage_or','lineage_ci_low','strict_statistical_replication'] if c in merged]
    report=['# Source-held-out whole-genome unitig GWAS','',*[f'- {k}: **{v}**' for k,v in summary.items() if k not in {'boundary','strict_candidate_ids'}],'','## Claim boundary','',summary['boundary'],'','## Top evidence-ranked candidates','',merged.head(40)[cols].to_markdown(index=False)]
    (out/'UNITIG_REPLICATION_REPORT.md').write_text('\n'.join(report)+'\n')
    hashes=[]
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name!='SHA256SUMS.txt': hashes.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}')
    (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
