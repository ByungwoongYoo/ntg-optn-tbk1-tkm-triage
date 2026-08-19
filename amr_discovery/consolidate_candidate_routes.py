#!/usr/bin/env python3
"""Consolidate strict candidate files across independently executed AMR routes.

This script never creates a candidate from a non-strict table. It records which
prespecified route(s) independently promoted the same marker and applies a
cross-method gate before any novelty audit. A statistical promotion remains only
an association candidate, never a causal resistance mechanism.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
import pandas as pd


def args():
    p=argparse.ArgumentParser(); p.add_argument('--artifacts-root',required=True); p.add_argument('--out',required=True); return p.parse_args()

def route(path: Path) -> str:
    s=str(path).lower()
    if 'score_gwas' in s or 'structure_adjusted_score' in s: return 'targeted_score'
    if 'targeted_portal_gwas' in s or 'final_evaluation_v3' in s: return 'targeted_pyseer'
    if 'lineage' in s or 'cmh' in s: return 'lineage_cmh'
    if 'rare' in s or 'burden' in s: return 'rare_burden'
    if 'unitig' in s: return 'unitig_gwas'
    return 'unknown'

def candidate_column(df: pd.DataFrame):
    for c in ['feature','variant','candidate_id','burden_id','canonical_sequence','gene']:
        if c in df.columns:return c
    return None

def strict_files(root: Path):
    for p in root.rglob('*.csv'):
        u=p.name.upper()
        if 'STRICT' in u and any(k in u for k in ['FEATURE','UNITIG','BURDEN','REPLICATE','CANDIDATE']): yield p

def main():
    a=args(); root=Path(a.artifacts_root); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    rows=[]; audit=[]
    for p in sorted(strict_files(root)):
        try: df=pd.read_csv(p)
        except Exception as e:
            audit.append({'path':str(p),'route':route(p),'read_error':repr(e),'n_rows':None}); continue
        c=candidate_column(df); sha=hashlib.sha256(p.read_bytes()).hexdigest(); audit.append({'path':str(p),'route':route(p),'sha256':sha,'n_rows':len(df),'candidate_column':c})
        if c is None: continue
        for _,r in df.iterrows():
            if pd.isna(r[c]): continue
            cand=str(r[c]).strip()
            if not cand: continue
            row={'candidate':cand,'route':route(p),'source_file':str(p),'source_sha256':sha}
            for key in ['gene','burden_id','candidate_id','canonical_sequence','component_features','validation_or','validation_ci_low','validation_ci_high','validation_q','source_ci_low','country_ci_low','q_max','p_max']:
                if key in r.index: row[key]=r[key]
            rows.append(row)
    pd.DataFrame(audit).to_csv(out/'STRICT_FILE_AUDIT.csv',index=False)
    raw=pd.DataFrame(rows)
    raw.to_csv(out/'STRICT_ROUTE_ROWS.csv',index=False)
    if raw.empty:
        summary={'status':'NO_STRICT_CANDIDATE_IN_ANY_COMPLETED_ROUTE','n_unique_candidates':0,'n_cross_method_promoted':0,'claim_boundary':'No genomic marker may be announced.'}
        pd.DataFrame(columns=['candidate','routes','n_routes','candidate_class','cross_method_promoted']).to_csv(out/'CONSOLIDATED_CANDIDATES.csv',index=False)
        (out/'CONSOLIDATION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2)); return
    grouped=[]
    for cand,g in raw.groupby('candidate'):
        routes=sorted(set(g.route)); sequence_like=bool(re.fullmatch('[ACGTNacgtn]+',cand))
        is_burden='rare_burden' in routes or '|' in cand and 'rare_le_' in cand
        if sequence_like or 'unitig_gwas' in routes: cls='unitig'
        elif is_burden: cls='gene_burden'
        else: cls='targeted_variant'
        if cls=='targeted_variant': promoted=('targeted_score' in routes and 'targeted_pyseer' in routes)
        elif cls=='unitig': promoted=('unitig_gwas' in routes)
        else: promoted=('rare_burden' in routes)
        # Lineage CMH is a robustness route, not a substitute for the independent targeted implementation.
        grouped.append({'candidate':cand,'routes':';'.join(routes),'n_routes':len(routes),'candidate_class':cls,'cross_method_promoted':promoted,'lineage_supported':'lineage_cmh' in routes,'source_files':';'.join(sorted(set(g.source_file.astype(str))))})
    consolidated=pd.DataFrame(grouped).sort_values(['cross_method_promoted','n_routes','candidate'],ascending=[False,False,True])
    consolidated.to_csv(out/'CONSOLIDATED_CANDIDATES.csv',index=False)
    promoted=consolidated[consolidated.cross_method_promoted].copy(); promoted.to_csv(out/'CROSS_METHOD_PROMOTED_CANDIDATES.csv',index=False)
    summary={'status':'CROSS_METHOD_CANDIDATES_REQUIRE_MECHANISM_NOVELTY_AND_CONTEXT_AUDIT' if len(promoted) else 'NO_CANDIDATE_SURVIVED_CROSS_METHOD_GATE','n_unique_candidates':len(consolidated),'n_cross_method_promoted':len(promoted),'promoted_candidates':promoted.candidate.astype(str).tolist(),'claim_boundary':'Cross-method promotion establishes neither novelty nor causality. Known-mechanism, lineage, genomic-context, independent-database and literature audits remain mandatory.'}
    (out/'CONSOLIDATION_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    report=['# Cross-method AMR candidate consolidation','',f"- Unique strict-route candidates: **{len(consolidated):,}**",f"- Cross-method promoted candidates: **{len(promoted):,}**",'',summary['claim_boundary'],'','## Consolidated route evidence','',consolidated.to_markdown(index=False)]
    (out/'CONSOLIDATION_REPORT.md').write_text('\n'.join(report)+'\n')
    hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']; (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n'); print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__': main()
