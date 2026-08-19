#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,statistics
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser();p.add_argument('--recovery',action='append',required=True,help='METHOD=per_genome_recovery.tsv');p.add_argument('--out',required=True);return p.parse_args()

def domain(gid):
    x=gid.lower()
    if x.startswith('vasv') or 'virus' in x or x.startswith('phage'):return 'virus'
    if x.startswith('pasv') or 'plasmid' in x:return 'plasmid'
    if x.startswith('fasv') or 'fung' in x or 'candida' in x:return 'fungus'
    if x.startswith('hasv') or 'host' in x or 'homo' in x:return 'host'
    if 'archae' in x or x.startswith('aasv'):return 'archaea_candidate'
    return 'bacteria_or_archaea_unresolved'

def main():
    a=parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);detail=[]
    for spec in a.recovery:
        method,path=spec.split('=',1)
        for r in csv.DictReader(open(path),delimiter='\t'):detail.append({'method':method,'genome_id':r['genome_id'],'domain':domain(r['genome_id']),'recovery_fraction':float(r['recovery_fraction']),'recovered_50':r['recovered_50'],'recovered_90':r['recovered_90']})
    summary=[]
    for method in sorted({r['method'] for r in detail}):
        for dom in sorted({r['domain'] for r in detail}):
            rows=[r for r in detail if r['method']==method and r['domain']==dom];vals=[r['recovery_fraction'] for r in rows];summary.append({'method':method,'domain':dom,'n_genomes':len(rows),'mean_recovery':sum(vals)/len(vals) if vals else None,'median_recovery':statistics.median(vals) if vals else None,'recovered_50':sum(r['recovered_50']=='true' for r in rows),'recovered_90':sum(r['recovered_90']=='true' for r in rows)})
    for name,rows in [('DOMAIN_RECOVERY_DETAIL.tsv',detail),('DOMAIN_RECOVERY_SUMMARY.tsv',summary)]:
        fields=list(rows[0]) if rows else ['method','domain','n_genomes','mean_recovery'];
        with open(out/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    report={'n_detail_rows':len(detail),'domains':sorted({r['domain'] for r in detail}),'boundary':'Domain labels are conservative BINID-prefix heuristics. ASV records without explicit prefixes cannot separate bacteria from archaea and are reported together.'};(out/'DOMAIN_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
