#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,statistics
from collections import defaultdict
from pathlib import Path

def args():
    p=argparse.ArgumentParser();p.add_argument('--abundance',action='append',required=True);p.add_argument('--recovery',action='append',required=True,help='METHOD=per_genome_recovery.tsv');p.add_argument('--low-percent',type=float,default=.01);p.add_argument('--high-percent',type=float,default=.1);p.add_argument('--out',required=True);return p.parse_args()

def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);by=defaultdict(dict)
    for path in a.abundance:
        for r in csv.DictReader(open(path),delimiter='\t'):by[r['genome_id']][r['sample_id']]=float(r['relative_read_percent'])
    profile={}
    for gid,x in by.items():
        vals=list(x.values());mean=sum(vals)/len(vals);mn=min(vals);mx=max(vals)
        if mean<a.low_percent:stratum='low'
        elif mean<a.high_percent:stratum='medium'
        else:stratum='high'
        if mx<a.low_percent:trajectory='persistent_low'
        elif mn<a.low_percent<=mx:trajectory='longitudinal_rescue_opportunity'
        else:trajectory='consistently_detectable'
        profile[gid]={'genome_id':gid,'n_samples':len(vals),'mean_read_percent':mean,'min_read_percent':mn,'max_read_percent':mx,'abundance_stratum':stratum,'trajectory_stratum':trajectory,**{f'{s}_read_percent':v for s,v in sorted(x.items())}}
    details=[];summaries=[]
    for spec in a.recovery:
        method,path=spec.split('=',1)
        for r in csv.DictReader(open(path),delimiter='\t'):
            gid=r['genome_id']
            if gid not in profile:continue
            details.append({'method':method,**profile[gid],'recovery_fraction':float(r['recovery_fraction']),'recovered_50':r['recovered_50'],'recovered_90':r['recovered_90']})
    groups=['low','medium','high','persistent_low','longitudinal_rescue_opportunity','consistently_detectable']
    methods=sorted({r['method'] for r in details})
    for method in methods:
        for group in groups:
            if group in {'low','medium','high'}:rows=[r for r in details if r['method']==method and r['abundance_stratum']==group]
            else:rows=[r for r in details if r['method']==method and r['trajectory_stratum']==group]
            vals=[r['recovery_fraction'] for r in rows]
            summaries.append({'method':method,'stratum':group,'n_genomes':len(rows),'mean_recovery':sum(vals)/len(vals) if vals else None,'median_recovery':statistics.median(vals) if vals else None,'minimum_recovery':min(vals) if vals else None,'recovered_50':sum(r['recovered_50']=='true' for r in rows),'recovered_90':sum(r['recovered_90']=='true' for r in rows)})
    for name,rows in [('ABUNDANCE_RECOVERY_DETAIL.tsv',details),('ABUNDANCE_STRATIFIED_SUMMARY.tsv',summaries)]:
        fields=list(rows[0]) if rows else ['method','stratum','n_genomes','mean_recovery'];
        with open(out/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    report={'n_abundance_genomes':len(profile),'n_joined_rows':len(details),'low_threshold_percent':a.low_percent,'high_threshold_percent':a.high_percent,'boundary':'Abundance and genome labels are gold-standard evaluation metadata accessed after truth-blind output freeze. Low-abundance strata do not influence selection.'}
    (out/'ABUNDANCE_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
