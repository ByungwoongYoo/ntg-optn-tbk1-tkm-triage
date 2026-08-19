#!/usr/bin/env python3
"""Search primary bibliographic indexes for exact candidate/gene/mutation terms.

A zero-hit automated search is not proof of novelty. Results are saved for manual
review and timestamped so the search can be repeated.
"""
from __future__ import annotations
import argparse,json,re,time,urllib.parse,urllib.request
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def args():
    p=argparse.ArgumentParser();p.add_argument('--classification',required=True);p.add_argument('--out',required=True);p.add_argument('--email',default='');return p.parse_args()

def get_json(url,attempts=4):
    headers={'User-Agent':'Kp-colistin-candidate-audit/1.0'}
    for i in range(attempts):
        try:
            req=urllib.request.Request(url,headers=headers)
            with urllib.request.urlopen(req,timeout=60) as r:return json.load(r)
        except Exception:
            if i==attempts-1:raise
            time.sleep(2**i)

def parse_mutation(candidate):
    m=re.search(r'([A-Z*])([0-9]+)([A-Z*])',str(candidate).upper())
    return m.group(0) if m else ''

def epmc(query):
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':query,'format':'json','pageSize':100,'resultType':'core'})
    d=get_json(url);result=d.get('resultList',{}).get('result',[])
    return {'url':url,'hit_count':int(d.get('hitCount',0)),'results':[{'id':r.get('id'),'source':r.get('source'),'title':r.get('title'),'authorString':r.get('authorString'),'journalTitle':r.get('journalTitle'),'pubYear':r.get('pubYear'),'doi':r.get('doi'),'pmid':r.get('pmid'),'pmcid':r.get('pmcid')} for r in result[:100]]}

def pubmed(query,email=''):
    params={'db':'pubmed','term':query,'retmode':'json','retmax':100,'tool':'kp_colistin_candidate_audit'}
    if email:params['email']=email
    search_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urllib.parse.urlencode(params)
    d=get_json(search_url);ids=d.get('esearchresult',{}).get('idlist',[]);count=int(d.get('esearchresult',{}).get('count',0))
    records=[]
    if ids:
        summary_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),'retmode':'json','tool':'kp_colistin_candidate_audit',**({'email':email} if email else {})})
        s=get_json(summary_url)
        for pmid in ids:
            r=s.get('result',{}).get(pmid,{})
            records.append({'pmid':pmid,'title':r.get('title'),'pubdate':r.get('pubdate'),'fulljournalname':r.get('fulljournalname'),'authors':r.get('authors')})
    return {'url':search_url,'hit_count':count,'results':records}

def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);df=pd.read_csv(a.classification,dtype=str) if Path(a.classification).exists() else pd.DataFrame();all_results=[];summary=[]
    for _,r in df.iterrows():
        candidate=str(r.get('candidate',''));gene=str(r.get('gene',''));mutation=parse_mutation(candidate);queries=[]
        if mutation:queries.append(f'"Klebsiella pneumoniae" AND colistin AND {gene} AND {mutation}')
        queries.append(f'"Klebsiella pneumoniae" AND colistin AND {gene}')
        if str(r.get('candidate_class'))=='unitig':queries.append('"Klebsiella pneumoniae" AND colistin AND genomic marker')
        candidate_records=[]
        for q in list(dict.fromkeys(queries)):
            er=epmc(q);time.sleep(.4);pr=pubmed(q,a.email);time.sleep(.4);candidate_records.append({'query':q,'europe_pmc':er,'pubmed':pr})
        exact_hits=0
        if mutation:
            for record in candidate_records[:1]:exact_hits=max(record['europe_pmc']['hit_count'],record['pubmed']['hit_count'])
        all_results.append({'candidate':candidate,'gene':gene,'mutation':mutation,'searches':candidate_records})
        summary.append({'candidate':candidate,'gene':gene,'mutation':mutation,'exact_query_hit_count':exact_hits,'gene_query_epmc_hits':candidate_records[-1]['europe_pmc']['hit_count'],'gene_query_pubmed_hits':candidate_records[-1]['pubmed']['hit_count'],'automated_exact_hit_found':exact_hits>0})
    payload={'searched_at_utc':datetime.now(timezone.utc).isoformat(),'boundary':'Zero hits do not establish novelty; spelling, numbering, reference-strain and indexing differences require manual review.','candidates':all_results}
    (out/'CANDIDATE_LITERATURE_SEARCH_FULL.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False,default=str)+'\n')
    pd.DataFrame(summary).to_csv(out/'CANDIDATE_LITERATURE_SEARCH_SUMMARY.csv',index=False)
    (out/'SEARCH_BOUNDARY.txt').write_text(payload['boundary']+'\n');print(json.dumps({'n_candidates':len(summary),'n_exact_query_hits':sum(x['automated_exact_hit_found'] for x in summary),'boundary':payload['boundary']},indent=2))
if __name__=='__main__':main()
