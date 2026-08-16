#!/usr/bin/env python3
from __future__ import annotations

import json, math, os, time
from pathlib import Path
from urllib.parse import urlencode
import pandas as pd
import requests

OUT=Path(os.environ.get('OUT_DIR','safety_v2_output')); OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({'User-Agent':'OpenProblemTournamentSafetyV2/0.1 (yoonge3@gmail.com)'})
END='https://api.fda.gov/drug/event.json'

def phrase(field, value):
    value=value.replace('"','')
    return f'{field}:"{value}"'

def any_of(field, aliases):
    return '('+' OR '.join(phrase(field,x) for x in aliases)+')'

def count(search):
    r=S.get(END,params={'search':search,'limit':1},timeout=120)
    if r.status_code==404:return 0, r.url
    r.raise_for_status()
    return int(r.json()['meta']['results']['total']),r.url

def ror_table(pair,period,drange):
    drug=any_of('patient.drug.medicinalproduct',pair['drug'])
    herb=any_of('patient.drug.medicinalproduct',pair['herb'])
    event=any_of('patient.reaction.reactionmeddrapt',pair['event'])
    date=f'receivedate:[{drange[0]} TO {drange[1]}]'
    qs={
      'a':f'{drug} AND {herb} AND {event} AND {date}',
      'exposed':f'{drug} AND {herb} AND {date}',
      'drug_event':f'{drug} AND {event} AND {date}',
      'drug_total':f'{drug} AND {date}',
      'herb_total':f'{herb} AND {date}',
    }
    vals={}; urls={}
    for k,q in qs.items():
        vals[k],urls[k]=count(q); time.sleep(.15)
    a=vals['a']; b=max(vals['exposed']-a,0); c=max(vals['drug_event']-a,0); d=max(vals['drug_total']-vals['exposed']-c,0)
    aa,bb,cc,dd=a+.5,b+.5,c+.5,d+.5
    R=(aa/bb)/(cc/dd) if bb and cc and dd else None
    if R and R>0:
        se=math.sqrt(1/aa+1/bb+1/cc+1/dd); ci=(math.exp(math.log(R)-1.96*se),math.exp(math.log(R)+1.96*se))
    else:ci=(None,None)
    return {'pair':pair['name'],'period':period,**vals,'b':b,'c':c,'d':d,'ror':R,'ror_lo':ci[0],'ror_hi':ci[1],'queries':urls}

pairs=[
 {'name':'ginkgo_warfarin_haemorrhage','herb':['GINKGO','GINKGO BILOBA'],'drug':['WARFARIN','COUMADIN'],'event':['HAEMORRHAGE','GASTROINTESTINAL HAEMORRHAGE','INTERNATIONAL NORMALISED RATIO INCREASED']},
 {'name':'st_johns_wort_ssri_serotonin','herb':['ST JOHNS WORT','ST. JOHNS WORT','HYPERICUM'],'drug':['SERTRALINE','FLUOXETINE','PAROXETINE','CITALOPRAM','ESCITALOPRAM'],'event':['SEROTONIN SYNDROME']},
 {'name':'licorice_diuretic_hypokalaemia','herb':['LICORICE','LIQUORICE','GLYCYRRHIZA'],'drug':['FUROSEMIDE','HYDROCHLOROTHIAZIDE'],'event':['HYPOKALAEMIA','BLOOD POTASSIUM DECREASED']},
]
periods={'early':('20040101','20151231'),'late':('20160101','20260428')}
rows=[]; errors=[]
for p in pairs:
  for period,dr in periods.items():
    try: rows.append(ror_table(p,period,dr))
    except Exception as e: errors.append({'pair':p['name'],'period':period,'error':repr(e)})
flat=[{k:v for k,v in r.items() if k!='queries'} for r in rows]
pd.DataFrame(flat).to_csv(OUT/'safety_v2.csv',index=False)
rep=[]
for p in pairs:
  rr={x['period']:x for x in rows if x['pair']==p['name']}
  if set(rr)==set(periods):
    e,l=rr['early'],rr['late']
    if e['a']>=5 and l['a']>=5 and e['ror'] and l['ror'] and e['ror']>1 and l['ror']>1: rep.append(p['name'])
res={'rows':rows,'errors':errors,'replicated_directional_signals':rep,'interpretation':'Technical pharmacovigilance smoke only. FAERS reports do not establish causality or incidence; listed products and reactions are not linked one-to-one.'}
(OUT/'safety_v2.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
(OUT/'SUMMARY.md').write_text('# Herb–drug safety smoke v2\n\n'+f'- Successful queries: {len(rows)}/6\n- Errors: {len(errors)}\n- Directionally replicated prespecified signals: {rep}\n\nFAERS co-reporting is a signal-detection source, not causal or incidence evidence.\n',encoding='utf-8')
print(json.dumps({'n_rows':len(rows),'errors':errors,'replicated':rep},indent=2))
