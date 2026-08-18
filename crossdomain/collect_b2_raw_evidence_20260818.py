#!/usr/bin/env python3
import base64, concurrent.futures, hashlib, json, os, pathlib, subprocess, time, urllib.request, urllib.parse, zipfile
R=os.environ['GITHUB_REPOSITORY']; T=os.environ['GH_TOKEN']; A='https://api.github.com'; O=pathlib.Path('B2_RAW_COLLECTED')
H={'Authorization':f'Bearer {T}','Accept':'application/vnd.github+json','X-GitHub-Api-Version':'2022-11-28','User-Agent':'b2-raw-collector'}
RUNS={'g1_v10':(32040699080,256),'g1_v8_head':(32038183046,9),'g1_v8_tail':(32038657803,35),'g23_v9':(32040180627,67)}
DEST={'g1_v10':'evidence/g1/v10_final_t2_t4','g1_v8_head':'evidence/g1/v8_head_t2_t10','g1_v8_tail':'evidence/g1/v8_tail_t11_t45','g23_v9':'evidence/g2_g3/v9'}
DIRECT={'summaries/g1_v10':9292041311,'summaries/g23_v9':9291987535,'evidence/g4_g7_drat/g4':9279813096,'evidence/g4_g7_drat/g5':9291076270,'evidence/g4_g7_drat/g6':9291072827,'evidence/g4_g7_drat/g7':9291073263}

def get(url):
  last=None
  for i in range(15):
    try:
      with urllib.request.urlopen(urllib.request.Request(url,headers=H),timeout=180) as r:return json.load(r)
    except Exception as e:last=e;time.sleep(min(90,5*(i+1)))
  raise RuntimeError((url,last))
def pages(url,key):
  out=[];p=1
  while 1:
    d=get(f'{url}{"&" if "?" in url else "?"}per_page=100&page={p}');x=d.get(key,[]);out+=x
    if len(x)<100:return out
    p+=1
def sha(p):
  h=hashlib.sha256()
  with open(p,'rb') as f:
    for c in iter(lambda:f.read(1048576),b''):h.update(c)
  return h.hexdigest()
def dl(art,dest):
  dest.mkdir(parents=True,exist_ok=True);z=dest.parent/f'.{dest.name}.zip';tmp=pathlib.Path(str(z)+'.part')
  cmd=['curl','-fL','--retry','30','--retry-all-errors','--retry-delay','5','--connect-timeout','60','--max-time','7200','-H',f'Authorization: Bearer {T}','-H','Accept: application/vnd.github+json',art['archive_download_url'],'-o',str(tmp)]
  subprocess.run(cmd,check=True);tmp.replace(z)
  with zipfile.ZipFile(z) as q:
    bad=q.testzip();assert bad is None,(art['name'],bad);q.extractall(dest)
  m={'id':art['id'],'name':art['name'],'api_size':art.get('size_in_bytes'),'api_digest':art.get('digest'),'archive_bytes':z.stat().st_size,'archive_sha256':sha(z),'workflow_run':art.get('workflow_run')}
  (dest/'ARTIFACT_METADATA.json').write_text(json.dumps(m,indent=2));z.unlink();return m

def provenance(label,run,jobs):
  d=O/'environment'/f'{label}_{run["id"]}';d.mkdir(parents=True,exist_ok=True)
  (d/'RUN.json').write_text(json.dumps(run,indent=2));(d/'JOBS.json').write_text(json.dumps(jobs,indent=2))
  if run.get('path') and run.get('head_sha'):
    c=get(f'{A}/repos/{R}/contents/{urllib.parse.quote(run["path"],safe="/")}?ref={run["head_sha"]}')
    (d/'WORKFLOW_AT_RUN.yml').write_bytes(base64.b64decode(c['content']))

idx=[]
for label,(rid,n) in RUNS.items():
  run=get(f'{A}/repos/{R}/actions/runs/{rid}');jobs=pages(f'{A}/repos/{R}/actions/runs/{rid}/jobs','jobs');arts=pages(f'{A}/repos/{R}/actions/runs/{rid}/artifacts','artifacts');assert len(arts)==n,(label,len(arts),n);provenance(label,run,jobs)
  root=O/DEST[label]
  with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
    fs=[ex.submit(dl,a,root/a['name']) for a in arts]
    for f in concurrent.futures.as_completed(fs):idx.append(f.result())
for rel,aid in DIRECT.items():
  a=get(f'{A}/repos/{R}/actions/artifacts/{aid}');idx.append(dl(a,O/rel/a['name']))
(O/'manifest').mkdir(parents=True,exist_ok=True);(O/'manifest'/'ARTIFACT_INDEX.json').write_text(json.dumps(idx,indent=2))
for p in ['crossdomain/b2_gap1_u_split_v8.cpp','crossdomain/b2_gap1_v_split_v10.cpp','crossdomain/b2_gap_g_u_split_v9.cpp','crossdomain/b2_gap_sat_v4.py','crossdomain/B2_Z100_FINAL_RESOLUTION_20260817.md','crossdomain/B2_Z100_FINAL_RESULTS_20260817.json','crossdomain/B2_Z100_VERIFY_FINAL_20260817.py']:
  q=pathlib.Path(p)
  if q.exists():
    d=O/'repository_snapshot'/q.parent;d.mkdir(parents=True,exist_ok=True);(d/q.name).write_bytes(q.read_bytes())
for q in pathlib.Path('.github/workflows').glob('b2*.yml'):
  d=O/'repository_snapshot/.github/workflows';d.mkdir(parents=True,exist_ok=True);(d/q.name).write_bytes(q.read_bytes())
subprocess.run('uname -a; cat /etc/os-release; g++ --version; python --version; lscpu',shell=True,text=True,stdout=open(O/'environment'/'PACKAGING_ENV.txt','w'),stderr=subprocess.STDOUT)
print(json.dumps({'artifacts':len(idx),'files':sum(x.is_file() for x in O.rglob('*')),'bytes':sum(x.stat().st_size for x in O.rglob('*') if x.is_file())},indent=2))
