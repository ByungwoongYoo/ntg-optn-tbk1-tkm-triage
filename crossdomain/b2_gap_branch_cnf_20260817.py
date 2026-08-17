#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from crossdomain.b2_gap_sat_v4 import build

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--gap',type=int,required=True)
    ap.add_argument('--third',type=int,required=True)
    ap.add_argument('--out-dir',required=True)
    a=ap.parse_args(); g=a.gap; t=a.third
    lo=2*g; hi=50-5*g
    if not (1<=g<=7 and lo<=t<=hi):
        raise SystemExit(f'invalid normalized branch: g={g}, third={t}, expected {lo}..{hi}')
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    cnf,x,_=build(g)
    # t is the first selected residue strictly after g.
    for i in range(g+1,t): cnf.append([-x[i]])
    cnf.append([x[t]])
    # Reflection x -> g-x swaps the next gap h=t-g and the wrap gap w=100-last.
    # WLOG h<=w, equivalently last<=100-t+g.
    upper=100-t+g
    for i in range(upper+1,100): cnf.append([-x[i]])
    cnf_path=out/f'g{g}_t{t}.cnf'; cnf.to_file(cnf_path)
    meta={'problem':'14-element B_2[2] subset of Z_100','gap':g,'third_selected':t,
          'third_range':[lo,hi],'reflection_upper_selected_index':upper,
          'coverage':'For any set in normalized gap g, let h be the gap after (0,g) and w the wrap gap before 0. Reflection x->g-x swaps h and w, so choose orientation h<=w. With the other eleven non-wrap gaps at least g, 100>=12g+2h, hence t=g+h lies in [2g,50-5g]. The listed branches are exhaustive up to reflection.',
          'cnf_vars':cnf.nv,'cnf_clauses':len(cnf.clauses),'cnf_file':cnf_path.name}
    (out/'BRANCH_METADATA.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps(meta,indent=2))
if __name__=='__main__': main()
