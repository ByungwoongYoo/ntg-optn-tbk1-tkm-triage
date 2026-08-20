#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--report',required=True);p.add_argument('--method',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 rows=list(csv.reader(open(a.report),delimiter='\t'));header=rows[0];method_col=1 if len(header)>1 else None;data={}
 for row in rows[1:]:
  if not row:continue
  key=row[0].strip();val=row[method_col].strip() if method_col is not None and len(row)>method_col else '';data[key]=val
 aliases={'genome_fraction':'Genome fraction (%)','duplication_ratio':'Duplication ratio','#_misassemblies':'# misassemblies','misassembled_contigs_length':'Misassembled contigs length','n50':'N50','nga50':'NGA50','largest_contig':'Largest contig','total_length':'Total length (>= 0 bp)','#_contigs':'# contigs (>= 0 bp)','mismatches_per_100kb':'# mismatches per 100 kbp','indels_per_100kb':'# indels per 100 kbp'}
 out={'method':a.method}
 for k,label in aliases.items():
  v=data.get(label,'')
  try:out[k]=float(v.replace(',',''))
  except:out[k]=None
 Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
