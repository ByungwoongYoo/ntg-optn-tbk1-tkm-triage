#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re
from collections import Counter,defaultdict
from pathlib import Path

FORMS='丸|圓|散|湯|丹|膏|飲|飮|煎|酒|餅|粉|劑|方'
DETAIL_RE=re.compile(r'^\s*([\u3400-\u9fff𠀀-𫠝]{1,12}(?:'+FORMS+r'))[，,。:：\s]*(?=(?:治|療|主|$))')
NAME_TOKEN_RE=re.compile(r'([\u3400-\u9fff𠀀-𫠝]{1,12}(?:'+FORMS+r'))')
OTHER_SOURCES=['壽域神方','衛生易簡方','施圓端效方','吳氏集驗方','神效名方','簡方','金翼方','玉機微義','急救仙方','壽親養老書','醫林方','新效方']

def norm(s):
    return re.sub(r'[\s，,。:：；;、「」『』（）()]+','',s)

def clean_name(s):
    s=norm(s)
    # conservative: no character correction; only trim obvious leading list punctuation/noise digits
    s=re.sub(r'^[0-9一二三四五六七八九十卜山木金火水土]+(?=.{2,})','',s)
    return s

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--source-blocks',required=True); ap.add_argument('--out-dir',required=True); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    blocks=json.loads(Path(a.source_blocks).read_text(encoding='utf-8'))
    detail=[]; summary=[]; source_flags=[]
    for bi,b in enumerate(blocks):
        lines=[x.strip() for x in b['text'].splitlines() if x.strip()]
        for li,line in enumerate(lines):
            m=DETAIL_RE.match(line)
            if m:
                raw=m.group(1); name=clean_name(raw)
                flags=[s for s in OTHER_SOURCES if s in line or (li>0 and s in lines[li-1])]
                detail.append({'name':name,'raw_name':raw,'block':bi,'line_index':li,'line':line,'context':' | '.join(lines[max(0,li-1):min(len(lines),li+3)]),'url':b['url'],'other_source_flags':flags})
            # summary-line heuristic: two or more dosage-form names in a short list-like paragraph
            toks=[clean_name(x) for x in NAME_TOKEN_RE.findall(line)]
            toks=[x for x in toks if 2<=len(x)<=14]
            if len(set(toks))>=2 and len(line)<=120 and not line.startswith('右'):
                for x in dict.fromkeys(toks):
                    summary.append({'name':x,'block':bi,'line_index':li,'line':line,'url':b['url']})
            for s in OTHER_SOURCES:
                if s in line:
                    source_flags.append({'source':s,'block':bi,'line_index':li,'line':line,'url':b['url']})
    dcount=Counter(x['name'] for x in detail); scount=Counter(x['name'] for x in summary)
    names=sorted(set(dcount)|set(scount))
    inv=[]
    for n in names:
        examples=[x for x in detail if x['name']==n][:2]
        inv.append({'name':n,'detail_heading_count':dcount[n],'summary_list_count':scount[n],'strong_detail':dcount[n]>0,'cross_source_flagged_occurrences':sum(bool(x['other_source_flags']) for x in detail if x['name']==n),'examples':examples})
    inv.sort(key=lambda x:(not x['strong_detail'],-x['detail_heading_count'],-x['summary_list_count'],x['name']))
    strong=[x for x in inv if x['strong_detail']]
    result={
      'target':'煙霞聖效方',
      'input_blocks':len(blocks),
      'strong_unique_formula_headings':len(strong),
      'detail_heading_occurrences':len(detail),
      'summary_name_occurrences':len(summary),
      'summary_unique_names':len(set(x['name'] for x in summary)),
      'strong_names_also_seen_in_summary':sum(1 for x in strong if x['summary_list_count']>0),
      'other_source_flag_occurrences':len(source_flags),
      'top_strong_names':[{'name':x['name'],'detail_heading_count':x['detail_heading_count'],'summary_list_count':x['summary_list_count'],'cross_source_flagged_occurrences':x['cross_source_flagged_occurrences']} for x in strong[:100]],
      'claim_boundary':'This is a conservative machine inventory of candidate formula headings inside source-attributed Yanxia Shengxiaofang blocks. It does not correct OCR, establish exact historical boundaries, or prove that every listed formula originated uniquely in the lost book. Other-source flags identify places needing philological adjudication.'
    }
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'formula_inventory.json').write_text(json.dumps(inv,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'detail_headings.json').write_text(json.dumps(detail,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'summary_names.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'other_source_flags.json').write_text(json.dumps(source_flags,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
