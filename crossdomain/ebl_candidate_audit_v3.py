#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from collections import Counter
from pathlib import Path

ID_PATTERNS=[
    re.compile(r'\bK[.\s-]?(\d+)(?:[.\s-]?([A-Z]))?\b',re.I),
    re.compile(r'\bSm[.\s-]?(\d+)\b',re.I),
    re.compile(r'\bBM[.\s-]?(\d+)\b',re.I),
]


def canon(x):
    return re.sub(r'[^A-Z0-9]','',str(x or '').upper())


def pubnorm(x):
    return re.sub(r'\s+',' ',str(x or '')).strip()


def note_text(x):
    if isinstance(x,dict):
        return str(x.get('text') or '')
    return str(x or '')


def refs_in_notes(text):
    refs=set()
    for m in ID_PATTERNS[0].finditer(text):
        refs.add('K'+m.group(1)+(m.group(2) or '').upper())
    for m in ID_PATTERNS[1].finditer(text): refs.add('SM'+m.group(1))
    for m in ID_PATTERNS[2].finditer(text): refs.add('BM'+m.group(1))
    return refs


def dimval(d):
    try: return float((d or {}).get('value'))
    except Exception: return None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--fragments-json',required=True)
    ap.add_argument('--ranked-v2',required=True)
    ap.add_argument('--out-dir',required=True)
    a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)

    frags=json.loads(Path(a.fragments_json).read_text(encoding='utf-8'))
    if not isinstance(frags,list): frags=frags.get('fragments',[])
    by={f['_id']:f for f in frags}
    ranked=json.loads(Path(a.ranked_v2).read_text(encoding='utf-8'))

    audited=[]
    for c in ranked:
        fa,fb=by[c['a']],by[c['b']]
        pa,pb=pubnorm(fa.get('publication')),pubnorm(fb.get('publication'))
        na,nb=note_text(fa.get('notes')),note_text(fb.get('notes'))
        ca,cb=canon(c['a']),canon(c['b'])
        ra,rb=refs_in_notes(na),refs_in_notes(nb)
        counterpart_noted=(cb in ra or ca in rb)
        join_language=bool(re.search(r'\bjoin(?:s|ed)?\b|indirect join|direct join|\(\+\)',na+' '+nb,re.I))
        transfer_language=bool(re.search(r'edition transferred from|transferred from',na+' '+nb,re.I))
        ta,tb=dimval(fa.get('thickness')),dimval(fb.get('thickness'))
        td=abs(ta-tb) if ta is not None and tb is not None else None

        if pa and pb: publication_state='BOTH_HAVE_EDITION'
        elif pa or pb: publication_state='ONE_HAS_EDITION'
        else: publication_state='NO_EDITION_FIELD'

        if counterpart_noted:
            metadata_status='EXCLUDE_COUNTERPART_ALREADY_REFERENCED_IN_NOTES'
        else:
            metadata_status='SURVIVES_METADATA_SCREEN'

        pub_factor={'NO_EDITION_FIELD':1.0,'ONE_HAS_EDITION':0.9,'BOTH_HAVE_EDITION':0.7}[publication_state]
        if td is None: geom_factor=0.85
        elif td<=0.30: geom_factor=1.10
        elif td<=1.00: geom_factor=1.0
        else: geom_factor=0.25
        note_factor=0.85 if join_language else 1.0
        if transfer_language: note_factor*=0.75
        if counterpart_noted: note_factor=0.0
        score=float(c['sign_tfidf_similarity'])*pub_factor*geom_factor*note_factor

        audited.append({
            **c,
            'a_publication_v3':pa,'b_publication_v3':pb,
            'publication_state_v3':publication_state,
            'a_note_text':na[:2500],'b_note_text':nb[:2500],
            'a_note_refs':sorted(ra),'b_note_refs':sorted(rb),
            'counterpart_explicitly_referenced_in_notes':counterpart_noted,
            'other_join_language_in_notes':join_language,
            'edition_transfer_language':transfer_language,
            'thickness_delta_cm_v3':td,
            'metadata_status_v3':metadata_status,
            'metadata_priority_score_v3':score,
        })

    audited.sort(key=lambda x:(x['metadata_status_v3']!='SURVIVES_METADATA_SCREEN',-x['metadata_priority_score_v3']))
    survivors=[x for x in audited if x['metadata_status_v3']=='SURVIVES_METADATA_SCREEN']
    result={
        'n_input':len(ranked),
        'n_survive_metadata_screen':len(survivors),
        'publication_state_counts':dict(Counter(x['publication_state_v3'] for x in audited)),
        'counterpart_note_reference_exclusions':sum(x['counterpart_explicitly_referenced_in_notes'] for x in audited),
        'top10_survivors':[{
            k:x.get(k) for k in ['a','b','sign_tfidf_similarity','publication_state_v3','thickness_delta_cm_v3','other_join_language_in_notes','edition_transfer_language','metadata_priority_score_v3']
        } for x in survivors[:10]],
        'correction_from_v2':'v2 mislabeled pairs with the same non-empty generic publication string as B_ONE_UNPUBLISHED. v3 classifies publication-field presence independently for each fragment and never treats equal non-empty strings as unpublished.',
        'claim_boundary':'This is still metadata/textual triage, not a physical join. A publishable join candidate needs tablet images or 3D/edge geometry, orientation/curvature compatibility, sign-by-sign philological continuity, and an explicit eBL/CDLI/British Museum/prior-literature audit. Publication-field absence is not proof of unpublished status.'
    }
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'audited_candidates_v3.json').write_text(json.dumps(audited,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'survivors_v3.json').write_text(json.dumps(survivors,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
