#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,io,re,time
from pathlib import Path
import requests
from PIL import Image,ImageOps,ImageDraw,ImageFont

def cdli_num(rec,side):
    x=str(rec.get(side+'_cdli') or '');m=re.search(r'(\d+)',x);return m.group(1) if m else ''
def get_photo(s,num):
    urls=[f'https://cdli.earth/dl/photo/P{num}.jpg',f'https://cdli.earth/dl/tn_photo/P{num}.jpg',f'https://cdli.earth/dl/lineart/P{num}_l.jpg',f'https://cdli.earth/dl/tn_lineart/P{num}_l.jpg']
    errors=[]
    for u in urls:
        try:
            r=s.get(u,timeout=60)
            if r.ok and r.headers.get('content-type','').lower().startswith('image') and len(r.content)>1000:
                return Image.open(io.BytesIO(r.content)).convert('RGB'),u,len(r.content)
            errors.append((u,r.status_code,len(r.content)))
        except Exception as e:errors.append((u,repr(e)))
    return None,None,errors
def fit(im,maxh=550,maxw=700):
    im=ImageOps.exif_transpose(im); im.thumbnail((maxw,maxh));return im
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--candidates',required=True);ap.add_argument('--out-dir',required=True);a=ap.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True);s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 CDLI-medical-join-photo-audit/1.0'})
    cand=json.loads(Path(a.candidates).read_text(encoding='utf-8'))[:20];manifest=[];thumbs=[]
    for rank,r in enumerate(cand,1):
        panels=[];meta={'rank':rank,'a':r['a'],'b':r['b'],'model_score':r.get('physical_join_model_score'),'sign_similarity':r.get('sign_tfidf_similarity'),'downloads':{}}
        for side in ('a','b'):
            num=cdli_num(r,side);im,url,info=get_photo(s,num) if num else (None,None,'no CDLI number');meta['downloads'][side]={'cdli_number':num,'url':url,'info':info}
            if im is None:
                im=Image.new('RGB',(450,300),'white');ImageDraw.Draw(im).text((15,130),f'No public image\n{r[side]}',fill='black')
            panels.append(fit(im));time.sleep(.08)
        h=max(x.height for x in panels)+80;w=sum(x.width for x in panels)+30
        canvas=Image.new('RGB',(w,h),'white');d=ImageDraw.Draw(canvas);x=10
        d.text((10,10),f'#{rank} {r["a"]}  <->  {r["b"]}   model={r.get("physical_join_model_score",0):.3f} text={r.get("sign_tfidf_similarity",0):.3f}',fill='black')
        for im in panels:
            canvas.paste(im,(x,55));x+=im.width+10
        fp=out/f'pair_{rank:02d}_{re.sub("[^A-Za-z0-9]+","_",r["a"])}__{re.sub("[^A-Za-z0-9]+","_",r["b"])}.jpg';canvas.save(fp,quality=92);meta['packet_file']=fp.name;manifest.append(meta)
        t=canvas.copy();t.thumbnail((800,400));thumbs.append(t)
    # contact sheet
    if thumbs:
        W=max(i.width for i in thumbs); H=sum(i.height+10 for i in thumbs)
        sheet=Image.new('RGB',(W,H),'white');y=0
        for im in thumbs:sheet.paste(im,(0,y));y+=im.height+10
        sheet.save(out/'CONTACT_SHEET_TOP20.jpg',quality=90)
    (out/'PHOTO_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'pairs':len(manifest),'with_two_images':sum(bool(m['downloads']['a']['url'] and m['downloads']['b']['url']) for m in manifest),'contact_sheet':bool(thumbs),'claim_boundary':'Photo packets are for visual falsification only. Front/line-art images alone cannot establish a physical join; edge/surface fit or curator inspection is required.'},indent=2))
if __name__=='__main__':main()
