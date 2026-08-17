#!/usr/bin/env python3
from __future__ import annotations
import io,json,os,re,time
from pathlib import Path
import requests
from PIL import Image,ImageOps,ImageDraw

MAN='https://gallica.bnf.fr/iiif/ark:/12148/btv1b9059908w/manifest.json'

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/bnf_fr2988'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 historical-cipher-IIIF-research/1.0'})
    r=s.get(MAN,timeout=90);r.raise_for_status();man=r.json();(out/'manifest.json').write_text(json.dumps(man,ensure_ascii=False,indent=2),encoding='utf-8')
    canv=(man.get('sequences') or [{}])[0].get('canvases') or []
    rows=[];thumbs=[]
    for i,c in enumerate(canv[:35],1):
        label=c.get('label'); imurl=None
        try:
            res=c['images'][0]['resource']; imurl=res.get('@id')
            if imurl: imurl=re.sub(r'/full/full/0/native\.jpg$',r'/full/1200,/0/native.jpg',imurl)
        except Exception: pass
        rec={'index':i,'label':label,'canvas_id':c.get('@id'),'image_url':imurl}
        if imurl:
            try:
                q=s.get(imurl,timeout=90);rec['status']=q.status_code;rec['bytes']=len(q.content)
                if q.ok and q.headers.get('content-type','').startswith('image'):
                    img=Image.open(io.BytesIO(q.content)).convert('RGB');img=ImageOps.exif_transpose(img);img.thumbnail((450,600));
                    fp=out/f'view_{i:03d}.jpg';img.save(fp,quality=88);rec['file']=fp.name
                    tile=Image.new('RGB',(470,650),'white');tile.paste(img,((470-img.width)//2,35));ImageDraw.Draw(tile).text((8,8),f'#{i} label={label}',fill='black');thumbs.append(tile)
            except Exception as e:rec['error']=repr(e)
        rows.append(rec);time.sleep(.05)
    if thumbs:
        cols=4;tw,th=470,650;nr=(len(thumbs)+cols-1)//cols;sheet=Image.new('RGB',(cols*tw,nr*th),'white')
        for i,t in enumerate(thumbs):sheet.paste(t,((i%cols)*tw,(i//cols)*th))
        sheet.save(out/'CONTACT_SHEET_FIRST35.jpg',quality=88)
    result={'status':'BNF_FR2988_IIIF_ACQUIRED','manifest_label':man.get('label'),'canvas_count':len(canv),'first35':rows,
            'metadata':man.get('metadata'),'claim_boundary':'This step only acquires public BnF IIIF images and audits folio/page alignment for DECRYPT record 2294. No decryption claim is made.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'status':result['status'],'canvas_count':len(canv),'first_labels':[x['label'] for x in rows]},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
