#!/usr/bin/env python3
from __future__ import annotations
import io,json,os,re,time
from pathlib import Path
import requests
from PIL import Image,ImageOps,ImageDraw
MANUSCRIPTS={
 'fr3019':('https://gallica.bnf.fr/iiif/ark:/12148/btv1b9059994n/manifest.json',range(142,158)),
 'fr20506':('https://gallica.bnf.fr/iiif/ark:/12148/btv1b525047581/manifest.json',range(265,285)),
}
def main():
 out=Path(os.environ.get('OUT_DIR','artifact/ranzo_parallel'));out.mkdir(parents=True,exist_ok=True);s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 Ranzo-cipher-IIIF-research/2.0'})
 result={}
 for name,(url,idxs) in MANUSCRIPTS.items():
  r=s.get(url,timeout=90);r.raise_for_status();man=r.json();canv=(man.get('sequences') or [{}])[0].get('canvases') or [];sub=out/name;sub.mkdir(exist_ok=True);rows=[];thumbs=[]
  for idx in idxs:
   if idx<1 or idx>len(canv):continue
   c=canv[idx-1]; imurl=None
   try:imurl=c['images'][0]['resource'].get('@id')
   except Exception:pass
   if imurl:imurl=re.sub(r'/full/full/0/native\.jpg$',r'/full/1400,/0/native.jpg',imurl)
   rec={'view_index_1based':idx,'label':c.get('label'),'canvas_id':c.get('@id'),'image_url':imurl}
   if imurl:
    try:
     q=s.get(imurl,timeout=90);rec['status']=q.status_code;rec['bytes']=len(q.content)
     if q.ok and q.headers.get('content-type','').startswith('image'):
      im=Image.open(io.BytesIO(q.content)).convert('RGB');im=ImageOps.exif_transpose(im);im.thumbnail((500,700));fp=sub/f'view_{idx:03d}.jpg';im.save(fp,quality=90);rec['file']=fp.name
      tile=Image.new('RGB',(520,745),'white');tile.paste(im,((520-im.width)//2,35));ImageDraw.Draw(tile).text((8,8),f'{name} view {idx}',fill='black');thumbs.append(tile)
    except Exception as e:rec['error']=repr(e)
   rows.append(rec);time.sleep(.04)
  if thumbs:
   cols=4;tw,th=520,745;nr=(len(thumbs)+cols-1)//cols;sheet=Image.new('RGB',(cols*tw,nr*th),'white')
   for i,t in enumerate(thumbs):sheet.paste(t,((i%cols)*tw,(i//cols)*th))
   sheet.save(sub/'CONTACT_SHEET.jpg',quality=88)
  result[name]={'manifest':url,'canvas_count':len(canv),'rows':rows}
 (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:{'canvas_count':v['canvas_count'],'views':[x['view_index_1based'] for x in v['rows']]} for k,v in result.items()},indent=2))
if __name__=='__main__':main()
