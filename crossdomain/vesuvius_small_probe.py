#!/usr/bin/env python3
from __future__ import annotations
import json, os
from pathlib import Path
import numpy as np
import fsspec, zarr, tifffile

ZARR='vesuvius-challenge-open-data/PHerc0332/volumes/20231201141544-3.240um-70keV-masked.zarr/'
TIFF='vesuvius-challenge-open-data/PHerc0139/segments/20250731185658-z_dbg_gen_09900/surface-volumes/9.362um-1.2m-113keV-volume-20250728140407.tifs/00.tif'

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/vesuvius_small'));out.mkdir(parents=True,exist_ok=True)
    s3=fsspec.filesystem('s3',anon=True)
    res={'zarr':{},'surface_tiff':{}}
    # Zarr metadata and one small multiresolution crop only.
    store=s3.get_mapper('s3://'+ZARR)
    root=zarr.open(store,mode='r')
    res['zarr']['keys']=list(root.array_keys())
    for key in res['zarr']['keys'][:6]:
        a=root[key];res['zarr'].setdefault('levels',{})[key]={'shape':list(a.shape),'chunks':list(a.chunks),'dtype':str(a.dtype)}
    level='2' if '2' in root else res['zarr']['keys'][-1]
    a=root[level]; z=min(250,a.shape[0]-1); y0=max(0,a.shape[1]//2-128);x0=max(0,a.shape[2]//2-128)
    crop=np.asarray(a[z,y0:y0+256,x0:x0+256])
    np.save(out/'volume_crop.npy',crop)
    res['zarr']['crop']={'level':level,'z':z,'y0':y0,'x0':x0,'shape':list(crop.shape),'min':float(crop.min()),'max':float(crop.max()),'mean':float(crop.mean()),'std':float(crop.std())}
    # Surface TIFF metadata. Read only one public layer; tifffile may issue ranged reads through fsspec.
    info=s3.info('s3://'+TIFF);res['surface_tiff']['object_bytes']=info.get('Size',info.get('size'))
    with s3.open('s3://'+TIFF,'rb') as fh:
        im=tifffile.imread(fh)
    res['surface_tiff'].update({'shape':list(im.shape),'dtype':str(im.dtype),'min':float(np.min(im)),'max':float(np.max(im)),'mean':float(np.mean(im)),'std':float(np.std(im))})
    # Save a 512x512 center crop rather than full source layer.
    yy=max(0,im.shape[-2]//2-256);xx=max(0,im.shape[-1]//2-256);small=np.asarray(im[yy:yy+512,xx:xx+512])
    np.save(out/'surface_crop.npy',small)
    res['claim_boundary']='This proves small random-access feasibility only. It is not ink detection and contains no new readable-text claim. A survivor must reproduce a labeled known-ink region with an official/public baseline before any new-scroll search.'
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
    print(json.dumps(res,indent=2))
if __name__=='__main__':main()
