#!/usr/bin/env python3
from __future__ import annotations
import json, os
from pathlib import Path
import numpy as np
import fsspec, zarr, tifffile

BASE='vesuvius-challenge-open-data/PHerc0332/volumes/20231201141544-3.240um-70keV-masked.zarr'
TIFF='vesuvius-challenge-open-data/PHerc0139/segments/20250731185658-z_dbg_gen_09900/surface-volumes/9.362um-1.2m-113keV-volume-20250728140407.tifs/00.tif'

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/vesuvius_small_v2'));out.mkdir(parents=True,exist_ok=True)
    s3=fsspec.filesystem('s3',anon=True); res={'object_checks':{},'zarr':{},'surface_tiff':{}}
    # Some S3 public buckets permit GetObject but not ListBucket; directly request known metadata objects.
    for obj in [BASE+'/.zgroup',BASE+'/.zattrs',BASE+'/2/.zarray',TIFF]:
        try:
            info=s3.info('s3://'+obj);res['object_checks'][obj]={'exists':True,'size':info.get('Size',info.get('size'))}
        except Exception as e:res['object_checks'][obj]={'exists':False,'error':repr(e)}
    # Open level 2 directly, avoiding root listing/consolidation.
    mapper=s3.get_mapper('s3://'+BASE+'/2/',check=False,create=False)
    arr=zarr.open_array(mapper,mode='r')
    z=min(250,arr.shape[0]-1);y0=max(0,arr.shape[1]//2-128);x0=max(0,arr.shape[2]//2-128)
    crop=np.asarray(arr[z,y0:y0+256,x0:x0+256]);np.save(out/'volume_crop.npy',crop)
    res['zarr']={'shape':list(arr.shape),'chunks':list(arr.chunks),'dtype':str(arr.dtype),'crop':{'z':z,'y0':y0,'x0':x0,'shape':list(crop.shape),'mean':float(crop.mean()),'std':float(crop.std()),'min':float(crop.min()),'max':float(crop.max())}}
    # Read the one documented surface TIFF if reachable.
    try:
        with s3.open('s3://'+TIFF,'rb') as fh: im=tifffile.imread(fh)
        yy=max(0,im.shape[-2]//2-256);xx=max(0,im.shape[-1]//2-256);small=np.asarray(im[yy:yy+512,xx:xx+512]);np.save(out/'surface_crop.npy',small)
        res['surface_tiff']={'shape':list(im.shape),'dtype':str(im.dtype),'crop_shape':list(small.shape),'mean':float(small.mean()),'std':float(small.std())}
    except Exception as e:res['surface_tiff']={'error':repr(e)}
    res['claim_boundary']='Successful direct random access is only infrastructure validation, not ink detection or text recovery.'
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
