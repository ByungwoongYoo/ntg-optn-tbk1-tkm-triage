#!/usr/bin/env python3
from __future__ import annotations
import json, os, subprocess
from pathlib import Path
import numpy as np
import fsspec, zarr

VOL='vesuvius-challenge-open-data/PHerc0332/volumes/20251211183505-2.399um-0.2m-78keV-masked.zarr'

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/vesuvius_current'));out.mkdir(parents=True,exist_ok=True)
    s3=fsspec.filesystem('s3',anon=True)
    # Root listing is known to work now; capture current levels and small crop.
    store=s3.get_mapper('s3://'+VOL+'/',check=False)
    root=zarr.open(store,mode='r')
    keys=list(root.array_keys())
    res={'volume':VOL,'levels':{}}
    for k in keys:
        a=root[k];res['levels'][k]={'shape':list(a.shape),'chunks':list(a.chunks),'dtype':str(a.dtype)}
    level='2' if '2' in keys else keys[-1]
    a=root[level];z=min(250,a.shape[0]-1);y=max(0,a.shape[1]//2-128);x=max(0,a.shape[2]//2-128)
    crop=np.asarray(a[z,y:y+256,x:x+256]);np.save(out/'current_volume_crop.npy',crop)
    res['crop']={'level':level,'z':z,'y0':y,'x0':x,'shape':list(crop.shape),'min':float(crop.min()),'max':float(crop.max()),'mean':float(crop.mean()),'std':float(crop.std())}
    res['claim_boundary']='Current public data access is confirmed. This is infrastructure validation only, not ink detection or unread-text recovery.'
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
