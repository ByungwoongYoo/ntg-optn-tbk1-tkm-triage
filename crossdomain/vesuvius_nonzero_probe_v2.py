#!/usr/bin/env python3
from __future__ import annotations
import json, os
from pathlib import Path
import numpy as np, fsspec, zarr
VOL='vesuvius-challenge-open-data/PHerc0332/volumes/20251211183505-2.399um-0.2m-78keV-masked.zarr'

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/vesuvius_v2'));out.mkdir(parents=True,exist_ok=True)
    fs=fsspec.filesystem('s3',anon=True);root=zarr.open(fs.get_mapper('s3://'+VOL+'/',check=False),mode='r')
    coarse=root['5']; candidates=[]
    # Frozen deterministic grid, chosen before seeing values.
    for z in np.linspace(64,coarse.shape[0]-65,12,dtype=int):
      for y in np.linspace(64,coarse.shape[1]-65,5,dtype=int):
       for x in np.linspace(64,coarse.shape[2]-65,5,dtype=int):
        a=np.asarray(coarse[z-16:z+16,y-32:y+32,x-32:x+32])
        candidates.append({'z':int(z),'y':int(y),'x':int(x),'mean':float(a.mean()),'std':float(a.std()),'nonzero':float(np.mean(a>0)),'max':int(a.max())})
    candidates.sort(key=lambda r:(r['std'],r['nonzero']),reverse=True);best=candidates[0]
    # Level 2 is 8x finer than level 5. Fetch a slice around the corresponding location.
    fine=root['2'];zf=min(fine.shape[0]-1,best['z']*8);yf=best['y']*8;xf=best['x']*8
    crop=np.asarray(fine[zf,max(0,yf-192):min(fine.shape[1],yf+192),max(0,xf-192):min(fine.shape[2],xf+192)])
    np.save(out/'nonzero_fine_slice.npy',crop)
    res={'volume':VOL,'grid_points_tested':len(candidates),'best_coarse_location':best,'fine_location':{'level':'2','z':int(zf),'y':int(yf),'x':int(xf)},'fine_crop':{'shape':list(crop.shape),'mean':float(crop.mean()),'std':float(crop.std()),'min':int(crop.min()),'max':int(crop.max()),'nonzero_fraction':float(np.mean(crop>0))},'top10_coarse':candidates[:10],
         'status':'NONZERO_CT_REGION_CONFIRMED' if crop.std()>0 and np.mean(crop>0)>0 else 'NO_NONZERO_REGION_FOUND',
         'claim_boundary':'This fixes the prior all-zero infrastructure crop by deterministic coarse-grid sampling. It confirms usable CT signal only. No surface segmentation, ink detection, text recovery, title identification, or discovery claim is made.'}
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
