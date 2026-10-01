"""Freeze provider implementation evidence and exact station source windows."""
from pathlib import Path
import hashlib
import json
import zipfile
import sys
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from pyproj import Geod
from source_inventory import RemoteZip
from independent_gnss_transfer import REGIONS

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / 'inputs'
OUT = IN / 'gnss_reference'


def main():
    meta = json.loads((IN/'provider_metadata/zenodo_19493073.json').read_text())
    f = next(x for x in meta['files'] if x['key']=='software.zip')
    names = ['software/MintPy/src/mintpy/utils/time_func.py',
             'software/MintPy/LICENSE']
    archived = {}
    with zipfile.ZipFile(RemoteZip(f['links']['self'], f['size'])) as z:
        for name in names:
            if name not in z.namelist():
                matches = [n for n in z.namelist() if n.endswith('/'+Path(name).name) and '/MintPy/' in n]
                if len(matches)!=1: raise RuntimeError(matches)
                name = matches[0]
            content=z.read(name)
            dest=OUT/('provider_time_func.py' if name.endswith('.py') else 'provider_MintPy_LICENSE.txt')
            dest.write_bytes(content)
            archived[name]={'sha256':hashlib.sha256(content).hexdigest(),'bytes':len(content)}
    code=(OUT/'provider_time_func.py').read_text()
    ix=code.index('def get_design_matrix4time_func')
    section=code[ix:ix+6500]
    (OUT/'provider_time_function_default_excerpt.txt').write_bytes(section.encode())
    audit={'archive':f, 'members':archived,
           'gnss_notebook_model_argument':'Not supplied: get_los_obs(...), model=None',
           'gnss_default':'See retained source excerpt: polynomial degree 1; no periodic or step terms.',
           'insar_model':'Provider paper: linear trend, annual/semiannual terms and earthquake steps.',
           'claim_limit':'Periods and spatial LOS frame are matched, but temporal models are not identical. Archived GNSS screening and fits are retained; this is not a fresh raw-daily refit.'}
    (OUT/'provider_default_model_audit.json').write_bytes((json.dumps(audit,indent=2)+'\n').encode())
    geod=Geod(ellps='WGS84')
    for region,track,bbox,reference,local_cal,start,end in REGIONS:
        data=pd.read_csv(OUT/f'{track}_archived_rates_with_coordinates.csv')
        roi=data.longitude.between(bbox[0],bbox[2])&data.latitude.between(bbox[1],bbox[3])
        data=data[roi|data.station.eq(reference)].copy().sort_values('station').reset_index(drop=True)
        path=IN/'aria_california'/f'{track}_velocity_ref{reference}.tif'
        windows=[]; metadata=[]
        with rasterio.open(path) as ds:
            for _,s in data.iterrows():
                row,col=ds.index(s.longitude,s.latitude)
                assert 3<=row<ds.height-3 and 3<=col<ds.width-3
                windows.append(ds.read(1,window=Window(col-3,row-3,7,7)))
                xx,yy=ds.xy(row,col)
                metadata.append({'station':s.station,'sample_row':row,'sample_col':col,
                                 'pixel_center_longitude':xx,'pixel_center_latitude':yy})
            protocol={'region':region,'track':track,'source_name':path.name,
                      'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                      'shape':[ds.height,ds.width],'transform':list(ds.transform)[:6],
                      'crs':str(ds.crs),'unit':'m/yr','window_radius':3,
                      'selection':'All archived stations inside fixed ROI plus provider reference; no outcome screening.',
                      'purpose':'Lossless exact source windows for offline comparison; full selected GeoTIFFs supplied separately.'}
        np.savez_compressed(OUT/f'{track}_station_windows.npz',stations=data.station.to_numpy(dtype=str),
                            source_windows=np.asarray(windows,dtype=np.float32))
        pd.DataFrame(metadata).to_csv(OUT/f'{track}_station_window_metadata.csv',index=False)
        (OUT/f'{track}_station_window_provenance.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
        print(track, 'frozen source windows',len(data),flush=True)
    print(section[:2200],flush=True)

if __name__=='__main__': main()
