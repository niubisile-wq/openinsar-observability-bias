from pathlib import Path
import json
import hashlib
from datetime import datetime
import urllib.parse
import numpy as np
import rasterio
from source_inventory import fetch

ROOT=Path(__file__).resolve().parent.parent
STUDY=ROOT.parent/'01_当前返修工作区/InSAR_revision_20260926/strengthening'
OUT=ROOT/'inputs/fresno_real_blocks'
OUT.mkdir(parents=True,exist_ok=True)
SERVICE='https://gis.water.ca.gov/arcgisimg/rest/services/SAR/Vertical_Displacement_TRE_ALTAMIRA_Annual_Rate_Mosaic/ImageServer'
features=json.loads((STUDY/'external/dwr/roi_features.json').read_text(encoding='utf-8'))['features']
points=np.array([p for f in features for ring in f['geometry']['rings'] for p in ring])
bbox=[points[:,0].min(),points[:,1].min(),points[:,0].max(),points[:,1].max()]
sample=np.array(features[0]['geometry']['rings'][0]);dx=sample[:,0].max()-sample[:,0].min();dy=sample[:,1].max()-sample[:,1].min()
width=int(round((bbox[2]-bbox[0])/dx));height=int(round((bbox[3]-bbox[1])/dy))
items=json.loads((ROOT/'inputs/provider_metadata/dwr_annual_items.json').read_text(encoding='utf-8'))
matches=[f['attributes'] for f in items['features'] if f['attributes']['Name'].endswith('20201001_20211001')]
if len(matches)!=1: raise RuntimeError('Expected one unambiguous annual image, found '+repr(matches))
item=matches[0]
query={'f':'image','bbox':','.join(map(str,bbox)),'bboxSR':'4326','imageSR':'4326','size':f'{width},{height}',
       'format':'tiff','pixelType':'F32','noData':'-9999','interpolation':'RSP_NearestNeighbor','adjustAspectRatio':'false',
       'renderingRule':json.dumps({'rasterFunction':'None'}),
       'mosaicRule':json.dumps({'mosaicMethod':'esriMosaicLockRaster','lockRasterIds':[item['OBJECTID']],
                              'mosaicOperation':'MT_FIRST'})}
raw,_,_=fetch(SERVICE+'/exportImage?'+urllib.parse.urlencode(query))
target=OUT/'DWR_20201001_20211001.tif';target.write_bytes(raw)
with rasterio.open(target) as ds:
    info={'shape':[ds.height,ds.width],'crs':str(ds.crs),'transform':list(ds.transform)[:6],
          'nodata':ds.nodata,'dtype':ds.dtypes[0]}
    arr=ds.read(1,masked=True)
    assert (ds.width,ds.height)==(width,height) and ds.dtypes[0]=='float32'
    info.update(valid_pixels=int((~np.ma.getmaskarray(arr)).sum()),minimum=float(arr.min()),maximum=float(arr.max()))
metadata={'source_service':SERVICE,'provider_image':item,'request':query,'raster':info,
          'sha256':hashlib.sha256(raw).hexdigest(),'units':'Provider feet over the stated annual period; convert to millimetres/year using exact period duration.',
          'duration_years':(item['DateTo']-item['DateFrom'])/1000/86400/365.2425,
          'footprint_source':'Frozen 38,738 DWR observation-cell polygons, cross-vintage support benchmark',
          'purpose':'Controlled masking on known published DWR values; not independent validation of InSAR or original missing cells.'}
(OUT/'DWR_rate_provenance.json').write_bytes((json.dumps(metadata,indent=2)+'\n').encode())
print(json.dumps(metadata['raster']), 'feature_count',len(features),'duration',metadata['duration_years'],flush=True)
