"""Download compact raw DWR interval rasters with fixed windows and geometry."""
from pathlib import Path
import argparse,json,hashlib
from concurrent.futures import ThreadPoolExecutor
import requests
import numpy as np
import rasterio
from dwr_units import WINDOWS,duration_years

BASE='https://gis.water.ca.gov/arcgisimg/rest/services/SAR/'
ANNUAL=BASE+'Vertical_Displacement_TRE_ALTAMIRA_Annual_Rate_Mosaic/ImageServer'
CUMULATIVE=BASE+'Vertical_Displacement_TRE_ALTAMIRA_Total_Since_20150613_Mosaic/ImageServer'
BBOX='-119.90046966824156,36.59977083241381,-119.6998252715033,36.85082062245584'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    catalogs={}
    for kind,url in [('annual',ANNUAL),('cumulative',CUMULATIVE)]:
        r=requests.get(url+'/query',params={'f':'json','where':'1=1','outFields':'*','returnGeometry':'false'},timeout=60);r.raise_for_status()
        catalogs[kind]=r.json()
        if 'features' not in catalogs[kind]:raise ValueError(catalogs[kind])
        (out/f'{kind}_catalog.json').write_text(json.dumps(catalogs[kind],indent=2),encoding='utf-8')
        metadata=requests.get(url,params={'f':'pjson'},timeout=60).json()
        (out/f'{kind}_service_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    def download(item):
        name,(start,end)=item;kind='cumulative' if name.startswith('total_') else 'annual'
        date_suffix=start.replace('-','')+'_'+end.replace('-','')
        features=[f['attributes'] for f in catalogs[kind]['features'] if f['attributes']['Name'].endswith(date_suffix)]
        if len(features)!=1:raise RuntimeError((name,'ambiguous catalog match',features))
        feature=features[0];url=ANNUAL if kind=='annual' else CUMULATIVE
        params={'f':'image','bbox':BBOX,'bboxSR':4326,'imageSR':4326,'size':'177,279',
            'format':'tiff','pixelType':'F32','noData':'-9999','interpolation':'RSP_NearestNeighbor',
            'adjustAspectRatio':'false','renderingRule':json.dumps({'rasterFunction':'None'}),
            'mosaicRule':json.dumps({'mosaicMethod':'esriMosaicLockRaster','lockRasterIds':[feature['OBJECTID']],'mosaicOperation':'MT_FIRST'})}
        data=requests.get(url+'/exportImage',params=params,timeout=90);data.raise_for_status()
        path=out/name;path.write_bytes(data.content)
        with rasterio.open(path) as src:
            values=src.read(1);valid=np.isfinite(values)&(values!=src.nodata)
            assert src.shape==(279,177) and str(src.crs)=='EPSG:4326'
            detail={'shape':list(src.shape),'transform':list(src.transform)[:6],'valid_cells':int(valid.sum()),
                'minimum_feet':float(values[valid].min()),'maximum_feet':float(values[valid].max())}
        return {'file':name,'sha256':hashlib.sha256(data.content).hexdigest(),'source_service':url,
            'source_item':feature,'request':params,'raster':detail,'interval_years':duration_years(name),
            'conversion':'mm/year = raw feet * 304.8 / ((DateTo-DateFrom in days)/365.2425)'}
    with ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(download,WINDOWS.items()))
    (out/'raw_input_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print(json.dumps([{'file':r['file'],'item':r['source_item']['OBJECTID'],'sha256':r['sha256']} for r in records],indent=2))

if __name__=='__main__':main()
