"""Download actual DWR published observation-cell polygons in a fixed ROI."""
from pathlib import Path
import json,time,zipfile,requests
from common import ROOT,dump,sha256
from fetch_external import fetch

DEST=ROOT/'external/dwr';DEST.mkdir(exist_ok=True)
BASE='https://gis.water.ca.gov/arcgis/rest/services/Elevation/Vertical_Displacement_Point_Data_Locations_2026Q1/MapServer/0'
ROI=[-119.9,36.6,-119.7,36.85]

def query(params,path):
    if path.exists():return json.loads(path.read_text())
    for attempt in range(4):
        try:
            r=requests.get(BASE+'/query',params={'f':'json',**params},timeout=100);r.raise_for_status();j=r.json()
            if 'error' in j:raise ValueError(j['error'])
            dump(path,j);return j
        except Exception:
            if attempt==3:raise
    raise RuntimeError('query')

def main():
    ids=query(dict(where='1=1',geometry=','.join(map(str,ROI)),geometryType='esriGeometryEnvelope',inSR=4326,spatialRel='esriSpatialRelIntersects',returnIdsOnly='true'),DEST/'roi_ids.json')['objectIds']
    ids=sorted(ids);features=[]
    for i in range(0,len(ids),800):
        selected=ids[i:i+800]
        j=query(dict(objectIds=','.join(map(str,selected)),outFields='OBJECTID,CODE,Longitude,Latitude,DataTable,StartDate',returnGeometry='true',outSR=4326),DEST/f'features_{i:06d}.json')
        batch=j.get('features',[])
        if {f['attributes']['OBJECTID'] for f in batch}!=set(selected):raise ValueError('Incomplete feature batch')
        features.extend(batch);print('DWR',len(features),'/',len(ids),flush=True)
    dump(DEST/'roi_features.json',dict(geometryType='esriGeometryPolygon',spatialReference=dict(wkid=4326),features=features))
    name='GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0_R6_C7'
    url='https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_4326_3ss/V1-0/tiles/'+name+'.zip'
    pop=fetch(url,DEST/(name+'.zip'))
    with zipfile.ZipFile(DEST/(name+'.zip')) as archive:
        for member in archive.namelist():
            if member.endswith('.tif'):
                target=DEST/Path(member).name
                if not target.exists():
                    with archive.open(member) as src,target.open('wb') as out:
                        import shutil
                        shutil.copyfileobj(src,out)
    dump(DEST/'manifest.json',dict(layer_url=BASE,query_roi=ROI,features=len(features),features_sha256=sha256(DEST/'roi_features.json'),population_download=pop,interpretation='Published polygons of 100 m measurement support, not buffers invented around point coordinates; no velocity is inferred in uncovered areas'))
    print('DWR DOWNLOAD COMPLETE',len(features),flush=True)

if __name__=='__main__':main()
