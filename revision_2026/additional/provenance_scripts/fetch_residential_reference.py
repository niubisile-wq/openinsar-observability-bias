"""Download official existing-land-use geometry without owner/address fields."""
from pathlib import Path
import json
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from source_inventory import fetch

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'inputs/residential_reference'
OUT.mkdir(parents=True,exist_ok=True)
SOURCES=[
    'https://services2.arcgis.com/WkBUojyNPhsWOk1W/arcgis/rest/services/Existing_Land_Use/FeatureServer/19',
    'https://services2.arcgis.com/FS7YxIXpWoaR2sAe/ArcGIS/rest/services/Existing_Land_Use/FeatureServer/0',
    'https://gis4u.fresno.gov/arcgis/rest/services/PublicInfoServices/AddrParcelStreet/FeatureServer/1',
    'https://gis4u.fresno.gov/exa/rest/services/PublicInfoServices/AddrParcelStreet/FeatureServer/1',
]

def query(url,params):
    encoded=urllib.parse.urlencode(params)
    if len(encoded)>1700:
        request=urllib.request.Request(url,data=encoded.encode(),headers={'User-Agent':'InSAR-revision-reproducibility/1.0','Content-Type':'application/x-www-form-urlencoded'})
        with urllib.request.urlopen(request,timeout=120) as r: raw=r.read()
    else: raw,_,_=fetch(url+'?'+encoded)
    obj=json.loads(raw)
    if 'error' in obj: raise RuntimeError(json.dumps(obj['error']))
    return obj,raw

def main():
    attempts=[]
    for source in SOURCES:
        try:
            meta,raw=query(source,{'f':'pjson'})
            (OUT/'layer_metadata.json').write_bytes(raw)
            print('SOURCE',source,'TYPE',meta.get('geometryType'),'FIELDS', [v['name'] for v in meta.get('fields',[])],flush=True)
            if meta.get('geometryType')!='esriGeometryPolygon': raise RuntimeError('Need polygon land-use geometry')
            fields=[v['name'] for v in meta['fields']]
            oid=meta.get('objectIdField') or next(v['name'] for v in meta['fields'] if v['type']=='esriFieldTypeOID')
            selected=[n for n in fields if n.upper() in ['EXISTING_LAND_USE','EXISTING_LAND_USE_TEXT','LAND_USE','LANDUSE','DESCRIPTION','ELU','ELUTEXT','ELU_DESC','EXLU','EXLU_DESC','LANDUSE_DESC','LANDUSE_DESCRIPTION','LAYERREFRESHDATE']]
            if not selected: raise RuntimeError('No recognized existing-land-use classification fields')
            common={'where':'1=1','geometry':'-119.9,36.6,-119.7,36.85','geometryType':'esriGeometryEnvelope',
                    'inSR':'4326','spatialRel':'esriSpatialRelIntersects','f':'json'}
            ids,_=query(source+'/query',{**common,'returnIdsOnly':'true'})
            object_ids=sorted(ids['objectIds']);all_features=[]
            if not object_ids: raise RuntimeError('No features in the fixed Fresno ROI; reject incompatible source')
            print('Matched features',len(object_ids),flush=True)
            def retrieve(offset):
                chunk=object_ids[offset:offset+500]
                target=OUT/f'landuse_{offset:06d}.geojson'
                if target.exists():
                    raw=target.read_bytes();obj=json.loads(raw)
                    fields_present=set(obj['features'][0]['properties']) if obj.get('features') else set()
                    if not set([oid]+selected)<=fields_present: obj=None
                else: obj=None
                if obj is None:
                    obj,raw=query(source+'/query',{'objectIds':','.join(map(str,chunk)),
                                                 'outFields':','.join([oid]+selected),'outSR':'4326',
                                                 'returnGeometry':'true','f':'geojson'})
                got=obj.get('features',[])
                if len(got)!=len(chunk): raise RuntimeError('Truncated parcel query')
                return offset,raw,got
            with ThreadPoolExecutor(max_workers=4) as pool:
                futures=[pool.submit(retrieve,offset) for offset in range(0,len(object_ids),500)]
                for future in as_completed(futures):
                    offset,raw,got=future.result()
                    (OUT/f'landuse_{offset:06d}.geojson').write_bytes(raw)
                    all_features.extend(got)
                    if len(all_features)%5000==0: print('Downloaded',len(all_features),flush=True)
            recovered=[f['properties'][oid] for f in all_features]
            assert sorted(recovered)==object_ids and len(set(recovered))==len(object_ids)
            categories={}
            for f in all_features:
                key=json.dumps({k:f['properties'].get(k) for k in selected if k.upper()!='LAYERREFRESHDATE'},sort_keys=True)
                categories[key]=categories.get(key,0)+1
            metadata={'source_url':source,'retrieved_date':'2026-10-01','feature_count':len(all_features),
                      'selected_fields':[oid]+selected,'categories':categories,
                      'source_independence':'Municipal existing-land-use map, not an InSAR or population model. Unknown map epoch; do not call household truth.',
                      'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob('landuse_*.geojson'))},
                      'checks':{'all_provider_roi_ids_retrieved_once':True},'attempts':attempts}
            (OUT/'provenance.json').write_bytes((json.dumps(metadata,indent=2)+'\n').encode())
            print('Categories',json.dumps(categories,indent=2),flush=True)
            return
        except Exception as e:
            attempts.append({'source':source,'error':repr(e)})
            print('Attempt error',repr(e),flush=True)
    (OUT/'access_errors.json').write_bytes((json.dumps(attempts,indent=2)+'\n').encode())
    raise RuntimeError('No official residential polygons recovered')

if __name__=='__main__': main()
