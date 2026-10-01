from pathlib import Path
import json
import urllib.parse
from source_inventory import fetch
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'inputs/provider_metadata'
queries={
 'fresno_arcgis_search': 'https://www.arcgis.com/sharing/rest/search?'+urllib.parse.urlencode({'q':'Fresno (title:Parcels OR title:"Land Use")','num':100,'f':'json'}),
 'fresno_webmap': 'https://www.arcgis.com/sharing/rest/content/items/d7be2dd59704497ab6440db8f76ca71c/data?f=json',
 'dwr_annual': 'https://gis.water.ca.gov/arcgisimg/rest/services/SAR/Vertical_Displacement_TRE_ALTAMIRA_Annual_Rate_Mosaic/ImageServer?f=pjson',
 'dwr_annual_items':'https://gis.water.ca.gov/arcgisimg/rest/services/SAR/Vertical_Displacement_TRE_ALTAMIRA_Annual_Rate_Mosaic/ImageServer/query?'+urllib.parse.urlencode({'where':'1=1','outFields':'*','returnGeometry':'false','f':'json'}),
}
for name,url in queries.items():
    try:
        raw,_,_=fetch(url);(OUT/(name+'.json')).write_bytes(raw);obj=json.loads(raw)
        if name=='fresno_arcgis_search': print([(r['title'],r['owner'],r.get('url'),r['id']) for r in obj.get('results',[])],flush=True)
        elif name=='fresno_webmap': print(obj.get('operationalLayers'),flush=True)
        elif name=='dwr_annual': print('DWRpixel',obj.get('pixelSizeX'),obj.get('pixelSizeY'),'fields',[(f['name'],f['type']) for f in obj.get('fields',[])],flush=True)
        elif name=='dwr_annual_items': print([f['attributes'] for f in obj.get('features',[])][:8],flush=True)
    except Exception as e: print(name,repr(e),flush=True)
