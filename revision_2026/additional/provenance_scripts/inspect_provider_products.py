from pathlib import Path
import json
import shapefile
import urllib.parse
from source_inventory import fetch
ROOT=Path(__file__).resolve().parent.parent
SRC=ROOT/'inputs/aria_california'
OUT=ROOT/'inputs/provider_metadata'
for name in ['A064','A137']:
    notebook=SRC/f'run_{name}_completed.ipynb'
    if notebook.exists():
        obj=json.loads(notebook.read_text(encoding='utf-8'));blocks=[]
        for i,cell in enumerate(obj['cells']):
            blocks.append(f'CELL {i} ({cell["cell_type"]})\n'+''.join(cell.get('source',[])))
            for output in cell.get('outputs',[]):
                if 'text' in output: blocks.append(''.join(output['text']))
                txt=output.get('data',{}).get('text/plain')
                if txt: blocks.append(''.join(txt))
        (OUT/f'{name}_notebook_text.txt').write_bytes(('\n\n'.join(blocks)+'\n').encode())
    path=SRC/f'track_{name[1:]}_CA.shp'
    if path.exists():
        try:
            if path.read_bytes().lstrip().startswith(b'{'):
                print(name,'AOI GeoJSON named .shp',path.read_text(encoding='utf-8')[:1200],flush=True)
            else:
                reader=shapefile.Reader(shp=str(path),dbf=None)
                print(name,'AOI bounds',reader.bbox,'shapes',len(reader.shapes()),flush=True)
        except Exception as e: print(name,'AOI error',repr(e),flush=True)
url='https://services2.arcgis.com/WkBUojyNPhsWOk1W/arcgis/rest/services/Existing_Land_Use/FeatureServer?f=pjson'
raw,_,_=fetch(url);(OUT/'fresno_existing_landuse_service.json').write_bytes(raw)
print('Fresno official layers',json.loads(raw).get('layers'),flush=True)
url='https://www.arcgis.com/sharing/rest/content/items/bc57c53874e14b1e850d2324b7a568ca?f=pjson'
raw,_,_=fetch(url);(OUT/'fresno_existing_landuse_item.json').write_bytes(raw)
print('Fresno item', {k:json.loads(raw).get(k) for k in ['title','owner','created','modified','description','licenseInfo','accessInformation']},flush=True)
