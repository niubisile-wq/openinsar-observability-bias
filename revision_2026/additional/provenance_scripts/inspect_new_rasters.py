from pathlib import Path
import json
import rasterio
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'inputs/provider_metadata'
for path in (ROOT/'inputs/aria_california').glob('*.tif'):
    with rasterio.open(path) as ds:
        metadata={'file':path.name,'shape':[ds.height,ds.width],'bounds':list(ds.bounds),'crs':str(ds.crs),
                  'transform':list(ds.transform)[:6],'nodata':ds.nodata,'tags':ds.tags(),'namespaces':ds.tag_namespaces(),
                  'band_tags':ds.tags(1)}
    (OUT/(path.stem+'_metadata.json')).write_bytes((json.dumps(metadata,indent=2)+'\n').encode())
    print(json.dumps(metadata),flush=True)
for name in ['A064','A137']:
    path=ROOT/'inputs/aria_california'/f'run_{name}_completed.ipynb'
    obj=json.loads(path.read_text(encoding='utf-8'))
    code='\n\n'.join('CELL '+str(i)+'\n'+''.join(c.get('source',[])) for i,c in enumerate(obj['cells']) if c['cell_type']=='code')
    (OUT/f'{name}_provider_code.py.txt').write_bytes(code.encode())
    excerpts=[]
    for i,c in enumerate(obj['cells']):
        source=''.join(c.get('source',[]))
        if ('print' in source and any(k in source for k in ['metadata','sites_df','inc_angle','start_date','end_date','ref_site'])) or 'GPS' in source:
            excerpts.append({'cell':i,'source':source[:6000],
                             'text_outputs':[''.join(o.get('text',[]) or o.get('data',{}).get('text/plain',[]))[:5000] for o in c.get('outputs',[])]})
    (OUT/f'{name}_metadata_excerpts.json').write_bytes((json.dumps(excerpts,indent=2)+'\n').encode())
