"""Retain archived GNSS outputs as provenance, without executing provider notebooks."""
from pathlib import Path
import json
import re
import base64
import pandas as pd
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'inputs/gnss_reference'
OUT.mkdir(parents=True,exist_ok=True)
for track in ['A064','A137']:
    path=ROOT/'inputs/aria_california'/f'run_{track}_completed.ipynb'
    obj=json.loads(path.read_text(encoding='utf-8'))
    records=[];excerpts=[]
    for i,cell in enumerate(obj['cells']):
        source=''.join(cell.get('source',[]))
        if 'gnss_velocities = gnss.get_los_obs' in source:
            streams=''.join(''.join(o.get('text',[])) for o in cell.get('outputs',[]))
            matches=re.findall(r"\['([A-Z0-9]{4})'\s+'([-+0-9.eE]+)'\]",streams)
            assert len(matches)>100, len(matches)
            assert len(set(s for s,v in matches))==len(matches)
            for station,value in matches: records.append({'track':track,'station':station,'provider_gnss_los_velocity_mm_yr':float(value),'source_notebook_cell':i})
        if 'display(sites_df)' in source:
            for output in cell.get('outputs',[]):
                html=output.get('data',{}).get('text/html')
                if html: (OUT/f'{track}_provider_track_metadata.html').write_bytes(''.join(html).encode())
        if ('plot_insar_cartopy(' in source and 'ref_site =' in source) or 'plt.title(f"Velocities' in source:
            for j,output in enumerate(cell.get('outputs',[])):
                png=output.get('data',{}).get('image/png')
                if png:
                    content=''.join(png) if isinstance(png,list) else png
                    target=OUT/f'{track}_cell{i}_image{j}.png';target.write_bytes(base64.b64decode(content))
                    excerpts.append({'cell':i,'image':target.name})
    pd.DataFrame(records).to_csv(OUT/f'{track}_archived_gnss_projected_rates.csv',index=False)
    print(track,'archived projected station rates',len(records),'date-bearing figures',excerpts,flush=True)
