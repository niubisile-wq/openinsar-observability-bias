from pathlib import Path
import json
import io
import zipfile
import hashlib
import pandas as pd
from source_inventory import fetch,RemoteZip
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'inputs/gnss_reference'
url='https://geodesy.unr.edu/NGLStationPages/DataHoldings.txt'
raw,_,_=fetch(url);(OUT/'DataHoldings_20261001.txt').write_bytes(raw)
records=[]
for line in raw.decode().splitlines()[1:]:
    f=line.split()
    if len(f)<11:continue
    try:
        lon=float(f[2]);lon=lon-360 if lon>180 else lon
        records.append({'station':f[0],'latitude':float(f[1]),'longitude':lon,'first_date':f[7],'last_date':f[8]})
    except ValueError:continue
hold=pd.DataFrame(records).drop_duplicates('station')
for track,bbox,ref in [('A064',[-118.4,33.85,-118.1,34.1],'P597'),('A137',[-121.9,38.4,-121.6,38.7],'P266')]:
    rates=pd.read_csv(OUT/f'{track}_archived_gnss_projected_rates.csv')
    joined=rates.merge(hold,on='station',how='left',validate='one_to_one')
    assert not joined.latitude.isna().any()
    joined.to_csv(OUT/f'{track}_archived_rates_with_coordinates.csv',index=False)
    roi=joined.longitude.between(bbox[0],bbox[2])&joined.latitude.between(bbox[1],bbox[3])
    print(track,'fixed ROI archived stations',joined[roi][['station','latitude','longitude']].to_dict('records'),'reference',joined[joined.station.eq(ref)].to_dict('records'),flush=True)
# Freeze the exact default GNSS trend model from the author's archived MintPy snapshot.
record=json.loads((ROOT/'inputs/provider_metadata/zenodo_19493073.json').read_text(encoding='utf-8'))
entry=next(f for f in record['files'] if f['key']=='software.zip')
with zipfile.ZipFile(RemoteZip(entry['links']['self'],entry['size'])) as z:
    candidates=[v for v in z.infolist() if v.filename.endswith('mintpy/objects/gnss.py') or v.filename.endswith('mintpy/objects/gps.py')]
    for item in candidates:
        content=z.read(item);(OUT/('provider_'+Path(item.filename).name)).write_bytes(content)
        print('Frozen GNSS implementation',item.filename,len(content),flush=True)
(OUT/'metadata_source_hashes.json').write_bytes((json.dumps({'NGL_holdings_sha256':hashlib.sha256(raw).hexdigest()},indent=2)+'\n').encode())
