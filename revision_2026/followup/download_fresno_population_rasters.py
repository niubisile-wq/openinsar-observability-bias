import os
from pathlib import Path
import hashlib,requests

OUT=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))/'inputs'/'fresno_population'
OUT.mkdir(parents=True,exist_ok=True)
BASE='https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_54009_1000/V1-0/tiles/'
items={
 'ghsl_2020_1km_fresno_tile.zip':BASE+'GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0_R5_C8.zip',
 'worldpop_2020_1km_us_unadj.tif':'https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/2020/USA/usa_ppp_2020_1km_Aggregated_UNadj.tif',
}
for name,url in items.items():
 p=OUT/name
 if not p.exists():
  with requests.get(url,stream=True,timeout=(20,90)) as r:
   r.raise_for_status()
   length=int(r.headers.get('Content-Length','0'))
   limit=100_000_000
   if length and length>limit: raise RuntimeError(f'Refusing unexpected large download {name}: {length}')
   with p.open('wb') as f:
    for chunk in r.iter_content(1024*1024):
     if chunk:f.write(chunk)
 h=hashlib.sha256(p.read_bytes()).hexdigest()
 print(name,p.stat().st_size,h)
