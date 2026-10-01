from pathlib import Path
import hashlib
import json
import zipfile
from source_inventory import fetch
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'inputs/population'
OUT.mkdir(parents=True,exist_ok=True)
name='GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0_R6_C6'
url='https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_4326_3ss/V1-0/tiles/'+name+'.zip'
target=OUT/(name+'.zip')
if not target.exists():
    print('Fetching GHSL west coast tile',flush=True)
    raw,headers,status=fetch(url);target.write_bytes(raw)
with zipfile.ZipFile(target) as z:
    files=[]
    for item in z.infolist():
        if item.filename.lower().endswith('.tif'):
            dest=OUT/Path(item.filename).name
            if not dest.exists(): dest.write_bytes(z.read(item))
            files.append(dest.name)
metadata={'source':url,'population_epoch':2020,'release':'R2023A','units':'Count allocated per grid cell',
          'zip_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'bytes':target.stat().st_size,'extracted':files}
(OUT/'GHSL_west_provenance.json').write_bytes((json.dumps(metadata,indent=2)+'\n').encode())
print(metadata,flush=True)
