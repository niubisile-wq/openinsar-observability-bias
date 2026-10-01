from pathlib import Path
import json
import zipfile
import hashlib
import sys
from source_inventory import RemoteZip, fetch

ROOT=Path(__file__).resolve().parent.parent
META=ROOT/'inputs/provider_metadata'
OUT=ROOT/'inputs/aria_california'
OUT.mkdir(parents=True,exist_ok=True)

def main():
    record=json.loads((META/'zenodo_19493073.json').read_text(encoding='utf-8'))
    selected={
        'cmdline_cli.zip': lambda n: not n.endswith('/'),
        'validation_notebooks.zip': lambda n: 'run_A064_' in n or 'run_A137_' in n,
        'velocities.zip': lambda n: ('A064_velocity_' in n or 'A137_velocity_' in n),
        'coherence.zip': lambda n: 'A064_avg' in n or 'A137_avg' in n,
    }
    records=[]
    for key,predicate in selected.items():
        entry=next(f for f in record['files'] if f['key']==key)
        with zipfile.ZipFile(RemoteZip(entry['links']['self'],entry['size'])) as z:
            for item in z.infolist():
                if not predicate(item.filename): continue
                # File names from provider, flattened intentionally to selected-product directory.
                target=OUT/Path(item.filename).name
                if target.exists() and target.stat().st_size==item.file_size:
                    raw=target.read_bytes()
                else:
                    print('Fetching',item.filename,item.compress_size,flush=True)
                    raw=z.read(item)
                    target.write_bytes(raw)
                records.append({'source_record':'10.5281/zenodo.19493073','archive':key,
                                'member':item.filename,'source_crc32':item.CRC,'bytes':len(raw),
                                'sha256':hashlib.sha256(raw).hexdigest(),'output':target.name})
                print('Saved',target.name,len(raw),flush=True)
        (OUT/'selected_member_provenance.json').write_bytes((json.dumps(records,indent=2)+'\n').encode())

if __name__=='__main__': main()
