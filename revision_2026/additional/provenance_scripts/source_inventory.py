"""Read provider metadata and ZIP catalogues before choosing any outcomes."""
from pathlib import Path
import io
import json
import struct
import zipfile
import urllib.request
import urllib.parse
import hashlib
import time

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'inputs' / 'provider_metadata'
OUT.mkdir(parents=True, exist_ok=True)

def fetch(url, headers=None):
    req = urllib.request.Request(url, headers={'User-Agent': 'InSAR-revision-reproducibility/1.0', **(headers or {})})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read(), {k.lower():v for k,v in r.headers.items()}, r.status

def save(name, obj):
    (OUT / name).write_bytes((json.dumps(obj, ensure_ascii=False, indent=2)+'\n').encode('utf-8'))

class RemoteZip(io.RawIOBase):
    def __init__(self, url, size):
        self.url = url; self.size = size; self.pos = 0
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, offset, whence=0):
        self.pos = offset if whence==0 else (self.pos+offset if whence==1 else self.size+offset)
        return self.pos
    def read(self, n=-1):
        n = self.size-self.pos if n<0 else min(n,self.size-self.pos)
        if n<=0: return b''
        start = self.pos; end = start+n-1
        # Distinct query avoids an intermediary returning a cached, different range.
        sep = '&' if '?' in self.url else '?'
        url = self.url+sep+f'revision_range={start}-{end}'
        data, headers, status = fetch(url, {'Range': f'bytes={start}-{end}'})
        if status!=206 or len(data)!=n or not headers.get('content-range','').startswith(f'bytes {start}-{end}/'):
            raise RuntimeError(f'Provider did not honor requested range: {status} {headers.get("content-range")} {len(data)}/{n}')
        self.pos+=n
        return data

def main():
    for record in ['19493073', '19120798']:
        try:
            raw, headers, status = fetch(f'https://zenodo.org/api/records/{record}')
            (OUT / f'zenodo_{record}.json').write_bytes(raw)
            obj=json.loads(raw)
            print(record, 'license', obj['metadata'].get('license'), flush=True)
            for file in obj['files']:
                print(file['key'], file['size'], flush=True)
                if file['key'] in ['velocities.zip','coherence.zip','validation_notebooks.zip','cmdline_cli.zip',
                                   'NC_ALOS2_170_10_56_P345_results.zip','NC_ALOS2_068_10_56_P345_results.zip','figure_inputs.zip']:
                    url=file['links']['self']
                    with zipfile.ZipFile(RemoteZip(url,file['size'])) as z:
                        entries=[dict(name=f.filename,size=f.file_size,compressed=f.compress_size,crc=f.CRC,
                                      header_offset=f.header_offset,compression=f.compress_type) for f in z.infolist()]
                    save(f'{record}_{file["key"]}_catalogue.json',entries)
                    print('catalogue', len(entries), [v['name'] for v in entries[:8]], flush=True)
        except Exception as e:
            save(f'zenodo_{record}_error.json', {'error':repr(e)})
            print('ERROR', record, repr(e), flush=True)
    url='https://gis4u.fresno.gov/exa/rest/services/PublicInfoServices/AddrParcelStreet/MapServer/1'
    try:
        raw, _, _=fetch(url+'?f=pjson');(OUT/'fresno_city_parcel_layer.json').write_bytes(raw)
        fields=json.loads(raw).get('fields',[])
        print('Parcel fields', [f['name'] for f in fields], flush=True)
        query=urllib.parse.urlencode({'where':'1=1','outFields':'EXISTING_LAND_USE,EXISTING_LAND_USE_TEXT',
                                     'returnDistinctValues':'true','returnGeometry':'false','f':'json'})
        raw,_,_=fetch(url+'/query?'+query);(OUT/'fresno_city_landuse_categories.json').write_bytes(raw)
        print('Parcel categories bytes',len(raw),flush=True)
    except Exception as e:
        save('fresno_parcel_error.json',{'error':repr(e)});print('ERROR parcel',repr(e),flush=True)

if __name__=='__main__': main()
