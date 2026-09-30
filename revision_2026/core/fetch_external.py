"""Fetch source metadata and one external population layer, with integrity records."""
from pathlib import Path
from urllib.parse import urljoin
import re,json,hashlib,time,requests

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'external';OUT.mkdir(exist_ok=True)
session=requests.Session()
session.headers['User-Agent']='InSAR scientific revision reproducibility research'
session.headers['Accept-Encoding']='identity'

def fetch(url,path):
    path=Path(path)
    if not path.exists():
        for attempt in range(4):
            try:
                with session.get(url,stream=True,timeout=(25,90)) as r:
                    r.raise_for_status(); total=int(r.headers.get('Content-Length',0));n=0;last=time.time()
                    with path.with_suffix(path.suffix+'.part').open('wb') as f:
                        for block in r.iter_content(1024*1024):
                            f.write(block);n+=len(block)
                            if time.time()-last>20:
                                print(path.name,n,total,flush=True);last=time.time()
                    if total and not r.headers.get('Content-Encoding') and total!=n:raise IOError(f'Length mismatch {n}/{total}')
                path.with_suffix(path.suffix+'.part').replace(path)
                break
            except Exception as e:
                print('RETRY',attempt+1,str(e),flush=True)
                if attempt==3:raise
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return dict(url=url,file=path.name,bytes=path.stat().st_size,sha256=h.hexdigest())

def main():
    records=[]
    page='https://hub.worldpop.org/geodata/summary?id=29165'
    records.append(fetch(page,OUT/'worldpop_29165.html'))
    html=(OUT/'worldpop_29165.html').read_text(encoding='utf-8')
    links=re.findall(r'href=[\"\']([^\"\']+)[\"\']',html)
    urls=[urljoin(page,x).replace('&amp;','&') for x in links if 'data.worldpop.org' in x]
    print('WORLDPOP_DOWNLOAD_LINKS',urls,flush=True)
    tif=next(x for x in urls if '.tif' in x)
    records.append(fetch(tif,OUT/'worldpop_thailand_2020_UNadj.tif'))
    for repo in ['yumorishita/LiCSBAS2','comet-licsar/LiCSBAS','copernicus-land/egms-api']:
        safe=repo.replace('/','_')
        try:
            records.append(fetch(f'https://api.github.com/repos/{repo}',OUT/f'{safe}_repository.json'))
            branch=json.loads((OUT/f'{safe}_repository.json').read_text())['default_branch']
            records.append(fetch(f'https://api.github.com/repos/{repo}/git/trees/{branch}?recursive=1',OUT/f'{safe}_tree.json'))
        except Exception as e:print('METADATA_ERROR',repo,str(e),flush=True)
    (OUT/'download_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print('FINISHED',len(records),flush=True)

if __name__=='__main__':main()
