"""Predeclare stratified sampling and matched short-baseline processing inputs."""
from pathlib import Path
from datetime import datetime
from urllib.parse import urljoin
import argparse,concurrent.futures,json,re,time,threading
import numpy as np,pandas as pd,requests,rasterio
from common import ROOT,CFG,PROJECT,BACKUP,dump,sha256

DEST=ROOT/'external/licsar_expanded';DEST.mkdir(parents=True,exist_ok=True)
BASE=CFG['catalogue_url'];LOCK=threading.Lock()

def plan():
    p=DEST/'experiment_input_plan.csv'
    if p.exists():return pd.read_csv(p)
    r=requests.get(BASE,timeout=45);r.raise_for_status()
    (DEST/'catalogue.html').write_text(r.text,encoding='utf-8')
    pairs=sorted(set(re.findall(r'\b(\d{8}_\d{8})\b',r.text)))
    rows=[]
    for pair in pairs:
        a,b=pair.split('_')
        if a<'20160101' or b>'20221231':continue
        days=(datetime.strptime(b,'%Y%m%d')-datetime.strptime(a,'%Y%m%d')).days
        rows.append(dict(pair=pair,year=int(a[:4]),quarter=(int(a[4:6])-1)//3+1,days=days,baseline_group='short' if days<=24 else 'long'))
    df=pd.DataFrame(rows);df['sampling_rank']=0;df['timeseries']=False
    rng=np.random.default_rng(CFG['seed'])
    for _,part in df.groupby(['year','quarter']):
        short=rng.permutation(part[part.baseline_group.eq('short')].index).tolist()
        long=rng.permutation(part[part.baseline_group.eq('long')].index).tolist()
        order=[]
        while len(order)<8 and (short or long):
            for seq in [short,long]:
                if seq and len(order)<8:order.append(seq.pop())
        if len(order)<8:raise ValueError('Fewer than eight pairs in a year-quarter stratum')
        for rank,i in enumerate(order,1):
            # Balance baseline classes even in the 28-pair nested subset.
            balanced_rank=(rank+1 if rank%2 else rank-1) if int(part.quarter.iloc[0])%2==0 else rank
            df.loc[i,'sampling_rank']=balanced_rank
    # Two-year matched processing period, selected by time/baseline before examining displacement.
    df['timeseries']=(df.pair.str[:8]>='20190101')&(df.pair.str[9:]<='20201231')&(df.days<=96)
    df.to_csv(p,index=False)
    dump(DEST/'design.json',dict(seed=CFG['seed'],sampling='8 pairs per year-quarter, alternating <=24 and >24 day baselines; nested ranks 1,2,4,8 give 28,56,112,224 pairs; no selection by support or deformation result',timeseries='All listed pairs with both dates in 2019-2020 and temporal baseline <=96 days',sampling_count=int((df.sampling_rank>0).sum()),timeseries_count=int(df.timeseries.sum()),catalogue_sha256=sha256(DEST/'catalogue.html')))
    return df

def fetch_pair(pair):
    result={'pair':pair,'files':{}}
    index=None
    try:
        for ext in ['cc','unw']:
            filename=f'{pair}.geo.{ext}.tif'
            original=(PROJECT/'raw/licsar'/filename) if ext=='cc' else (BACKUP/'revision_scientific_reports_v1/inputs_needed/licsar_chao_28_unw'/filename)
            path=original if original.exists() else DEST/filename
            url=None
            if not path.exists():
                if index is None:
                    r=requests.get(BASE+pair,timeout=45);r.raise_for_status();index=r.text
                hrefs=re.findall(r'href=[\"\']([^\"\']+)[\"\']',index)
                url=urljoin(BASE+pair+'/',next(x for x in hrefs if x.endswith(f'.geo.{ext}.tif')))
                for attempt in range(4):
                    try:
                        with requests.get(url,stream=True,timeout=(30,100)) as r:
                            r.raise_for_status();n=0;expected=int(r.headers.get('Content-Length',0))
                            with path.with_suffix('.tif.part').open('wb') as f:
                                for chunk in r.iter_content(2**20):f.write(chunk);n+=len(chunk)
                        if expected and expected!=n:raise IOError('Download length mismatch')
                        with rasterio.open(path.with_suffix('.tif.part')) as d:d.read(1,window=((0,1),(0,1)))
                        path.with_suffix('.tif.part').replace(path)
                        break
                    except Exception:
                        if attempt==3:raise
            with rasterio.open(path) as d:
                result['files'][ext]=dict(path=str(path),url=url,sha256=sha256(path),bytes=path.stat().st_size,shape=list(d.shape),transform=list(d.transform)[:6],crs=str(d.crs))
        result['status']='ok'
    except Exception as e:result.update(status='error',error=str(e))
    with LOCK:
        with (DEST/'transfer_log.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(result)+'\n')
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--group',choices=['sampling','timeseries'],default='sampling');ap.add_argument('--plan-only',action='store_true');args=ap.parse_args()
    df=plan();selected=df[df.sampling_rank>0] if args.group=='sampling' else df[df.timeseries]
    print('PLAN',args.group,len(selected),flush=True)
    if args.plan_only:return
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        fs={pool.submit(fetch_pair,p):p for p in selected.pair}
        for n,f in enumerate(concurrent.futures.as_completed(fs),1):
            row=f.result();results.append(row);print(f'{n}/{len(fs)} {row["pair"]} {row["status"]}',row.get('error',''),flush=True)
            dump(DEST/(args.group+'_manifest.json'),sorted(results,key=lambda x:x['pair']))
    errors=[x for x in results if x['status']!='ok']
    if errors:raise RuntimeError(f'{len(errors)} pairs failed; inspect manifest and rerun')
    print('COMPLETE',args.group,len(results),flush=True)

if __name__=='__main__':main()
