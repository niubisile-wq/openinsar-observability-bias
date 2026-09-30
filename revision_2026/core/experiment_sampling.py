"""E4 temporal robustness and expansion, with conditional sample variation."""
from datetime import datetime
import json
import numpy as np
import pandas as pd
import networkx as nx
import rasterio
from common import OUT,ROOT,CFG,dump,exact_joint_support,resolve_recorded_path

DEST=OUT/'sampling';DEST.mkdir(exist_ok=True)
PAIR=DEST/'pairs';PAIR.mkdir(exist_ok=True)

def network(pairs):
    graph=nx.Graph();graph.add_edges_from(p.split('_') for p in pairs)
    comps=list(nx.connected_components(graph));dates=sorted(graph.nodes)
    return dict(pairs=len(pairs),dates=len(dates),components=len(comps),largest_component=max(map(len,comps)),cycle_rank=len(pairs)-len(dates)+len(comps),start=dates[0],end=dates[-1])

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--frozen-summaries',action='store_true',help='Replay E4 from per-pair integrals and frozen nested support fields; raw pixel transfer is separate.')
    frozen=parser.parse_args().frozen_summaries
    z=np.load(OUT/'fine/grid.npz');tx,ty=z['x_edges'],z['y_edges']
    weights={e:z[e]*z['domain'] for e in ['population','builtup_m2','area_m2']}
    records=json.loads((ROOT/'external/licsar_expanded/sampling_manifest.json').read_text(encoding='utf-8'))
    if len(records)!=224 or any(r['status']!='ok' for r in records):raise ValueError('Sampling transfer incomplete')
    design=pd.read_csv(ROOT/'external/licsar_expanded/experiment_input_plan.csv').set_index('pair')
    sums={n:{k:np.zeros(weights['population'].shape,dtype='float64') for k in ['unw','q20','q30','q40']} for n in [1,2,4,8]}
    rows=[]
    for i,record in enumerate([] if frozen else sorted(records,key=lambda r:r['pair'])):
        pair=record['pair'];path=PAIR/(pair+'.npz')
        legacy=OUT/'fine/pairs'/(pair+'.npz')
        if legacy.exists():path=legacy
        if path.exists():
            with np.load(path) as data:support={k:data[k] for k in ['unw','q20','q30','q40']}
        else:
            with rasterio.open(resolve_recorded_path(record['files']['cc']['path'])) as cd,rasterio.open(resolve_recorded_path(record['files']['unw']['path'])) as ud:
                cc=cd.read(1);unw=ud.read(1)
                support=exact_joint_support(cc,(cd.read_masks(1)>0)&np.isfinite(cc),cd.transform,(ud.read_masks(1)>0)&np.isfinite(unw)&(unw!=0),ud.transform,tx,ty,CFG['coherence_thresholds'])
            np.savez_compressed(path,**support)
        row=design.loc[pair]
        for k,s in support.items():
            for n in sums:
                if row.sampling_rank<=n:sums[n][k]+=s
            for e,w in weights.items():rows.append(dict(pair=pair,support=k,exposure=e,year=int(row.year),quarter=int(row.quarter),days=int(row.days),baseline_group=row.baseline_group,sampling_rank=int(row.sampling_rank),total=float(w.sum()),unsupported=float((w*(1-s.astype('float64'))).sum())))
        if (i+1)%10==0:print('EXPANSION',i+1,'/224',flush=True)
    if frozen:
        df=pd.read_csv(DEST/'expanded_pair_exposure.csv')
        if set(df.pair)!=set(r['pair'] for r in records):raise ValueError('Frozen pair-integral catalogue mismatch')
    else:
        df=pd.DataFrame(rows);df.to_csv(DEST/'expanded_pair_exposure.csv',index=False)
    summary=[];curves=[]
    for n,fields in sums.items():
        if frozen:
            # Explicit float64 comparison matches the original NumPy 2 scalar
            # promotion at support-distribution thresholds under NumPy 1.26 too.
            with np.load(DEST/f'mean_{28*n}.npz') as z:fields={k:z[k].astype('float64') for k in z.files}
        else:
            fields={k:(s/(28*n)).astype('float32') for k,s in fields.items()}
            np.savez_compressed(DEST/f'mean_{28*n}.npz',**fields)
        for k,s in fields.items():
            for e,w in weights.items():
                summary.append(dict(pairs=28*n,support=k,exposure=e,total=float(w.sum()),unsupported=float((w*(1-s.astype('float64'))).sum())))
                for t in np.linspace(0,1,101):curves.append(dict(pairs=28*n,support=k,exposure=e,threshold=t,below_share=float(w[s<t].sum()/w.sum())))
    pd.DataFrame(summary).to_csv(DEST/'nested_summary.csv',index=False)
    pd.DataFrame(curves).to_csv(DEST/'nested_curves.csv',index=False)
    # Leave-one-year, season and baseline comparisons for both legacy and pool.
    legacy=pd.read_csv(OUT/'fine/pair_exposure.csv');legacy=legacy[legacy.domain.eq('geographic')].copy()
    legacy['year']=legacy.pair.str[:4].astype(int);legacy['quarter']=(legacy.pair.str[4:6].astype(int)-1)//3+1
    legacy['days']=legacy.pair.map(lambda p:(datetime.strptime(p[9:],'%Y%m%d')-datetime.strptime(p[:8],'%Y%m%d')).days)
    robustness=[]
    for name,data in [('legacy28',legacy),('expanded224',df)]:
        for support,part in data.groupby('support'):
            for exposure,part2 in part.groupby('exposure'):
                filters=[('all',np.ones(len(part2),dtype=bool))]
                filters += [('leave_'+str(y),part2.year.ne(y).values) for y in range(2016,2023)]
                filters += [('year_'+str(y),part2.year.eq(y).values) for y in range(2016,2023)]
                filters += [('quarter_'+str(q),part2.quarter.eq(q).values) for q in range(1,5)]
                filters += [('baseline_le24',part2.days.le(24).values),('baseline_gt24',part2.days.gt(24).values),('baseline_gt96',part2.days.gt(96).values)]
                for label,keep in filters:
                    sub=part2.loc[keep]
                    if len(sub):robustness.append(dict(sample=name,support=support,exposure=exposure,subset=label,pairs=len(sub),mean_unsupported=float(sub.unsupported.mean()),share_unsupported=float((sub.unsupported/sub.total).mean())))
    pd.DataFrame(robustness).to_csv(DEST/'temporal_robustness.csv',index=False)
    rng=np.random.default_rng(CFG['seed']);replicates=[]
    pool=df[(df.support=='q30')&(df.exposure=='population')].set_index('pair')
    strata=[part.index.values for _,part in pool.groupby(['year','quarter'])]
    for rep in range(100):
        permutations=[rng.permutation(s) for s in strata]
        for n in [1,2,4,8]:
            ids=np.concatenate([s[:n] for s in permutations]);sub=pool.loc[ids]
            replicates.append(dict(replicate=rep,pairs=len(ids),unsupported=float(sub.unsupported.mean()),unsupported_share=float((sub.unsupported/sub.total).mean())))
    pd.DataFrame(replicates).to_csv(DEST/'conditional_sampling_replicates.csv',index=False)
    nets={name:network(list(data.pair.unique())) for name,data in [('legacy28',legacy),('expanded224',df)]}
    nets['catalogue2016_2022']=network(list(design.index))
    dump(DEST/'network.json',nets)
    dump(DEST/'method.json',dict(primary='Equal area-pair mean, identical fine allocation and geographic footprint',nested='28/56/112/224 pairs; four calendar-quarter strata per year, ranks prescribed before outcome extraction',resampling='100 equal-quarter samples drawn without replacement from the fixed 224-pair pool; intervals quantify conditional sampling variation, not population parameter confidence',temporal='Year and quarter refer to first acquisition date; long-baseline interferograms can span other seasons',representativeness='Design deliberately includes long baselines; equal-pair means do not estimate final-product validity or a catalogue-frequency-weighted average',networks=nets))
    print('E4 COMPLETE',flush=True)

if __name__=='__main__':main()
