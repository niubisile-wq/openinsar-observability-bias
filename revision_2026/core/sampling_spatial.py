"""Descriptive spatial rank stability on fixed geographic blocks, not utility validation."""
import json, math
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from common import OUT,dump

def main():
    g=np.load(OUT/'fine/grid.npz');d=g['domain'];p=g['population']*d;a=g['area_m2']*d
    x=(g['x_edges'][1:]+g['x_edges'][:-1])/2;y=(g['y_edges'][1:]+g['y_edges'][:-1])/2
    # Geographic 0.05-degree blocks anchored at integer multiples of 0.05 degrees.
    # Assign native pixels by their centers; this is aggregation, not reallocation.
    bx=np.floor(x/.05+1e-9).astype(int);by=np.floor(y/.05+1e-9).astype(int)
    nx=int(bx.max()-bx.min()+1);ny=int(by.max()-by.min()+1)
    ids=((by[:,None]-by.min())*nx+(bx[None,:]-bx.min())).ravel()
    total=np.bincount(ids,weights=p.ravel(),minlength=nx*ny)
    area=np.bincount(ids,weights=a.ravel(),minlength=nx*ny)
    # Fixed numeric inclusion rule recorded before looking at rankings.
    eligible=(total>=1000)&(area>=5e6)
    rows=[];scores={}
    fields={'legacy28':OUT/'fine/mean_support.npz',**{f'nested{n}':OUT/f'sampling/mean_{n}.npz' for n in [28,56,112,224]}}
    for label,path in fields.items():
        s=np.load(path)['q30'];u=np.bincount(ids,weights=(p*(1-s.astype('float64'))).ravel(),minlength=nx*ny)
        share=np.divide(u,total,out=np.zeros_like(u),where=total>0);scores[label]={'deficit':u,'deficit_share':share}
        for i in np.flatnonzero(area>0):
            rows.append(dict(sample=label,block_id=int(i),west=float((bx.min()+i%nx)*.05),south=float((by.min()+i//nx)*.05),area_m2=float(area[i]),population=float(total[i]),deficit=float(u[i]),deficit_share=float(share[i]),eligible=bool(eligible[i])))
    pd.DataFrame(rows).to_csv(OUT/'sampling/spatial_blocks.csv',index=False)
    comparison=[];n=int(eligible.sum());k=max(1,math.ceil(.1*n));indices=np.flatnonzero(eligible)
    for label,v in scores.items():
        for metric,values in v.items():
            ref=scores['nested224'][metric];aa=set(indices[np.argsort(values[eligible],kind='stable')[-k:]]);bb=set(indices[np.argsort(ref[eligible],kind='stable')[-k:]])
            comparison.append(dict(sample=label,reference='nested224',metric=metric,blocks=n,top_count=k,spearman=float(spearmanr(values[eligible],ref[eligible]).statistic),top_decile_jaccard=len(aa&bb)/len(aa|bb)))
    pd.DataFrame(comparison).to_csv(OUT/'sampling/spatial_rank_stability.csv',index=False)
    dump(OUT/'sampling/spatial_method.json',dict(blocks='0.05-degree longitude/latitude, integer-multiple origin, native pixels assigned by center; all blocks retained in CSV',eligibility='At least 1,000 GHSL population units and 5 square kilometers in primary domain; fixed before extracting rankings',eligible_blocks=n,interpretation='Descriptive ranking changes across selected samples; nested224 is a comparator, not truth or evidence of better measurement priorities. Spatial blocks are not treated as independent replicates.'))
    print(pd.DataFrame(comparison).to_string(index=False),flush=True)

if __name__=='__main__':main()
