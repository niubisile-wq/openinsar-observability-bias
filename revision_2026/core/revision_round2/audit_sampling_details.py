"""Independent temporal group and spatial-ranking arithmetic checks."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
from scipy.stats import rankdata
H=Path(__file__).resolve().parent;R=H.parent
def main():
    old=pd.read_csv(R/'results/fine/pair_exposure.csv');old=old[old.domain=='geographic'].copy()
    new=pd.read_csv(R/'results/sampling/expanded_pair_exposure.csv')
    for d in [old,new]:
        first=pd.to_datetime(d.pair.str[:8],format='%Y%m%d');last=pd.to_datetime(d.pair.str[9:],format='%Y%m%d')
        for k,v in dict(year=first.dt.year,quarter=first.dt.quarter,days=(last-first).dt.days).items():
            if k in d:assert np.array_equal(d[k],v)
            d[k]=v
    plan=pd.read_csv(R/'external/licsar_expanded/experiment_input_plan.csv').set_index('pair')
    manifest=json.loads((R/'external/licsar_expanded/sampling_manifest.json').read_text())
    assert len(manifest)==224 and all(r['status']=='ok' for r in manifest)
    pairs=set(r['pair'] for r in manifest);assert pairs==set(new.pair)
    selected=plan.loc[sorted(pairs)]
    for _,d in selected.groupby(['year','quarter']):assert sorted(d.sampling_rank)==list(range(1,9))
    temporal=pd.read_csv(R/'results/sampling/temporal_robustness.csv');errors=[]
    for r in temporal.itertuples():
        d=(old if r.sample=='legacy28' else new);d=d[(d.support==r.support)&(d.exposure==r.exposure)]
        s=r.subset
        if s.startswith('leave_'):d=d[d.year!=int(s[6:])]
        elif s.startswith('year_'):d=d[d.year==int(s[5:])]
        elif s.startswith('quarter_'):d=d[d.quarter==int(s[8:])]
        elif s=='baseline_le24':d=d[d.days<=24]
        elif s=='baseline_gt24':d=d[d.days>24]
        elif s=='baseline_gt96':d=d[d.days>96]
        else:assert s=='all'
        assert len(d)==r.pairs
        err=abs(d.unsupported.sum()/len(d)-r.mean_unsupported);errors.append(err)
        assert err<1e-5
        assert abs((d.unsupported/d.total).sum()/len(d)-r.share_unsupported)<1e-12
    blocks=pd.read_csv(R/'results/sampling/spatial_blocks.csv');stability=pd.read_csv(R/'results/sampling/spatial_rank_stability.csv')
    ref=blocks[(blocks['sample']=='nested224')&blocks.eligible].sort_values('block_id')
    for r in stability.itertuples():
        d=blocks[(blocks['sample']==r.sample)&blocks.eligible].sort_values('block_id');assert list(d.block_id)==list(ref.block_id)
        a=d[r.metric].to_numpy();b=ref[r.metric].to_numpy()
        rho=np.corrcoef(rankdata(a),rankdata(b))[0,1]
        topa=set(d.block_id.iloc[np.argsort(a,kind='stable')[-r.top_count:]])
        topb=set(ref.block_id.iloc[np.argsort(b,kind='stable')[-r.top_count:]])
        assert abs(rho-r.spearman)<1e-12 and abs(len(topa&topb)/len(topa|topb)-r.top_decile_jaccard)<1e-12
    v=dict(temporal_rows=len(temporal),max_absolute_arithmetic_difference=max(errors),ranking_rows=len(stability),
       manifest_all_224_ok=True,year_quarter_rank_groups=28,scope='Temporal and ranking arithmetic independently recomputed from retained pair/block tables; spatial blocks not independently regenerated.')
    (H/'evidence/sampling_details_audit.json').write_text(json.dumps(v,indent=2));print(v)
if __name__=='__main__':main()
