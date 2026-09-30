"""Compare the same Fresno census-block totals allocated by GHSL 100 m and CA-POP 100 m."""
from pathlib import Path
import os
from collections import defaultdict
import json, zipfile, hashlib
import numpy as np
import pandas as pd
import shapefile
import rasterio
from rasterio.windows import from_bounds, Window
from rasterio.warp import reproject, Resampling
from affine import Affine
from pyproj import Transformer
import shapely
from shapely.geometry import box, shape
from shapely.ops import transform as geom_transform
from shapely.strtree import STRtree
import matplotlib.pyplot as plt

ROOT=Path(os.environ.get('INSAR_STRENGTHENING_ROOT', Path(__file__).resolve().parent.parent))
INPUT=ROOT/'external'/'round3'
OUT=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))/'fresno_population_allocation'
OUT.mkdir(parents=True,exist_ok=True)
ROI=(-119.9,36.6,-119.7,36.85)
ZIP=INPUT/'ca2020.pl.zip'
BLOCKS=INPUT/'census_blocks'/'tl_2020_06019_tabblock20.shp'
CAPOP=INPUT/'CAPOP_2020_100m_TOTAL.tif'
DWR=ROOT/'results'/'dwr'/'grid.npz'

def census_counts():
    z=zipfile.ZipFile(ZIP)
    geography={}
    with z.open('cageo2020.pl') as src:
        for line in src:
            f=line.decode('latin1').rstrip('\n').split('|')
            if f[2]=='750' and f[8].startswith('7500000US06019'):
                geography[f[7]]=f[8].split('US',1)[1]
    counts={}
    with z.open('ca000012020.pl') as src:
        for line in src:
            f=line.decode('latin1').rstrip('\n').split('|')
            geoid=geography.get(f[4])
            if geoid is not None: counts[geoid]=int(f[5])
    if len(geography)!=len(counts): raise RuntimeError('Census geography/table join is incomplete.')
    return counts

def complete_blocks(counts):
    r=shapefile.Reader(str(BLOCKS))
    # Census TIGER/Line county shapefile is NAD83 longitude/latitude.
    to_wgs=Transformer.from_crs('EPSG:4269','EPSG:4326',always_xy=True).transform
    to_laea=Transformer.from_crs('EPSG:4269','EPSG:3310',always_xy=True).transform
    roi=box(*ROI);geoms=[];records=[]
    for rec,raw in zip(r.records(),r.shapes()):
        attrs=rec.as_dict();geoid=attrs['GEOID20']
        if geoid not in counts: raise RuntimeError('Census TIGER block missing a P1 record: '+geoid)
        g=shape(raw.__geo_interface__)
        if roi.covers(geom_transform(to_wgs,g)):
            geoms.append(geom_transform(to_laea,g));records.append((geoid,counts[geoid]))
    if len(records)!=6787 or sum(p for _,p in records)!=523953:
        raise RuntimeError(f'Unexpected complete-block baseline: {len(records)} / {sum(p for _,p in records)}')
    return geoms,records

def target_grid(blocks):
    with rasterio.open(CAPOP) as src:
        projected=shapely.union_all(blocks).bounds
        w=from_bounds(*projected,transform=src.transform)
        c0=max(0,int(np.floor(w.col_off))-2);r0=max(0,int(np.floor(w.row_off))-2)
        c1=min(src.width,int(np.ceil(w.col_off+w.width))+2);r1=min(src.height,int(np.ceil(w.row_off+w.height))+2)
        win=Window(c0,r0,c1-c0,r1-r0)
        cap=src.read(1,window=win,masked=True).astype('float64')
        cap=np.where(np.ma.getmaskarray(cap)|(~np.isfinite(cap.data))|(cap.data<0),0,cap.data)
        tr=src.window_transform(win);shape2=cap.shape
    with np.load(DWR) as z:
        x=z['x_edges'];y=z['y_edges'];pop=z['population'];area=z['area_m2'];obs=z['support']*area
    source_tr=Affine(x[1]-x[0],0,x[0],0,-(y[1]-y[0]),y[-1])
    warps=[]
    for arr in (pop,area,obs):
        out=np.zeros(shape2,dtype='float64')
        reproject(arr[::-1].copy(),out,src_transform=source_tr,src_crs='EPSG:4326',
                  dst_transform=tr,dst_crs='EPSG:3310',resampling=Resampling.sum,
                  src_nodata=None,dst_nodata=0,num_threads=2)
        warps.append(out)
    ghs,are,obsa=warps
    p0=float(pop.sum());a0=float(area.sum())
    mass_error=abs(float(ghs.sum())-p0)/p0
    area_error=abs(float(are.sum())-a0)/a0
    if mass_error>0.005 or area_error>0.005:
        raise RuntimeError(f'Conservative reprojection check failed: pop={mass_error:.4%}, area={area_error:.4%}')
    return cap,ghs,are,obsa,tr,{'input_GHSL_100m_domain_mass':p0,'reprojected_GHSL_domain_mass':float(ghs.sum()),'GHSL_mass_relative_error':mass_error,'input_area_m2':a0,'reprojected_area_m2':float(are.sum()),'area_mass_relative_error':area_error}

def intersect_and_score(blocks,records,grids,tr):
    cap,ghs,are,obs=grids
    h,w=cap.shape;dx=tr.a;dy=-tr.e;x0=tr.c;y0=tr.f
    bx=shapely.union_all(blocks).bounds
    c0=max(0,int(np.floor((bx[0]-x0)/dx))-1);c1=min(w,int(np.ceil((bx[2]-x0)/dx))+1)
    r0=max(0,int(np.floor((y0-bx[3])/dy))-1);r1=min(h,int(np.ceil((y0-bx[1])/dy))+1)
    rows,cols=np.mgrid[r0:r1,c0:c1]
    xx=x0+(cols.ravel()+.5)*dx; yy=y0-(rows.ravel()+.5)*dy
    cells=shapely.box(xx-dx/2,yy-dy/2,xx+dx/2,yy+dy/2)
    tree=STRtree(blocks)
    pairs=tree.query(cells,predicate='intersects')
    ci,bi=pairs
    inter=shapely.intersection(cells[ci],np.asarray(blocks,dtype=object)[bi])
    ia=shapely.area(inter); keep=ia>1e-5;ci=ci[keep];bi=bi[keep];ia=ia[keep]
    rr=rows.ravel()[ci];cc=cols.ravel()[ci]
    cell_area=dx*dy
    frac=ia/cell_area
    graw=ghs[rr,cc]*frac; craw=cap[rr,cc]*frac
    n=len(records);p=np.array([v for _,v in records],dtype='float64')
    area_sum=np.bincount(bi,weights=ia,minlength=n)
    gsum=np.bincount(bi,weights=graw,minlength=n)
    csum=np.bincount(bi,weights=craw,minlength=n)
    positive=p>0
    mg=(gsum>0)&positive;mc=(csum>0)&positive
    common=positive&mg&mc
    common_pop=float(p[common].sum())
    if common_pop < .95*float(p.sum()):
        raise RuntimeError(f'GHSL and CA-POP common allocation covers only {common_pop/p.sum():.2%} of the census population.')
    support=np.divide(obs,are,out=np.zeros_like(obs),where=are>0)
    support=np.clip(support,0,1)
    per={}
    for name,raw,den in [('Area_uniform',ia,area_sum),('GHSL_2020_100m',graw,gsum),('CA_POP_2020_100m',craw,csum)]:
        den_pair=den[bi]
        ok=common[bi] & (den_pair>0)
        alloc=np.zeros(ok.shape,dtype='float64')
        alloc[ok]=p[bi[ok]]*raw[ok]/den_pair[ok]
        total=float(alloc.sum());inside=float(np.sum(alloc*support[rr[ok],cc[ok]])) if False else 0.0
        # Compute in-place only for eligible pairs; no source pixel is counted twice.
        inside=float(np.sum(alloc[ok]*support[rr[ok],cc[ok]]))
        low=float(np.sum(alloc[ok]*(support[rr[ok],cc[ok]]<.5)))
        per[name]={'census_units_reallocated':total,'population_weighted_support_share':inside/total,
                   'support_deficit_units':total-inside,'below_50pct_support_units':low,
                   'blocks':int(common.sum()),'census_units_denominator':common_pop}
    if any(abs(v['census_units_reallocated']-common_pop)>1e-5 for v in per.values()):
        raise RuntimeError('Census block totals do not conserve within comparison domain.')
    # Per-block outputs identify local sign reversals and blocks that fail coverage.
    out=[]
    support_pair=support[rr,cc]
    cellsvalid=common[bi]
    pairfrac=np.zeros_like(ia);pairfrac[cellsvalid]=p[bi[cellsvalid]]/np.maximum(area_sum[bi[cellsvalid]],1e-30)*ia[cellsvalid]
    pairg=np.zeros_like(ia);pairc=np.zeros_like(ia)
    goodg=cellsvalid&(gsum[bi]>0);goodc=cellsvalid&(csum[bi]>0)
    pairg[goodg]=p[bi[goodg]]*graw[goodg]/gsum[bi[goodg]]
    pairc[goodc]=p[bi[goodc]]*craw[goodc]/csum[bi[goodc]]
    bsum=lambda v:np.bincount(bi,weights=v,minlength=n)
    block_def_area=bsum(pairfrac*(1-support_pair));block_def_g=bsum(pairg*(1-support_pair));block_def_c=bsum(pairc*(1-support_pair))
    # Scale, center-sampling and conditional covariance-bound experiment.
    scale_rows=[]
    px=x0+(cc+.5)*dx;py=y0-(rr+.5)*dy
    model_inputs=[('GHSL_2020_100m',graw,gsum),('CA_POP_2020_100m',craw,csum)]
    for model,raw,den in model_inputs:
        valid=common[bi]&(den[bi]>0)
        alloc=np.zeros_like(ia);alloc[valid]=p[bi[valid]]*raw[valid]/den[bi[valid]]
        for size in (250,500,1000,2000):
            for sx,sy in ((0,0),(size/2,0),(0,size/2),(size/2,size/2)):
                gx=np.floor((px+sx)/size).astype(np.int64);gy=np.floor((py+sy)/size).astype(np.int64)
                unique_grid,inv=np.unique(np.column_stack((gx,gy)),axis=0,return_inverse=True,return_index=False)
                ng=len(unique_grid)
                aa=np.bincount(inv,weights=ia,minlength=ng)
                pp=np.bincount(inv,weights=alloc,minlength=ng)
                ss=np.bincount(inv,weights=ia*support_pair,minlength=ng)
                hh=np.bincount(inv,weights=alloc*support_pair,minlength=ng)
                native=pp-hh;uniform=pp*(1-np.divide(ss,aa,out=np.zeros_like(ss),where=aa>0))
                delta=uniform-native
                mr=np.divide(pp,aa,out=np.zeros_like(pp),where=aa>0)
                ms=np.divide(ss,aa,out=np.zeros_like(ss),where=aa>0)
                mr2=np.bincount(inv,weights=np.divide(alloc*alloc,ia,out=np.zeros_like(alloc),where=ia>0),minlength=ng)/np.maximum(aa,1e-30)
                ms2=np.bincount(inv,weights=ia*support_pair*support_pair,minlength=ng)/np.maximum(aa,1e-30)
                vr=np.maximum(0,mr2-mr*mr);vs=np.maximum(0,ms2-ms*ms)
                bound=aa*np.sqrt(vr*vs)
                if np.any(np.abs(delta)>bound+1e-7):raise RuntimeError('Per-tile Cauchy bound does not cover exact within-tile error.')
                cx=(unique_grid[:,0]+.5)*size-sx;cy=(unique_grid[:,1]+.5)*size-sy
                pr=np.floor((y0-cy)/dy).astype(int);pc=np.floor((cx-x0)/dx).astype(int)
                ok=(pr>=0)&(pr<h)&(pc>=0)&(pc<w)
                center_s=np.zeros(ng);center_s[ok]=support[pr[ok],pc[ok]]
                center=pp*(1-center_s);dcenter=center-native
                n_top=max(1,int(np.ceil(.1*ng)))
                top_ref=set(np.argpartition(native,-n_top)[-n_top:].tolist())
                top_uni=set(np.argpartition(uniform,-n_top)[-n_top:].tolist())
                top_center=set(np.argpartition(center,-n_top)[-n_top:].tolist())
                base={'population_model':model,'scale_m':size,'origin_x_m':sx,'origin_y_m':sy,
                    'reporting_cells':ng,'domain_allocation_units':float(pp.sum()),
                    'native_deficit_units':float(native.sum()),'uniform_deficit_units':float(uniform.sum()),
                    'uniform_minus_native_pp':float(100*delta.sum()/common_pop),
                    'uniform_local_L1_pp':float(100*np.abs(delta).sum()/common_pop),
                    'uniform_local_cauchy_bound_pp':float(100*bound.sum()/common_pop),
                    'bound_to_observed_L1_ratio':float(bound.sum()/max(np.abs(delta).sum(),1e-30)),
                    'uniform_error_positive_cells':int(np.count_nonzero(delta>1e-10)),
                    'uniform_error_negative_cells':int(np.count_nonzero(delta<-1e-10)),
                    'uniform_top_decile_overlap':len(top_ref&top_uni)/n_top,
                    'center_minus_native_pp':float(100*dcenter.sum()/common_pop),
                    'center_local_L1_pp':float(100*np.abs(dcenter).sum()/common_pop),
                    'center_top_decile_overlap':len(top_ref&top_center)/n_top}
                scale_rows.append(base)
    scale_df=pd.DataFrame(scale_rows)
    for i,(gid,pop) in enumerate(records):
        out.append({'GEOID20':gid,'census_population_2020':pop,'included_common_population_domain':bool(common[i]),
                    'GHSL_native_cell_mass_in_block':gsum[i],'CA_POP_cell_mass_in_block':csum[i],
                    'support_deficit_area_uniform_units':block_def_area[i],
                    'support_deficit_GHSL_normalized_within_block_units':block_def_g[i],
                    'support_deficit_CA_POP_normalized_within_block_units':block_def_c[i]})
    byblock=pd.DataFrame(out)
    tract_df=byblock.copy();tract_df['tract_id']=tract_df['GEOID20'].str.slice(5,11)
    tract_df=(tract_df[tract_df['included_common_population_domain']]
       .groupby('tract_id',as_index=False).agg(census_population_units=('census_population_2020','sum'),
       GHSL_deficit=('support_deficit_GHSL_normalized_within_block_units','sum'),
       CA_POP_deficit=('support_deficit_CA_POP_normalized_within_block_units','sum'),
       area_uniform_deficit=('support_deficit_area_uniform_units','sum')))
    tract_df['CA_POP_minus_GHSL_deficit_units']=tract_df.CA_POP_deficit-tract_df.GHSL_deficit
    delta=tract_df['CA_POP_minus_GHSL_deficit_units'].to_numpy()
    tract_audit={'tract_count':int(len(tract_df)),'tracts_CA_POP_deficit_higher':int((delta>1e-8).sum()),
       'tracts_CA_POP_deficit_lower':int((delta<-1e-8).sum()),'net_CA_POP_minus_GHSL_units':float(delta.sum()),
       'summed_absolute_tract_difference_units':float(np.abs(delta).sum()),
       'cross_tract_sign_cancellation_fraction':float(1-abs(delta.sum())/max(np.abs(delta).sum(),1e-30)),
       'tract_resampling_inference':'Not performed: the full fixed study region is a finite descriptive domain, and no random sampling design is defined.'}
    return per,byblock,scale_df,tract_df,tract_audit,{'complete_blocks':n,'complete_block_census_population':float(p.sum()),
         'positive_population_blocks':int(positive.sum()),'blocks_with_GHSL_mass':int((mg|~positive).sum()),
         'blocks_with_CA_POP_mass':int((mc|~positive).sum()),'common_block_count':int(common.sum()),
         'common_census_population_units':common_pop,'common_population_share_of_all_complete_blocks':common_pop/float(p.sum()),
         'population_in_positive_blocks_without_GHSL_mass':float(p[positive&~mg].sum()),
         'population_in_positive_blocks_without_CA_POP_mass':float(p[positive&~mc].sum()),
         'candidate_cell_intersections':int(len(ia))}

def main():
    counts=census_counts();blocks,records=complete_blocks(counts)
    cap,ghs,are,obs,tr,warp_audit=target_grid(blocks)
    methods,byblock,scale_df,tract_df,tract_audit,domain=intersect_and_score(blocks,records,(cap,ghs,are,obs),tr)
    summary={'experiment':'Within-block normalization of alternative fine-resolution population distributions against the same 2020 Census block totals and DWR support field.',
       'scope':'Fresno fixed rectangle; complete 2020 Census blocks only; DWR 100-m measurement-cell polygon-union support.',
       'support':'Fractional observed area in 100-m projected target cells, reprojected from the area-weighted DWR mask field.',
       'models':methods,'domain':domain,'reprojection':warp_audit,'tract_decomposition':tract_audit,
       'interpretation':'CA-POP is an alternative dasymetric allocation model built from the same 2020 Census block counts and residential parcels/building footprints. Census block totals are normalized by construction. This is population-allocation sensitivity, not independent validation or household-level truth. Four partial block totals are not inferred; only blocks wholly inside the fixed ROI are retained.'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    byblock.to_csv(OUT/'per_block.csv',index=False)
    scale_df.to_csv(OUT/'scale_method_sensitivity.csv',index=False)
    tract_df.to_csv(OUT/'tract_decomposition.csv',index=False)
    pd.DataFrame([{'method':k,**v} for k,v in methods.items()]).to_csv(OUT/'method_summary.csv',index=False)
    labels=['Area uniform','GHSL 100 m','CA-POP 100 m']
    keys=['Area_uniform','GHSL_2020_100m','CA_POP_2020_100m']
    values=[1000*methods[k]['support_deficit_units']/methods[k]['census_units_denominator'] for k in keys]
    fig,ax=plt.subplots(figsize=(7.1,4.6),layout='constrained')
    bars=ax.bar(labels,values,color=['#8da0cb','#1b9e77','#d95f02'])
    ax.set_ylabel('Unsupported allocation units per 1,000 Census units')
    ax.set_ylim(0,max(values)*1.28);ax.grid(axis='y',alpha=.25);ax.set_axisbelow(True)
    for b,v in zip(bars,values):ax.text(b.get_x()+b.get_width()/2,v+max(values)*.025,f'{v:.2f}',ha='center',va='bottom',fontsize=9)
    ax.set_title('Fresno: same 2020 Census block totals, alternative 100-m allocations')
    ax.text(.01,-.20,f'Common denominator: {domain["common_census_population_units"]:,.0f} Census allocation units in {domain["common_block_count"]:,} complete blocks. Census block totals are enforced; CA-POP is a model sensitivity, not household truth.',transform=ax.transAxes,fontsize=8,va='top',wrap=True)
    fig.savefig(OUT/'support_share_comparison.png',dpi=220);fig.savefig(OUT/'support_share_comparison.pdf');plt.close(fig)
    print(json.dumps({'models':methods,'domain':domain,'reprojection':warp_audit,'tract_decomposition':tract_audit,'scale_rows':len(scale_df),'out':str(OUT)},indent=2))

if __name__=='__main__':main()
