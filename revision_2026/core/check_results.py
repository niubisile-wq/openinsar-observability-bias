"""Numerical evidence checks; no statistical accuracy claims from these tests."""
import json
import numpy as np,pandas as pd
from affine import Affine
from rasterio.warp import reproject,transform_bounds,Resampling
from scipy.stats import spearmanr
from common import OUT,ROOT,dump,overlaps

def main():
    checks={}
    base=pd.read_csv(OUT/'fine/allocation_summary.csv')
    checks['max_baseline_conservation_relative']=float(np.max(np.abs(base.supported+base.unsupported-base.total)/base.total))
    scale=pd.read_csv(OUT/'allocation/scale_origin_comparison.csv')
    uni=scale[scale.method.eq('area_uniform')]
    checks['max_covariance_identity_absolute']=float(np.max(np.abs(uni.signed_error-uni.covariance_integral)))
    checks['allocation_rows']=len(scale)
    assert checks['max_baseline_conservation_relative']<1e-10
    assert checks['max_covariance_identity_absolute']<1e-5
    # Analytic counterexamples: known fine populations [9,1], equal-area cells.
    e=np.array([9.,1.]);support=np.array([1.,0.]);fine=np.dot(e,1-support);uniform=e.sum()*(1-support.mean())
    assert fine==1 and uniform==5
    assert np.dot(e,support)-e.sum()*support.mean()==uniform-fine
    # Complementing support reverses the sign; uniform density removes bias.
    assert np.dot(e,support)-e.sum()*support.mean()==-(np.dot(e,1-support)-e.sum()*(1-support).mean())
    assert np.dot(np.ones(2),support)-2*support.mean()==0
    checks['analytic_covariance_positive_negative_zero']='pass'
    # Full empirical distribution: exact integral of F(t) is computed by its
    # steps, not by declaring a coarse plotted curve to be exact.
    z=np.load(OUT/'fine/grid.npz');s=np.load(OUT/'fine/mean_support.npz')['q30'].astype('float64');p=z['population']*z['domain'];a=z['area_m2']*z['domain']
    order=np.argsort(s.ravel());sorted_s=s.ravel()[order];mass=p.ravel()[order];cum=np.cumsum(mass)
    integ=float(np.dot(cum[:-1],np.diff(sorted_s))+cum[-1]*(1-sorted_s[-1]))
    direct=float((p*(1-s)).sum());checks['exact_cdf_integral_minus_deficit']=integ-direct
    assert abs(integ-direct)<1e-4
    # Integration-grid convergence: recompute intensive/cross moments from the
    # native source at 25 m, rather than interpolating the existing 50 m grid.
    x,y=z['x_edges'],z['y_edges'];st=Affine(x[1]-x[0],0,x[0],0,-(y[1]-y[0]),y[-1])
    bounds=np.array(transform_bounds('EPSG:4326','EPSG:32647',x[0],y[0],x[-1],y[-1],densify_pts=100))
    left,bottom=np.floor(bounds[:2]/25)*25;right,top=np.ceil(bounds[2:]/25)*25
    mx=np.arange(left,right+1,25);my=np.arange(bottom,top+1,25);dt=Affine(25,0,left,0,-25,top)
    partition=[]
    for sx,sy in [(0,0),(500,0),(0,500),(500,500)]:
        tx=np.arange(np.floor((left-sx)/1000)*1000+sx,np.ceil((right-sx)/1000)*1000+sx+1,1000)
        ty=np.arange(np.floor((bottom-sy)/1000)*1000+sy,np.ceil((top-sy)/1000)*1000+sy+1,1000)
        partition.append((sx,sy,overlaps(mx,tx,'source'),overlaps(my,ty,'source')))
    aggregate={i:{} for i in range(4)}
    for name,v in [('A',a),('S',a*s),('P',p),('H',p*s)]:
        dst=np.zeros((len(my)-1,len(mx)-1),dtype='float64')
        reproject(v[::-1].copy(),dst,src_transform=st,src_crs='EPSG:4326',dst_transform=dt,dst_crs='EPSG:32647',resampling=Resampling.sum,src_nodata=None,dst_nodata=0,num_threads=2)
        # GDAL's sum warp is approximate under nonlinear reprojection. The
        # production 50 m transfer closes to <1e-8 for every stored field.
        # Record the diagnostic 25 m drift without renormalizing it; require
        # <0.01% here and bound its combined effect on the reported estimator.
        drift=float(dst.sum()/v.sum()-1)
        checks['diagnostic_25m_'+name+'_relative_mass_drift']=drift
        if abs(drift)>=1e-4:raise ValueError('25 m diagnostic mass drift exceeds 0.01%')
        for i,(_,_,wx,wy) in enumerate(partition):aggregate[i][name]=(wx@(wy@dst[::-1]).T).T
        print('CHECK 25m',name,flush=True)
        del dst
    convergence=[]
    for i,(sx,sy,_,_) in enumerate(partition):
        q=aggregate[i];f=np.divide(q['S'],q['A'],out=np.zeros_like(q['S']),where=q['A']>0);u=float((q['P']*(1-f)).sum())
        old=uni[(uni.domain=='geographic')&(uni.support=='q30')&(uni.exposure=='population')&(uni.cell_size_m==1000)&(uni.shift_x_m==sx)&(uni.shift_y_m==sy)].estimated_unsupported.iloc[0]
        convergence.append(dict(shift_x=sx,shift_y=sy,integration25=u,integration50=float(old),difference=u-float(old),difference_fraction_total=(u-old)/p.sum()))
    pd.DataFrame(convergence).to_csv(OUT/'allocation/integration_convergence.csv',index=False)
    checks['integration_grid_max_difference_fraction_total']=max(abs(r['difference_fraction_total']) for r in convergence)
    assert checks['integration_grid_max_difference_fraction_total']<.001
    for region in ['masking','dwr/masking']:
        df=pd.read_csv(OUT/region/'replicates.csv');h=pd.read_csv(OUT/region/'synthetic_hazard_bounds.csv')
        assert len(df)==21600 and df.groupby(['mechanism','missing_pixel_fraction','cell_size_fine_pixels','origin_fine_pixels','method']).size().eq(100).all()
        assert ((h.truth>=h.bound_lower-1e-6)&(h.truth<=h.bound_upper+1e-6)).all()
        checks[region+'_replicate_rows']=len(df)
    dwr=pd.read_csv(OUT/'dwr/summary.csv');assert np.allclose(dwr.supported+dwr.unsupported,dwr.total,rtol=1e-12)
    checks['dwr_mass_conservation']='pass'
    pop=pd.read_csv(OUT/'population/model_epoch_comparison.csv');joint=pop[pop.comparison.eq('joint_valid')]
    assert np.allclose(joint.unsupported_common_total,joint.unsupported_share*joint.common_total,rtol=1e-12)
    checks['population_common_total_normalization']='pass'
    dump(OUT/'numerical_checks.json',checks)
    print(json.dumps(checks,indent=2),flush=True)

if __name__=='__main__':main()
