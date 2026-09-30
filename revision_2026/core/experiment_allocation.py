"""E1: conservative same-target allocation across metric scales and origins.

GHSL is an allocation model, not observed household locations. Preserve the
exposure-support cross moment when projecting; do not manufacture a new fine
reference by multiplying separately projected averages.
"""
import json
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio.warp import reproject,transform_bounds,transform,Resampling
from common import OUT,CFG,dump,overlaps

DEST=OUT/'allocation';DEST.mkdir(exist_ok=True)
METRIC=DEST/'metric50';METRIC.mkdir(exist_ok=True)

def geometry(g):
    x,y=g['x_edges'],g['y_edges']
    st=Affine(x[1]-x[0],0,x[0],0,-(y[1]-y[0]),y[-1])
    left,bottom,right,top=transform_bounds('EPSG:4326','EPSG:32647',x[0],y[0],x[-1],y[-1],densify_pts=100)
    left,bottom=np.floor(np.array([left,bottom])/50)*50
    right,top=np.ceil(np.array([right,top])/50)*50
    w,h=int(round((right-left)/50)),int(round((top-bottom)/50))
    return st,Affine(50,0,left,0,-50,top),left+np.arange(w+1)*50,bottom+np.arange(h+1)*50

def project(a,name,st,dt,x,y):
    path=METRIC/(name+'.npy')
    if not path.exists():
        dst=np.zeros((len(y)-1,len(x)-1),dtype='float64')
        reproject(a[::-1].copy(),dst,src_transform=st,src_crs='EPSG:4326',dst_transform=dt,dst_crs='EPSG:32647',resampling=Resampling.sum,src_nodata=None,dst_nodata=0,num_threads=2)
        discrepancy=float(dst.sum()-a.sum())
        if abs(discrepancy)>max(1e-5,abs(a.sum())*1e-8):raise ValueError(f'Conservation {name} {discrepancy}')
        np.save(path,dst[::-1])
    return np.load(path,mmap_mode='r')

def partitions(x,y,size,sx,sy):
    tx=np.arange(np.floor((x[0]-sx)/size)*size+sx,np.ceil((x[-1]-sx)/size)*size+sx+size*.1,size)
    ty=np.arange(np.floor((y[0]-sy)/size)*size+sy,np.ceil((y[-1]-sy)/size)*size+sy+size*.1,size)
    wx=overlaps(x,tx,'source');wy=overlaps(y,ty,'source')
    return tx,ty,lambda a:(wx@(wy@a).T).T

def center_support(f,d,g,tx,ty):
    xx,yy=np.meshgrid((tx[:-1]+tx[1:])/2,(ty[:-1]+ty[1:])/2)
    lon,lat=transform('EPSG:32647','EPSG:4326',xx.ravel(),yy.ravel())
    ix=np.searchsorted(g['x_edges'],lon,side='right')-1
    iy=np.searchsorted(g['y_edges'],lat,side='right')-1
    valid=(ix>=0)&(ix<f.shape[1])&(iy>=0)&(iy<f.shape[0])
    result=np.zeros(ix.shape,dtype='float64')
    ids=np.where(valid)[0];inside=d[iy[ids],ix[ids]]
    ids=ids[inside];result[ids]=f[iy[ids],ix[ids]]
    return result.reshape(xx.shape)

def main():
    with np.load(OUT/'fine/grid.npz') as z:g={k:z[k] for k in z.files}
    with np.load(OUT/'fine/mean_support.npz') as z:fields={k:z[k].astype('float64') for k in z.files}
    st,dt,x,y=geometry(g);rows=[];ledger=[]
    np.savez(DEST/'metric_geometry.npz',x_edges=x,y_edges=y)
    for domain,key in [('geographic','domain'),('inherited_strong','screen')]:
        d=g[key];a=g['area_m2']*d
        ma=project(a,domain+'_area',st,dt,x,y)
        for support,f in fields.items():
            ms=project(a*f,domain+'_'+support+'_supported_area',st,dt,x,y)
            for exposure in ['population','builtup_m2']:
                p=g[exposure]*d
                mp=project(p,domain+'_'+exposure,st,dt,x,y)
                mh=project(p*f,domain+'_'+support+'_'+exposure+'_cross',st,dt,x,y)
                fine_total=float(p.sum());fine_u=float((p*(1-f)).sum())
                ledger.append(dict(domain=domain,support=support,exposure=exposure,native_total=fine_total,projected_total=float(mp.sum()),native_unsupported=fine_u,projected_unsupported=float(mp.sum()-mh.sum())))
                for size in CFG['metric_grid_sizes_m']:
                    for sx,sy in [(0,0),(size/2,0),(0,size/2),(size/2,size/2)]:
                        tx,ty,aggregate=partitions(x,y,size,sx,sy)
                        A,S,P,H=[aggregate(z) for z in [ma,ms,mp,mh]]
                        keep=A>1e-6;support_cell=np.divide(S,A,out=np.zeros_like(S),where=keep)
                        reference=P-H;uniform=P*(1-support_cell)
                        center=P*(1-center_support(f,d,g,tx,ty))
                        covariance_integral=H-P*support_cell
                        if not np.allclose(uniform-reference,covariance_integral,rtol=1e-9,atol=1e-7):raise ValueError('Covariance identity')
                        if not np.isclose(reference.sum(),fine_u,rtol=1e-8,atol=1e-5):raise ValueError('Reference changed with grid')
                        for method,estimate in [('area_uniform',uniform),('cell_center',center)]:
                            err=(estimate-reference)[keep];ref=reference[keep]
                            rows.append(dict(domain=domain,support=support,exposure=exposure,cell_size_m=size,shift_x_m=sx,shift_y_m=sy,method=method,cells=int(keep.sum()),total_exposure=fine_total,reference_unsupported=fine_u,estimated_unsupported=float(estimate.sum()),signed_error=float(err.sum()),error_fraction_total=float(err.sum()/fine_total),error_fraction_reference=float(err.sum()/fine_u),cell_mae=float(np.abs(err).mean()),cell_rmse=float(np.sqrt(np.mean(err**2))),cell_absolute_error_sum=float(np.abs(err).sum()),covariance_integral=float(covariance_integral.sum())))
                        if domain=='geographic' and support=='q30' and exposure=='population' and size==1000 and sx==sy==0:
                            yy,xx=np.meshgrid((ty[:-1]+ty[1:])/2,(tx[:-1]+tx[1:])/2,indexing='ij')
                            pd.DataFrame(dict(easting=xx[keep],northing=yy[keep],area=A[keep],population=P[keep],support=support_cell[keep],reference=reference[keep],uniform=uniform[keep],center=center[keep],covariance_integral=covariance_integral[keep])).to_csv(DEST/'primary_1000m_cells.csv',index=False)
                print(domain,support,exposure,'complete',flush=True)
    pd.DataFrame(rows).to_csv(DEST/'scale_origin_comparison.csv',index=False)
    pd.DataFrame(ledger).to_csv(DEST/'conservation.csv',index=False)
    dump(DEST/'method.json',dict(crs='EPSG:32647',integration_grid_m=50,scales_m=CFG['metric_grid_sizes_m'],origins='four combinations of zero and half-cell shift',transfer='GDAL sum conservative projection of extensive first and cross moments, then exact rectangular mass splitting on the projected integration grid',fine_reference='GHSL 3ss within-pixel uniform population times exact spherical support; an allocation-model reference, not ground truth',area_uniform='Redistribute the same cell exposure uniformly over its included-domain area',cell_center='Use support of the fine pixel containing the metric cell center; centers outside the fixed domain have zero support',identity='U_uniform - U_fine = sum_c [H_c - P_c S_c/A_c] = sum_c A_c Cov_area(density,support)',boundary='Retain partial domain cells and their exact projected extensive weights; no full-cell population extrapolation'))
    print('E1 COMPLETE',flush=True)

if __name__=='__main__':main()
