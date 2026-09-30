"""Prepare documented small-ROI native LiCSAR inputs for LiCSBAS 11--16."""
import json,shutil
import numpy as np,pandas as pd,rasterio,networkx as nx
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from rasterio.windows import from_bounds,Window
from rasterio.warp import reproject,Resampling
from common import ROOT,OUT,CFG,dump,sha256,resolve_recorded_path

DEST=OUT/'timeseries';DEST.mkdir(exist_ok=True)
ROI=[100.35,13.65,100.65,14.05]

def main():
    manifest=json.loads((ROOT/'external/licsar_expanded/timeseries_manifest.json').read_text(encoding='utf-8'))
    design=pd.read_csv(ROOT/'external/licsar_expanded/experiment_input_plan.csv')
    needed=design[design.timeseries].pair.tolist()
    records={r['pair']:r for r in manifest if r['status']=='ok'}
    if set(needed)!=set(records):raise ValueError('Complete timeseries input transfer required')
    first=resolve_recorded_path(records[sorted(records)[0]]['files']['unw']['path'])
    with rasterio.open(first) as ds:
        win=from_bounds(*ROI,ds.transform)
        c0,r0=int(np.floor(win.col_off)),int(np.floor(win.row_off))
        c1,r1=int(np.ceil(win.col_off+win.width)),int(np.ceil(win.row_off+win.height))
        win=Window(c0,r0,c1-c0,r1-r0);tr=ds.window_transform(win);shape=(r1-r0,c1-c0)
    full=DEST/'GEOC_full';train=DEST/'GEOC_train'
    for p in [full,train]:p.mkdir(exist_ok=True)
    graph=nx.Graph();graph.add_edges_from(p.split('_') for p in sorted(records))
    if not nx.is_connected(graph):raise ValueError('Predeclared time-series period disconnected')
    rng=np.random.default_rng(CFG['seed']);held=[]
    for pair in rng.permutation(sorted(records)):
        a,b=pair.split('_');graph.remove_edge(a,b)
        if nx.is_connected(graph):held.append(str(pair))
        else:graph.add_edge(a,b)
        if len(held)>=int(round(len(records)*.1)):break
    dump(DEST/'holdout_design.json',dict(seed=CFG['seed'],rule='Randomly remove 10% of acquisition-pair edges while retaining global date connectivity; selection uses dates only, before reading pixel quality',training=sorted(set(records)-set(held)),held_out=sorted(held)))
    for i,(pair,record) in enumerate(sorted(records.items())):
        dest=full/pair;dest.mkdir(exist_ok=True)
        for ext,dtype in [('unw','float32'),('cc','uint8')]:
            path=dest/(pair+'.'+ext)
            if not path.exists():
                with rasterio.open(resolve_recorded_path(record['files'][ext]['path'])) as ds:
                    data=np.zeros(shape,dtype=dtype)
                    reproject(rasterio.band(ds,1),data,src_transform=ds.transform,src_crs=ds.crs,dst_transform=tr,dst_crs='EPSG:4326',resampling=Resampling.nearest,src_nodata=0,dst_nodata=0,num_threads=1)
                data[~np.isfinite(data)]=0;data.tofile(path)
            if pair not in held:
                t=train/pair;t.mkdir(exist_ok=True)
                if not (t/path.name).exists():shutil.copy2(path,t/path.name)
        browse=dest/(pair+'.unw.png')
        if not browse.exists():
            phase=np.fromfile(dest/(pair+'.unw'),dtype='float32').reshape(shape)
            wrapped=np.angle(np.exp(1j*phase/3));wrapped[phase==0]=np.nan
            plt.imsave(browse,wrapped,cmap='twilight',vmin=-np.pi,vmax=np.pi)
        if pair not in held and not (train/pair/browse.name).exists():shutil.copy2(browse,train/pair/browse.name)
        if (i+1)%25==0:print('TS PREP',i+1,'/',len(records),flush=True)
    metadata=ROOT/'external/licsar_metadata'
    for name,filename in [('hgt','062D_07629_131313.geo.hgt.tif'),('U.geo','062D_07629_131313.geo.U.tif')]:
        with rasterio.open(metadata/filename) as ds:
            data=np.zeros(shape,dtype='float32')
            reproject(rasterio.band(ds,1),data,src_transform=ds.transform,src_crs=ds.crs,dst_transform=tr,dst_crs='EPSG:4326',resampling=Resampling.bilinear,dst_nodata=0)
        data.tofile(full/name);shutil.copy2(full/name,train/name)
    params=f'range_samples: {shape[1]}\nazimuth_lines: {shape[0]}\nradar_frequency: 5405000000 Hz\n'
    # GAMMA dem_par uses pixel-center registration, unlike GDAL corner affine.
    dem=f'width: {shape[1]}\nnlines: {shape[0]}\ncorner_lat: {tr.f+tr.e/2}\ncorner_lon: {tr.c+tr.a/2}\npost_lat: {tr.e}\npost_lon: {tr.a}\nellipsoid_ra: 6378137.000\nellipsoid_reciprocal_flattening: 298.2572236\n'
    for p in [full,train]:
        (p/'slc.mli.par').write_text(params);(p/'EQA.dem_par').write_text(dem)
        shutil.copy2(metadata/'baselines',p/'baselines')
    dump(DEST/'geometry.json',dict(predeclared_roi=ROI,shape=list(shape),transform=list(tr)[:6],crs='EPSG:4326',pairs=len(records),dates=len(graph),holdout_pairs=len(held),row_order='north-to-south',phase_unit='radians',frequency_hz=5405000000,velocity_sign='LiCSBAS convention: -wavelength/(4*pi)*1000 phase-to-mm, LOS positive toward satellite',registration='Actual source GDAL affine, nearest neighbor onto first unwrapped raster native grid; DEM parameters converted to pixel centers',no_atmospheric_correction='No GACOS/weather correction; infer relative LOS displacement, not validated vertical land motion',upstream_commit='2627e50a1c703becd4a0b18754b05c1cf8c564d5'))
    print('TIMESERIES PREPARED',shape,len(records),flush=True)

if __name__=='__main__':main()
