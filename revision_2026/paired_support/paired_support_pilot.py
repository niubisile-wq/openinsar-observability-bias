import os
"""Exploratory paired-support allocation study; no manuscript changes.

Uses the same 2019-2020 Bangkok processing grid, GHSL allocation and AOI for
the 437-input-pair mean support and final LiCSBAS validity mask. Preserves
A, P, A*s and P*s separately. Neither field is independent physical truth.
"""
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import pandas as pd
import rasterio
import scipy
from affine import Affine
from rasterio.warp import Resampling, reproject, transform, transform_bounds
from scipy.sparse import coo_matrix

HERE = Path(os.environ.get('INSAR_PAIRED_ROOT', Path(__file__).resolve().parent))
PROJECT = HERE.parent
WORK = Path(os.environ['INSAR_STRENGTHENING_ROOT'])
INPUT = WORK / 'results/timeseries/matched_quality.npz'
GEOMETRY = WORK / 'results/timeseries/geometry.json'
DEST = HERE / 'paired_support_pilot'
SCALES = [250, 500, 1000, 2000]
RESOLUTIONS = [50, 25, 12.5]
CRS = 'EPSG:32647'


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding='utf-8')


def overlap(source, target):
    """Source-fraction interval overlaps. Columns must sum to one."""
    rr, cc, vv = [], [], []
    i = j = 0
    while i < len(source)-1 and j < len(target)-1:
        lo = max(source[i], target[j]); hi = min(source[i+1], target[j+1])
        if hi > lo:
            rr.append(j); cc.append(i); vv.append((hi-lo)/(source[i+1]-source[i]))
        if source[i+1] <= target[j+1]:
            i += 1
        else:
            j += 1
    w = coo_matrix((vv, (rr, cc)), shape=(len(target)-1, len(source)-1)).tocsr()
    if not np.allclose(np.asarray(w.sum(axis=0)), 1, rtol=0, atol=1e-12):
        raise ValueError('Incomplete partition')
    return w


def partition(x, y, size, sx, sy):
    tx = np.arange(np.floor((x[0]-sx)/size)*size+sx,
                   np.ceil((x[-1]-sx)/size)*size+sx+size*.1, size)
    ty = np.arange(np.floor((y[0]-sy)/size)*size+sy,
                   np.ceil((y[-1]-sy)/size)*size+sy+size*.1, size)
    wx, wy = overlap(x, tx), overlap(y, ty)
    return tx, ty, lambda a: (wx @ (wy @ a).T).T


def center_values(field, source_x, source_y, tx, ty):
    xx, yy = np.meshgrid((tx[:-1]+tx[1:])/2, (ty[:-1]+ty[1:])/2)
    lon, lat = transform(CRS, 'EPSG:4326', xx.ravel(), yy.ravel())
    c = np.searchsorted(source_x, lon, side='right')-1
    r = np.searchsorted(source_y, lat, side='right')-1
    good = (c >= 0) & (c < field.shape[1]) & (r >= 0) & (r < field.shape[0])
    result = np.zeros(c.size)
    result[good] = field[field.shape[0]-1-r[good], c[good]]
    return result.reshape(xx.shape)


def main():
    DEST.mkdir(exist_ok=True)
    protocol = dict(status='exploratory, not integrated into manuscript',
        input=str(INPUT.relative_to(PROJECT)),
        sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        fields=['matched437_mean_support', 'final_quality_mask'],
        scales_m=SCALES, integration_m=RESOLUTIONS, primary_integration_m=25,
        refinement='Initial 50/25-m check failed at a 250-m half-grid origin (0.055095 pp). All cases retained; refine to 25/12.5 m under the unchanged 0.01-pp gate. See initial_50m_gate_failure.json.',
        origins='four combinations of 0 and half coarse-cell size',
        methods=['uniform_area', 'center'],
        reference='Fixed GHSL allocation on the processing grid; not household truth',
        primary_metric='signed difference in percentage points of fixed total population',
        secondary_metrics=['difference relative to native deficit', 'sum absolute cell difference'],
        projection_max_relative_drift=1e-4, convergence_max_pp_total=.01,
        center_rule='source pixel containing nominal coarse-cell center; zero outside domain',
        boundary='Partial cells retain their included area and population; no full-cell extrapolation',
        limitation='The two support definitions estimate different quantities; neither supplies hidden-motion validation',
        dependencies=dict(python=platform.python_version(), numpy=np.__version__,
                          rasterio=rasterio.__version__, scipy=scipy.__version__))
    write_json(DEST/'protocol.json', protocol)
    geom = json.loads(GEOMETRY.read_text(encoding='utf-8'))
    st = Affine(*geom['transform'])
    with np.load(INPUT) as z:
        source_x, source_y = z['x_edges'], z['y_edges']
        p, a = z['population'].astype(float), z['area_m2'].astype(float)
        fields = {'matched437_mean_support':z['matched_support'].astype(float),
                  'final_quality_mask':z['full_valid'].astype(float)}
    if not (np.isfinite(p).all() and np.isfinite(a).all() and (p >= 0).all() and (a > 0).all()):
        raise ValueError('Invalid extensive input')
    for f in fields.values():
        if not (np.isfinite(f).all() and (f >= 0).all() and (f <= 1).all()):
            raise ValueError('Invalid support')
    if not np.allclose([st.c, st.f], [source_x[0], source_y[-1]], atol=1e-10, rtol=0):
        raise ValueError('Source geometry mismatch')
    total = float(p.sum())
    left, bottom, right, top = transform_bounds('EPSG:4326', CRS,
        source_x[0], source_y[0], source_x[-1], source_y[-1], densify_pts=100)
    left, bottom = np.floor(np.array([left, bottom])/50)*50
    right, top = np.ceil(np.array([right, top])/50)*50
    rows, ledger = [], []
    for resolution in RESOLUTIONS:
        x = left+np.arange(round((right-left)/resolution)+1)*resolution
        y = bottom+np.arange(round((top-bottom)/resolution)+1)*resolution
        dt = Affine(resolution, 0, left, 0, -resolution, top)

        def project(value, label):
            out = np.zeros((len(y)-1, len(x)-1), dtype='float64')
            reproject(value, out, src_transform=st, src_crs='EPSG:4326',
                dst_transform=dt, dst_crs=CRS, resampling=Resampling.sum,
                src_nodata=None, dst_nodata=0, num_threads=2)
            expected, actual = float(value.sum()), float(out.sum())
            drift = (actual-expected)/expected if expected else actual
            if abs(drift) > protocol['projection_max_relative_drift']:
                raise ValueError(('Projection drift', label, drift))
            ledger.append(dict(resolution_m=resolution, quantity=label,
                source_sum=expected, projected_sum=actual, relative_drift=drift))
            return out[::-1]

        ma, mp = project(a, 'area'), project(p, 'population')
        partitions = []
        for size in SCALES:
            for sx, sy in [(0,0), (size/2,0), (0,size/2), (size/2,size/2)]:
                tx, ty, agg = partition(x,y,size,sx,sy)
                A, P = agg(ma), agg(mp)
                if not np.isclose(P.sum(), mp.sum(), atol=1e-6, rtol=1e-12):
                    raise ValueError('Population not conserved in partition')
                partitions.append((size,sx,sy,tx,ty,agg,A,P))
        for name, f in fields.items():
            ms, mh = project(a*f, name+'_area_support'), project(p*f, name+'_population_support')
            native = float((p*(1-f)).sum())
            for size,sx,sy,tx,ty,agg,A,P in partitions:
                S, H = agg(ms), agg(mh)
                keep = A > 1e-6
                sbar = np.divide(S,A,out=np.zeros_like(S),where=keep)
                reference = P-H
                uniform = P*(1-sbar)
                covariance = H-P*sbar
                if not np.allclose(uniform-reference,covariance,atol=1e-7,rtol=1e-9):
                    raise ValueError('Covariance check')
                center = P*(1-center_values(f,source_x,source_y,tx,ty))
                for method, estimate in [('uniform_area',uniform),('center',center)]:
                    difference = estimate-reference
                    rows.append(dict(integration_m=resolution, support=name, size_m=size,
                        origin_x_m=sx, origin_y_m=sy, method=method,
                        population=total, source_native_deficit=native,
                        source_native_deficit_pct=100*native/total,
                        projected_native_deficit=float(reference.sum()),
                        estimated_deficit=float(estimate.sum()),
                        signed_difference=float(difference.sum()),
                        difference_pp_total=100*float(difference.sum())/total,
                        difference_pct_native=100*float(difference.sum())/native,
                        cell_absolute_difference_pp_total=100*float(np.abs(difference[keep]).sum())/total))
            print('complete', resolution, name, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(DEST/'all_cases.csv',index=False)
    pd.DataFrame(ledger).to_csv(DEST/'projection_ledger.csv',index=False)
    primary = df[df.integration_m.eq(25)].copy()
    primary.to_csv(DEST/'primary_cases.csv',index=False)
    keys=['support','size_m','origin_x_m','origin_y_m','method']
    initial=df[df.integration_m.eq(50)].merge(primary,on=keys,suffixes=('_50','_25'))
    initial['integration_difference_pp_total']=initial.difference_pp_total_25-initial.difference_pp_total_50
    initial.to_csv(DEST/'initial_50m_convergence.csv',index=False)
    compare=primary.merge(df[df.integration_m.eq(12.5)],on=keys,suffixes=('_25','_12_5'))
    compare['integration_difference_pp_total']=compare.difference_pp_total_12_5-compare.difference_pp_total_25
    compare.to_csv(DEST/'convergence.csv',index=False)
    maxdelta=float(compare.integration_difference_pp_total.abs().max())
    report=[]
    for (name,method), d in primary[primary.size_m.eq(1000)].groupby(['support','method']):
        report.append(dict(support=name,method=method,native_deficit_pct=float(d.source_native_deficit_pct.iloc[0]),
            difference_pp_total_min=float(d.difference_pp_total.min()),
            difference_pp_total_max=float(d.difference_pp_total.max()),
            relative_difference_pct_min=float(d.difference_pct_native.min()),
            relative_difference_pct_max=float(d.difference_pct_native.max())))
    audit=dict(primary_cases=len(primary),input_population=total,shape=list(p.shape),
        max_projection_relative_drift=max(abs(r['relative_drift']) for r in ledger),
        max_integration_discrepancy_pp_total=maxdelta,summary_1000m=report,
        convergence_pass=bool(maxdelta <= protocol['convergence_max_pp_total']),
        primary_integration_m=25, diagnostic_integration_m=12.5,
        status='pilot computed; direct-source intersection audit and manuscript integration not performed',
        interpretation=protocol['limitation'])
    write_json(DEST/'audit.json',audit)
    print(json.dumps(audit,indent=2))
    if not audit['convergence_pass']:
        raise ValueError(('Integration discrepancy',maxdelta))


if __name__ == '__main__':
    main()
