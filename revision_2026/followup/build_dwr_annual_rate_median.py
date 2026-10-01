import os
from pathlib import Path
import hashlib, json
import numpy as np
import rasterio
from dwr_units import feet_to_feet_per_year,WINDOWS,DAYS_PER_YEAR

BASE = Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent)) / 'dwr_velocity'
inputs = [
    BASE / 'rate_2015_2016.tif',
    BASE / 'annual_oct_2016_2017.tif',
    BASE / 'annual_oct_2017_2018.tif',
    BASE / 'annual_oct_2018_2019.tif',
    BASE / 'annual_oct_2019_2020.tif',
    BASE / 'rate_2020_2021.tif',
]
arrays = []
profile = None
for path in inputs:
    with rasterio.open(path) as src:
        if profile is None:
            profile = src.profile.copy()
            transform, crs = src.transform, src.crs
            shape = (src.height, src.width)
        elif (src.transform != transform or src.crs != crs or (src.height,src.width) != shape):
            raise ValueError(f'Grid mismatch: {path}')
        arr = src.read(1).astype('float64')
        nodata = src.nodata if src.nodata is not None else -9999
        arr[~np.isfinite(arr) | (arr == nodata)] = np.nan
        # The 2015-01 to 2016-01 service was inspected and is zero-filled in Fresno;
        # it is not used. This composite uses six October-to-October annual products.
        arrays.append(feet_to_feet_per_year(arr,path.name))
stack = np.stack(arrays)
complete = np.isfinite(stack).all(axis=0)
out = np.full(shape, -9999, dtype='float32')
out[complete] = np.median(stack[:, complete], axis=0).astype('float32')
profile.update(driver='GTiff', dtype='float32', count=1, nodata=-9999, compress='deflate')
target = BASE / 'rate_median_2015_2021.tif'
with rasterio.open(target, 'w', **profile) as dst:
    dst.write(out, 1)
    dst.update_tags(temporal_statistic='pixelwise median of six October-to-October annual rates',
                    windows='2015-10 to 2016-10; 2016-10 to 2017-10; 2017-10 to 2018-10; 2018-10 to 2019-10; 2019-10 to 2020-10; 2020-10 to 2021-10',
                    units='feet per year', minimum_valid_windows='all six')
metadata = {
    'output': target.name,
    'valid_cells_all_six_windows': int(complete.sum()),
    'total_grid_cells': int(complete.size),
    'median_rate_feet_per_year': float(np.median(out[complete])),
    'median_rate_mm_per_year': float(np.median(out[complete]) * 304.8),
    'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    'output_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
    'algorithm': 'Annualize each raw interval displacement using its exact day count, then take the pixelwise median; all six windows required.',
    'annualization_days_per_year': DAYS_PER_YEAR,
    'input_windows': {p.name:WINDOWS[p.name] for p in inputs},
}
(BASE / 'rate_median_2015_2021_method.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(json.dumps(metadata, indent=2))
