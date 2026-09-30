"""Transparent Windows runner: only I/O portability changes, serial execution.

Upstream scientific routines and scripts are unmodified. GeoTIFF preparation
uses Rasterio separately; unsupported osgeo GeoTIFF helpers must fail if called.
"""
from pathlib import Path
import sys,os,runpy,shutil,difflib
ROOT=Path(__file__).resolve().parent
UPSTREAM=ROOT/'external/LiCSBAS2'
PATCH=ROOT/'external/licsbas_windows_io';PATCH.mkdir(exist_ok=True)
original=(UPSTREAM/'LiCSBAS_lib/LiCSBAS_io_lib.py').read_text(encoding='utf-8')
patched=original.replace('from osgeo import gdal, osr','try:\n    from osgeo import gdal, osr\nexcept ModuleNotFoundError:\n    gdal = osr = None')
patched=patched.replace('gdal.UseExceptions()','if gdal is not None: gdal.UseExceptions()',1)
old="    value = subp.check_output(['grep', field,mlipar]).decode().split()[1].strip()\n    return value"
new="    with open(mlipar, encoding='utf-8') as stream:\n        for line in stream:\n            tokens = line.split()\n            if tokens and tokens[0].rstrip(':') == field:\n                return tokens[1]\n    raise KeyError(field)"
if old not in original:raise RuntimeError('Unexpected upstream get_param_par')
patched=patched.replace(old,new)
(PATCH/'LiCSBAS_io_lib.py').write_text(patched,encoding='utf-8')
(PATCH/'portability.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile='upstream/LiCSBAS_io_lib.py',tofile='windows/LiCSBAS_io_lib.py')),encoding='utf-8')
sys.path[:0]=[str(PATCH),str(UPSTREAM/'LiCSBAS_lib')]
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')

def link_as_copy(source,destination,target_is_directory=False,**kwargs):
    # Upstream uses these links solely for duplicate quick-look PNGs.
    source=Path(source)
    if not source.is_absolute():source=Path(destination).parent/source
    if source.suffix.lower()!='.png':raise RuntimeError('Unexpected non-figure symlink')
    shutil.copy2(source,destination)

os.symlink=link_as_copy
import numpy as np
np.random.seed(20260927)
script=UPSTREAM/'bin'/sys.argv[1]
sys.argv=[str(script)]+sys.argv[2:]
runpy.run_path(str(script),run_name='__main__')
