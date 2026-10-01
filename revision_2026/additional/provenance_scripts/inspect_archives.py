from pathlib import Path
import zipfile
ROOT=Path(__file__).resolve().parent.parent.parent
for path in (ROOT/'13_GitHub复现发布_20261001').glob('reviewer-*.zip'):
    with zipfile.ZipFile(path) as z:
        print(path.name, [(v.filename,v.file_size) for v in z.infolist() if 'dwr_velocity' in v.filename or 'pixels_by_block' in v.filename],flush=True)
