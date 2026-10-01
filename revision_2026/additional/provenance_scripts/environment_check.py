import importlib.util
print({m: bool(importlib.util.find_spec(m)) for m in ['h5py','rasterio','requests','shapefile','networkx','pyproj']})
