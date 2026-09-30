"""Fetch optional upstream inputs and reject changed provider artifacts."""
from pathlib import Path
import hashlib
import json
import time
import zipfile
from io import BytesIO


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(url, path, expected_sha256=None):
    import requests
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and (not expected_sha256 or sha256(path) == expected_sha256):
        return
    temporary = path.with_name(path.name + '.download')
    last_error = None
    for attempt in range(3):
        try:
            with requests.get(url, stream=True, timeout=(30, 120)) as response:
                response.raise_for_status()
                with temporary.open('wb') as target:
                    for block in response.iter_content(1024 * 1024):
                        if block:
                            target.write(block)
            if expected_sha256 and sha256(temporary) != expected_sha256:
                raise RuntimeError('Provider file checksum differs from the frozen input: ' + path.name)
            temporary.replace(path)
            print('Downloaded and checked ' + path.name, flush=True)
            return
        except (requests.RequestException, OSError) as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError('Could not retrieve ' + path.name) from last_error


def extract_zip(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            target = (destination / name).resolve()
            if destination.resolve() not in target.parents and target != destination.resolve():
                raise RuntimeError('Unsafe upstream ZIP member: ' + name)
        z.extractall(destination)


def fetch_inputs(root, workspace, dataset):
    import requests
    analysis = workspace / 'analysis'
    strengthening = workspace / 'strengthening'
    if dataset == 'dwr-polygons':
        service='https://gis.water.ca.gov/arcgis/rest/services/Elevation/Vertical_Displacement_Point_Data_Locations_2026Q1/MapServer/0/query'
        out=strengthening/'external/dwr'
        out.mkdir(parents=True,exist_ok=True)
        def query(parameters):
            for attempt in range(3):
                try:
                    response=requests.get(service,params={'f':'json',**parameters},timeout=(20,90))
                    response.raise_for_status()
                    result=response.json()
                    if 'error' in result:raise RuntimeError(str(result['error']))
                    return result
                except requests.RequestException:
                    if attempt==2:raise
                    time.sleep(2*(attempt+1))
        ids=query({'where':'1=1','geometry':'-119.9,36.6,-119.7,36.85','geometryType':'esriGeometryEnvelope',
                   'inSR':4326,'spatialRel':'esriSpatialRelIntersects','returnIdsOnly':'true'})['objectIds']
        features=[]
        for start in range(0,len(ids),800):
            parameters={'objectIds':','.join(map(str,sorted(ids)[start:start+800])),
                        'outFields':'*','returnGeometry':'true','outSR':4326}
            page=out/f'roi_page_{start}.json'
            result=json.loads(page.read_text()) if page.is_file() else query(parameters)
            page.write_text(json.dumps(result),encoding='utf-8')
            features.extend(result['features'])
        if len(features)!=38738 or len({f['attributes']['CODE'] for f in features})!=38738:
            raise RuntimeError('DWR footprint count or unique-code check differs from the retained 38,738-cell domain.')
        (out/'roi_features.json').write_text(json.dumps({'features':features}),encoding='utf-8')
        print('Downloaded 38,738 unique DWR observation polygons.',flush=True)
        return
    if dataset == 'acs':
        manifest = json.loads((root / 'input_manifests/acs_population_source_manifest.json').read_text())
        for item in manifest['sources']:
            download(item['url'], analysis / 'inputs/acs_2021' / item['file'], item['sha256'])
        return
    if dataset == 's20':
        manifest = json.loads((root / 'followup/fresno_population_source_manifest.json').read_text())
        records = {x['file']: x for x in manifest['data']}
        external = strengthening / 'external/round3'
        for name in ['ca2020.pl.zip', 'tl_2020_06019_tabblock20.zip']:
            item = records[name]
            download(item['url'], external / name, item['sha256'])
        extract_zip(external / 'tl_2020_06019_tabblock20.zip', external / 'census_blocks')
        # The provider archive's actual filename includes its GitHub owner prefix.
        response = requests.get('https://zenodo.org/api/records/5874927', timeout=60)
        response.raise_for_status()
        entry = next(x for x in response.json()['files'] if x['key'].endswith('CA-POP-v1_0_0.zip'))
        download(entry['links']['self'], external / 'CA-POP-v1_0_0.zip',
                 records['CA-POP-v1_0_0.zip']['sha256'])
        target = external / 'CAPOP_2020_100m_TOTAL.tif'
        if not target.is_file() or sha256(target) != records[target.name]['sha256']:
            with zipfile.ZipFile(external / 'CA-POP-v1_0_0.zip') as archive:
                names=[x for x in archive.namelist() if x.endswith('CAPOP_2020_100m_TOTAL.tif')]
                if names:
                    container, name = archive, names[0]
                else:
                    outer=next(x for x in archive.namelist() if x.endswith('CAPOP_2020_100m_TOTAL.tif.zip'))
                    container=zipfile.ZipFile(BytesIO(archive.read(outer)))
                    name=next(x for x in container.namelist() if x.endswith('.tif'))
                with container.open(name) as source, target.open('wb') as destination:
                    while block := source.read(1024 * 1024):destination.write(block)
            if sha256(target) != records[target.name]['sha256']:
                raise RuntimeError('CA-POP extracted raster hash mismatch.')
        return
    if dataset == 'census-blocks':
        service = 'https://services.arcgis.com/P3ePLMYs2RVChkJx/ArcGIS/rest/services/USA_Census_2020_DHC_Blocks/FeatureServer/1/query'
        output = analysis / 'census_blocks/fresno_roi'
        output.mkdir(parents=True, exist_ok=True)
        seen = set()
        for offset in range(0, 7000, 1000):
            params = {'where': "GEOID LIKE '06019%'", 'geometry': '-119.9,36.6,-119.7,36.85',
                      'geometryType': 'esriGeometryEnvelope', 'inSR': 4326, 'outSR': 4326,
                      'spatialRel': 'esriSpatialRelIntersects', 'outFields': 'GEOID,P0010001,ALAND,AWATER',
                      'returnGeometry': 'true', 'orderByFields': 'GEOID', 'resultOffset': offset,
                      'resultRecordCount': 1000, 'f': 'geojson'}
            page=output/f'blocks_{offset}.geojson'
            if page.is_file():
                data=json.loads(page.read_text(encoding='utf-8'))
            else:
                for attempt in range(3):
                    try:
                        response=requests.get(service,params=params,timeout=(20,90))
                        response.raise_for_status()
                        data=response.json()
                        break
                    except requests.RequestException:
                        if attempt==2:raise
                        time.sleep(2*(attempt+1))
            features = data.get('features')
            if features is None:
                raise RuntimeError('Census service did not return GeoJSON features.')
            for feature in features:
                gid = feature['properties']['GEOID']
                if gid in seen:
                    raise RuntimeError('Duplicate Census GEOID: ' + gid)
                seen.add(gid)
            (output / f'blocks_{offset}.geojson').write_text(json.dumps(data), encoding='utf-8')
        if len(seen) != 6999:
            raise RuntimeError(f'Upstream Census query changed: expected 6,999 blocks, found {len(seen)}.')
        print('Retrieved 6,999 unique Census block geometries.', flush=True)
        return
    if dataset == 'worldpop-thailand':
        url = 'https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/THA/tha_ppp_2020_UNadj.tif'
        provenance=json.loads((root/'data/strengthening/results/timeseries/worldpop_matched_PROVENANCE.json').read_text())
        expected=next(x['sha256'] for x in provenance['parents'] if x['file']=='worldpop_thailand_2020_UNadj.tif')
        download(url, strengthening / 'external/worldpop_thailand_2020_UNadj.tif',expected)
        return
    raise ValueError(dataset)
