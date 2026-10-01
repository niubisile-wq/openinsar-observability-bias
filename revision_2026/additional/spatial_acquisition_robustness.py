"""Exploratory spatial calibration and shared-acquisition robustness analyses."""
from pathlib import Path
import hashlib
import json
import sys
import os

import networkx as nx
import numpy as np
import pandas as pd
import pyproj
import shapely
from shapely.geometry import shape
from shapely.ops import transform, unary_union

ROOT = Path(__file__).resolve().parent.parent
REVISION = Path(os.environ.get('INSAR_BASELINE_REVISION', ROOT.parent / 'github_original_source/revision_2026'))
OUT = Path(os.environ.get('INSAR_NEW_OUTPUT', ROOT / 'results/spatial_acquisition'))
OUT.mkdir(parents=True, exist_ok=True)
INPUT_HASHES = {}


def record(path):
    INPUT_HASHES[path.relative_to(REVISION).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def save_json(path, value):
    path.write_bytes((json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))


def metrics(prediction, target, moe):
    residual = prediction - target
    return dict(groups=len(target), target_total=float(target.sum()),
                predicted_total=float(prediction.sum()),
                MAE=float(np.mean(np.abs(residual))),
                RMSE=float(np.sqrt(np.mean(residual ** 2))),
                WAPE_pct=float(100 * np.abs(residual).sum() / target.sum()),
                signed_total_error_pct=float(100 * residual.sum() / target.sum()),
                groups_within_acs_90pct_moe=int((np.abs(residual) <= moe).sum()))


def spatial_calibration():
    counts_path = record(REVISION / 'reference_outputs/acs_population_benchmark/blockgroup_counts.csv')
    geometry_path = record(REVISION / 'data/workspace/inputs/acs_compact/fresno_block_groups.geojson')
    values_path = record(REVISION / 'data/workspace/inputs/acs_compact/fresno_acs_values.json')
    counts = pd.read_csv(counts_path, dtype={'GEOID': str}).sort_values('GEOID').reset_index(drop=True)
    counts['GEOID'] = counts.GEOID.str.zfill(12)
    geodata = json.loads(geometry_path.read_text(encoding='utf-8'))
    projection = pyproj.Transformer.from_crs('EPSG:4326', 'EPSG:3310', always_xy=True).transform
    geometries = {feature['properties']['GEOID']: transform(projection, shape(feature['geometry']))
                  for feature in geodata['features']}
    geo = np.asarray([geometries[identifier] for identifier in counts.GEOID], dtype=object)
    centers = np.array([[polygon.centroid.x, polygon.centroid.y] for polygon in geo])
    y = counts.ACS_2017_2021_B01003_estimate.to_numpy(float)
    moe = counts.ACS_2017_2021_B01003_MOE_90pct.to_numpy(float)
    estimates = json.loads(values_path.read_text(encoding='utf-8'))['acs']
    survey_replicates = np.array([estimates[identifier]['replicates'] for identifier in counts.GEOID], dtype=float)
    assert len(counts) == 319 and survey_replicates.shape == (319, 80)
    assert np.array_equal(y, [estimates[identifier]['estimate'] for identifier in counts.GEOID])
    models = {
        'GHSL': counts.GHSL_2020_1km_count_allocated_by_area.to_numpy(float),
        'WorldPop': counts.WorldPop_2020_1km_UNadjusted_count_allocated_by_area.to_numpy(float),
        'Census_2020_P1': counts.census_2020_P1.to_numpy(float),
    }
    raw = []
    for model, prediction in models.items():
        raw.append({'model': model, 'method': 'raw', **metrics(prediction, y, moe)})
        raw.append({'model': model, 'method': 'full_domain_total_calibration_in_sample',
                    **metrics(prediction * y.sum() / prediction.sum(), y, moe)})
    pd.DataFrame(raw).to_csv(OUT / 'spatial_calibration_baselines.csv', index=False)
    rows, assignments, summaries = [], [], []
    for axis, column in [('west_east', 0), ('south_north', 1)]:
        order = np.argsort(centers[:, column], kind='stable')
        folds = np.empty(len(counts), dtype=int)
        for fold, ids in enumerate(np.array_split(order, 5)):
            folds[ids] = fold
        for buffer_m in [0, 1000, 3000]:
            predictions = {model: np.full(len(y), np.nan) for model in models}
            for fold in range(5):
                test = folds == fold
                test_union = unary_union(geo[test])
                distances = shapely.distance(geo, test_union)
                # A 0-m design is unbuffered; buffer designs exclude touching/nearby polygons.
                train = (folds != fold) & ((distances >= 0) if buffer_m == 0 else (distances > buffer_m))
                if int(train.sum()) < 30:
                    raise RuntimeError('Prespecified minimum of 30 training groups not met')
                assert not np.any(train & test)
                if buffer_m:
                    assert np.all(distances[train] > buffer_m)
                for group in range(len(y)):
                    assignments.append({'axis': axis, 'buffer_m': buffer_m, 'heldout_fold': fold,
                                        'GEOID': counts.GEOID.iloc[group],
                                        'role': 'test' if test[group] else ('train' if train[group] else 'buffer_excluded'),
                                        'distance_to_test_polygon_m': float(distances[group])})
                for model, prediction in models.items():
                    factor = float(y[train].sum() / prediction[train].sum())
                    predictions[model][test] = factor * prediction[test]
                    value = metrics(factor * prediction[test], y[test], moe[test])
                    # Uncertainty from the ACS survey design only, conditional on released grid values.
                    factors_rep = survey_replicates[train].sum(axis=0) / prediction[train].sum()
                    errors_rep = (factors_rep * prediction[test].sum() - survey_replicates[test].sum(axis=0))
                    error_point = factor * prediction[test].sum() - y[test].sum()
                    se = float(np.sqrt(4 / 80 * np.sum((errors_rep - error_point) ** 2)))
                    rows.append({'axis': axis, 'buffer_m': buffer_m, 'fold': fold, 'model': model,
                                 'training_groups': int(train.sum()), 'calibration_factor_from_training_only': factor,
                                 'signed_total_error_acs_only_SE': se, **value})
            for model, prediction in predictions.items():
                if not np.isfinite(prediction).all():
                    raise RuntimeError('Each group must be predicted exactly once per design')
                summaries.append({'axis': axis, 'buffer_m': buffer_m, 'model': model,
                                  'method': 'spatially_heldout_total_calibration',
                                  **metrics(prediction, y, moe)})
    pd.DataFrame(rows).to_csv(OUT / 'spatial_calibration_folds.csv', index=False)
    pd.DataFrame(assignments).to_csv(OUT / 'spatial_calibration_membership.csv', index=False)
    pd.DataFrame(summaries).to_csv(OUT / 'spatial_calibration_summary.csv', index=False)
    save_json(OUT / 'spatial_calibration_protocol.json', {
        'status': 'Completed exploratory robustness analysis; holds out recalibration targets, not source-model training',
        'region': 'Same fixed 319 complete Fresno block groups as S24',
        'target': 'ACS 2017-2021 point estimates with published 80 variance replicates',
        'partition': 'Five coordinate-rank geographic bands; west-east and south-north; no outcome used to construct folds',
        'buffers_m': [0, 1000, 3000], 'buffer_rule': 'Exclude training polygons within buffer of heldout polygon union',
        'calibration': 'One multiplier fitted from training groups only; no heldout ACS totals used in that multiplier',
        'uncertainty': '4/80 replicate-squared-difference SE of fold signed error; ACS survey design component only',
        'limits': ['Frozen GHSL/WorldPop models may use related source information; source independence is unresolved.',
                   'This holds out recalibration targets, not the original population-model training pipeline.',
                   'Block-group agreement does not validate household placement or missing-cell deformation.',
                   'Spatial folds and buffers are sensitivity designs, not independent region replicates.',
                   'The target spans 2017-2021 while grids and Census are 2020.'],
        'checks': {'train_test_disjoint': True, 'buffer_exclusions_verified': True,
                   'every_group_predicted_once_per_design': True, 'survey_targets_match_retained_S24_table': True},
    })
    print('ACS spatial recalibration robustness analysis completed', flush=True)


def acquisition_influence():
    path = record(REVISION / 'data/strengthening/results/sampling/expanded_pair_exposure.csv')
    data = pd.read_csv(path, dtype={'pair': str})
    unique = data.drop_duplicates('pair').copy()
    unique['first_date'] = unique.pair.str[:8]
    unique['last_date'] = unique.pair.str[9:]
    all_pairs = unique.pair.tolist()
    graph = nx.Graph(pair.split('_') for pair in all_pairs)
    dates = sorted(graph.nodes)
    assert len(all_pairs) == 224 and len(dates) == 225
    rows = []
    group_rows = list(data.groupby(['support', 'exposure']))
    for date in dates:
        keep_pairs = unique.loc[(unique.first_date != date) & (unique.last_date != date), 'pair'].tolist()
        reduced = nx.Graph(pair.split('_') for pair in keep_pairs)
        for (support, exposure), part in group_rows:
            part = part.copy()
            part['share'] = part.unsupported / part.total
            retained = part[part.pair.isin(keep_pairs)]
            baseline_mean = float(part.share.mean())
            baseline_stratified = float(part.groupby(['year', 'quarter']).share.mean().mean())
            original_strata = set(zip(part.year, part.quarter))
            retained_strata = set(zip(retained.year, retained.quarter))
            if original_strata != retained_strata:
                raise RuntimeError('Date deletion erased a year-quarter stratum; target would change')
            assert len(original_strata) == 28
            share_mean = float(retained.share.mean())
            share_stratified = float(retained.groupby(['year', 'quarter']).share.mean().mean())
            rows.append({'deleted_acquisition_date': date, 'support': support, 'exposure': exposure,
                         'removed_pairs': 224 - len(keep_pairs), 'retained_pairs': len(keep_pairs),
                         'baseline_deficit_share_pct': 100 * baseline_mean,
                         'remaining_equal_pair_share_pct': 100 * share_mean,
                         'equal_pair_change_pp': 100 * (share_mean - baseline_mean),
                         'remaining_equal_quarter_share_pct': 100 * share_stratified,
                         'equal_quarter_change_pp': 100 * (share_stratified - baseline_stratified),
                         'retained_year_quarter_strata': len(retained_strata),
                         'retained_dates': len(reduced.nodes),
                         'retained_components': nx.number_connected_components(reduced),
                         'retained_cycle_rank': len(reduced.edges) - len(reduced.nodes) + nx.number_connected_components(reduced)})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / 'acquisition_date_influence.csv', index=False)
    summary = []
    for (support, exposure), part in table.groupby(['support', 'exposure']):
        summary.append({'support': support, 'exposure': exposure,
                        'deleted_dates': len(part),
                        'max_pairs_removed': int(part.removed_pairs.max()),
                        'equal_pair_min_change_pp': float(part.equal_pair_change_pp.min()),
                        'equal_pair_max_change_pp': float(part.equal_pair_change_pp.max()),
                        'equal_quarter_min_change_pp': float(part.equal_quarter_change_pp.min()),
                        'equal_quarter_max_change_pp': float(part.equal_quarter_change_pp.max()),
                        'retained_components_min': int(part.retained_components.min()),
                        'retained_components_max': int(part.retained_components.max())})
    pd.DataFrame(summary).to_csv(OUT / 'acquisition_date_influence_summary.csv', index=False)
    save_json(OUT / 'acquisition_date_protocol.json', {
        'status': 'Completed finite-pool shared-acquisition sensitivity analysis',
        'baseline': {'pairs': 224, 'dates': 225, 'components': nx.number_connected_components(graph),
                     'year_quarter_strata': 28},
        'deletion_unit': 'One acquisition date and every interferogram incident to it',
        'estimands': ['Equal retained-pair mean', 'Equal original year-quarter stratum mean after date deletion'],
        'checks': {'all_original_strata_retained': True, 'node_incident_edges_deleted_together': True},
        'limits': ['Finite-pool support sensitivity only; no velocity time-series refit.',
                   'No iid pair or date assumption and no jackknife confidence interval.',
                   'The 224-pair sample graph is disconnected; this analysis does not establish processing connectivity.',
                   'Acquisition deletion does not supply new observations or external physical validation.'],
    })
    print('Shared-acquisition influence analysis completed', flush=True)


if __name__ == '__main__':
    spatial_calibration()
    acquisition_influence()
    save_json(OUT / 'robustness_provenance.json', {
        'status': 'Completed exploratory robustness analyses; no independent source-model validation claim',
        'source_scientific_release': 'v0.4.1',
        'inputs_sha256': INPUT_HASHES,
        'environment': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__,
                        'shapely': shapely.__version__, 'pyproj': pyproj.__version__, 'networkx': nx.__version__},
    })
