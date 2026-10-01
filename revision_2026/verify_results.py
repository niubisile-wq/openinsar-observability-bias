"""Compare recomputed scientific tables with the frozen manuscript outputs."""
from pathlib import Path
import json

# Each pair is (runtime relative path, release reference relative path).
TABLES = {
    'allocation': [('strengthening/results/allocation/' + x, 'core/allocation/' + x)
                   for x in ['scale_origin_comparison.csv', 'conservation.csv', 'primary_1000m_cells.csv']],
    'synthetic-masking': [('strengthening/results/masking/' + x, 'core/masking/' + x)
                          for x in ['replicates.csv', 'summary.csv', 'synthetic_hazard_bounds.csv']],
    'dwr-masking': [('strengthening/results/dwr/masking/' + x, 'core/dwr/masking/' + x)
                    for x in ['replicates.csv', 'summary.csv', 'synthetic_hazard_bounds.csv']],
    'dwr-allocation': [('strengthening/results/dwr/' + x, 'core/dwr/' + x)
                      for x in ['summary.csv', 'scale_origin_comparison.csv']],
    'population-analysis': [('strengthening/results/population/' + x, 'core/population/' + x)
                            for x in ['model_epoch_comparison.csv', 'model_curves.csv', 'landcover_strata.csv']],
    'sampling-analysis': [('strengthening/results/sampling/' + x, 'core/sampling/' + x)
                          for x in ['nested_summary.csv', 'nested_curves.csv', 'temporal_robustness.csv',
                                    'conditional_sampling_replicates.csv']],
    'sampling-spatial': [('strengthening/results/sampling/' + x, 'core/sampling/' + x)
                        for x in ['spatial_blocks.csv','spatial_rank_stability.csv']],
    'cross-sensitivity': [('strengthening/results/cross_sensitivity/' + x, 'core/cross_sensitivity/' + x)
                          for x in ['cross_sensitivity.csv', 'integration_convergence.csv', 'conservation.csv']],
    'paired': [('analysis/paired_support_pilot/' + x, 'paired_support_pilot/' + x)
               for x in ['all_cases.csv', 'primary_cases.csv', 'convergence.csv']],
    'direct-geometry': [('analysis/direct_geometry_audit/' + x, 'direct_geometry_audit/' + x)
                        for x in ['summary.csv', 'cell_intersections.csv']],
    'observed-los-masking': [('analysis/actual_los_masking/' + x, 'actual_los_masking/' + x)
                             for x in ['replicates.csv', 'summary.csv']],
    'partial-identification': [('analysis/partial_identification/bounds.csv', 'partial_identification/bounds.csv')],
    'published-method-benchmark': [('analysis/published_method_benchmark/' + x, 'published_method_benchmark/' + x)
                                   for x in ['replicates.csv', 'summary.csv']],
    'reporting-tolerance': [('analysis/reporting_tolerance/tolerance_summary.csv',
                            'reporting_tolerance/tolerance_summary.csv')],
    'acs-population-benchmark': [('analysis/acs_population_benchmark/' + x, 'acs_population_benchmark/' + x)
                                 for x in ['blockgroup_counts.csv', 'model_summary.csv']],
    'fine-population': [('analysis/fresno_population_allocation/' + x,
                         'followup/fresno_population_allocation/' + x)
                        for x in ['per_block.csv', 'method_summary.csv', 'scale_method_sensitivity.csv',
                                  'tract_decomposition.csv']],
    'population-benchmark': [('analysis/fresno_population_benchmark/fresno_census_block_population_model_metrics.csv',
                              'fresno_population_benchmark/fresno_census_block_population_model_metrics.csv')],
    'followup': [('analysis/analysis_results/' + x, 'analysis_results/' + x)
                  for x in ['allocation_covariance_bounds.csv',
                            'population_model_common_domain.csv',
                            'block_median_baseline_summary.csv']],
    'spatial': [('analysis/bangkok_spatial_decomposition/bangkok_spatial_effect_summary.csv',
                'bangkok_spatial_decomposition/bangkok_spatial_effect_summary.csv')],
    'census-support': [('analysis/census_block_analysis/fresno_block_support_summary.csv',
                       'census_block_analysis/fresno_block_support_summary.csv')],
    'census-rates': [('analysis/census_block_velocity_analysis/fresno_census_block_velocity_summary.csv',
                      'census_block_velocity_analysis/fresno_census_block_velocity_summary.csv')],
    'temporal': [('analysis/census_block_velocity_analysis/fresno_population_temporal_stability_summary.csv',
                  'census_block_velocity_analysis/fresno_population_temporal_stability_summary.csv')],
    'spatial-dependence': [],
}


ADDITIONAL_TABLES = {
    'additional-real-blocks': ('real_block_masking', ['real_block_unmasked_baseline.csv', 'real_block_masking_replicates.csv', 'real_block_masking_summary.csv']),
    'additional-dependence': ('spatial_acquisition', ['spatial_calibration_baselines.csv', 'spatial_calibration_folds.csv', 'spatial_calibration_membership.csv', 'spatial_calibration_summary.csv', 'acquisition_date_influence.csv', 'acquisition_date_influence_summary.csv']),
    'additional-quality': ('training_quality', ['quality_coverage_error.csv', 'quality_coverage_error_by_pair.csv']),
    'additional-regions': ('new_regions', ['new_region_support_baselines.csv', 'new_region_allocation_summary.csv', 'new_region_allocation_by_tile.csv']),
    'additional-residential': ('residential_allocation', ['residential_allocation_summary.csv', 'residential_allocation_by_block.csv', 'residential_landuse_eligibility.csv']),
    'additional-gnss': ('gnss_transfer', ['gnss_transfer_summary.csv', 'gnss_transfer_by_station.csv', 'gnss_transfer_exclusions.csv', 'gnss_transfer_collocation_clusters.csv']),
}
for mode, (folder, names) in ADDITIONAL_TABLES.items():
    TABLES[mode] = [('additional/'+folder+'/'+name, 'additional/'+folder+'/'+name) for name in names]


def verify_results(root, workspace, modes=None):
    import numpy as np
    import pandas as pd
    from pandas.api.types import is_numeric_dtype
    rows = []
    selected = list(TABLES) if modes is None else modes
    for mode in selected:
        for actual, expected in TABLES.get(mode, []):
            actual_path, expected_path = workspace / actual, root / 'reference_outputs' / expected
            if modes is None and not actual_path.exists():
                continue
            if not actual_path.is_file() or not expected_path.is_file():
                raise FileNotFoundError(f'Missing comparison input for {mode}: {actual_path} / {expected_path}')
            dtype = {'GEOID': str, 'GEOID20': str, 'tract_id': str}
            got, ref = pd.read_csv(actual_path, dtype=dtype), pd.read_csv(expected_path, dtype=dtype)
            # Raw totals/areas allow tiny GDAL/NumPy rounding differences;
            # fractions and proportions retain a stricter absolute tolerance.
            pd.testing.assert_frame_equal(got, ref, check_dtype=False, check_exact=False,
                                          rtol=1e-8, atol=1e-6)
            for column in ref:
                if is_numeric_dtype(ref[column]) and any(token in column.lower() for token in
                        ['share', 'fraction', '_pct', '_pp', 'spearman', 'correlation']):
                    np.testing.assert_allclose(got[column].to_numpy(float), ref[column].to_numpy(float),
                                               rtol=1e-8, atol=1e-9, equal_nan=True,
                                               err_msg=f'{mode}: {column}')
            max_error = 0.0
            for column in ref:
                if is_numeric_dtype(ref[column]) and ref[column].dtype.kind != 'b':
                    delta = np.abs(got[column].to_numpy(float) - ref[column].to_numpy(float))
                    if np.any(np.isfinite(delta)):
                        max_error = max(max_error, float(np.nanmax(delta)))
            rows.append({'mode': mode, 'table': expected, 'rows': len(ref), 'columns': len(ref.columns),
                         'max_absolute_numeric_difference': max_error, 'pass': True})
    if 'spatial-dependence' in selected:
        got = json.loads((workspace / 'analysis/fresno_population_allocation/spatial_autocorrelation.json').read_text())
        ref = json.loads((root / 'reference_outputs/followup/fresno_population_allocation/spatial_autocorrelation.json').read_text())
        for key in ['n_blocks', 'undirected_neighbor_pairs', 'isolates']:
            assert got[key] == ref[key], key
        for key, value in ref['fields'].items():
            assert np.isclose(got['fields'][key], value, rtol=1e-10, atol=1e-12), key
        rows.append({'mode': 'spatial-dependence', 'table': 'spatial_autocorrelation.json', 'pass': True})
    if not rows:
        raise RuntimeError('No recomputed tables were found; prepare/check-code alone does not reproduce an analysis.')
    return {'status': 'PASS', 'relative_tolerance': 1e-8, 'absolute_tolerance': 1e-6,
            'fraction_absolute_tolerance': 1e-9,
            'scope': 'Scientific computation from frozen analysis inputs; not independent physical validation.',
            'tables': rows}
