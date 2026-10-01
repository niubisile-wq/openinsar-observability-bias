"""Reproduce the 2026 revision from a portable, hash-checked checkout."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
CORE = {
    'allocation': ('experiment_allocation.py', []),
    'synthetic-masking': ('experiment_masking.py', []),
    'dwr-masking': ('experiment_masking.py', ['--region', 'dwr']),
    'population-analysis': ('experiment_population.py', ['--frozen-models']),
    'sampling-analysis': ('experiment_sampling.py', ['--frozen-summaries']),
    'sampling-spatial': ('sampling_spatial.py', []),
    'cross-sensitivity': ('experiment_cross_sensitivity.py', []),
    'dwr-allocation': ('experiment_dwr.py', ['--frozen-grid']),
}
FOLLOWUP = {
    'followup': 'followup_quality.py',
    'census-support': 'fresno_census_block_analysis.py',
    'annual-median': 'build_dwr_annual_rate_median.py',
    'census-rates': 'dwr_census_velocity_analysis.py',
    'temporal': 'fresno_temporal_class_stability.py',
    'persistence': 'fresno_six_year_persistence.py',
    'spatial': 'bangkok_spatial_effect_decomposition.py',
    'population-benchmark': 'fresno_census_population_model_benchmark.py',
    'fine-population': 'fresno_capop_support_sensitivity.py',
    'spatial-dependence': 'spatial_dependence_audit.py',
    'partial-identification': 'partial_identification.py',
    'published-method-benchmark': 'published_method_benchmark.py',
    'reporting-tolerance': 'reporting_tolerance_audit.py',
    'acs-population-benchmark': 'acs_population_benchmark.py',
}
PAIRED = {
    'paired': 'paired_support_pilot.py',
    'direct-geometry': 'direct_geometry_audit.py',
    'observed-los-masking': 'actual_los_masking.py',
}
ADDITIONAL = {
    'additional-real-blocks': ('real_block_masking.py', 'fresno_real_blocks/real_block_masking_inputs.npz', 'real_block_masking'),
    'additional-dependence': ('spatial_acquisition_robustness.py', None, 'spatial_acquisition'),
    'additional-quality': ('quality_coverage_error.py', 'training_quality/quality_curve_inputs.npz', 'training_quality'),
    'additional-regions': ('new_region_allocation.py', 'new_regions', 'new_regions'),
    'additional-residential': ('residential_allocation.py', 'residential_reference/residential_allocation_inputs.npz', 'residential_allocation'),
    'additional-gnss': ('independent_gnss_transfer.py', '.', 'gnss_transfer'),
}
QUICK = ['partial-identification', 'published-method-benchmark',
         'reporting-tolerance', 'acs-population-benchmark']
STANDARD = [*QUICK, 'dwr-allocation', 'population-analysis', 'sampling-analysis', 'sampling-spatial',
            'paired', 'direct-geometry', 'observed-los-masking',
            'synthetic-masking', 'dwr-masking', 'allocation', 'cross-sensitivity', 'followup', 'spatial', 'population-benchmark', *ADDITIONAL]


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def check_code():
    manifest = ROOT / 'MANIFEST.sha256.json'
    if not manifest.is_file():
        raise FileNotFoundError('Missing MANIFEST.sha256.json; use the published release checkout.')
    records = json.loads(manifest.read_text(encoding='utf-8'))
    for record in records:
        path = ROOT / record['path']
        if not path.is_file() or sha256(path) != record['sha256']:
            raise RuntimeError('Missing or changed release file: ' + record['path'])
    sources = sorted(ROOT.rglob('*.py'))
    for path in sources:
        if '__pycache__' not in path.parts and 'reproduction_output' not in path.parts:
            ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path))
    report = {'manifest_hashes_verified': len(records), 'python_sources_parsed': len(sources)}
    print(json.dumps(report), flush=True)
    return report


def prepare(workspace):
    workspace.mkdir(parents=True, exist_ok=True)
    for source_root, destination in [(ROOT / 'data/strengthening', workspace / 'strengthening'),
                                     (ROOT / 'data/workspace', workspace / 'analysis')]:
        for source in source_root.rglob('*'):
            if not source.is_file():
                continue
            target = destination / source.relative_to(source_root)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                shutil.copy2(source, target)
    (workspace / 'analysis').mkdir(exist_ok=True)
    matched=workspace/'analysis/inputs/bangkok_matched'
    matched.mkdir(parents=True,exist_ok=True)
    for name in ['matched_quality.npz','geometry.json']:
        target=matched/name
        if not target.exists():shutil.copy2(workspace/'strengthening/results/timeseries'/name,target)
    return workspace / 'strengthening', workspace / 'analysis'


def run_mode(mode, workspace, strengthening_override=None):
    strengthening, analysis = prepare(workspace)
    if strengthening_override:
        strengthening = strengthening_override.resolve()
    env = os.environ.copy()
    env.update({'PYTHONUTF8': '1', 'MPLBACKEND': 'Agg',
                'INSAR_CORE_ROOT': str(strengthening),
                'INSAR_STRENGTHENING_ROOT': str(strengthening),
                'INSAR_PROJECT': str(strengthening), 'INSAR_PROJECT_ROOT': str(strengthening.parent),
                'INSAR_FOLLOWUP_ROOT': str(analysis), 'INSAR_PAIRED_ROOT': str(analysis),
                'INSAR_BASE': str(analysis), 'BASE': str(analysis),
                'PROJROOT': str(strengthening.parent)})
    if mode in ADDITIONAL:
        name, input_name, folder = ADDITIONAL[mode]
        script = ROOT / 'additional' / name
        output = workspace / 'additional' / folder
        env.update({'INSAR_BASELINE_REVISION': str(ROOT), 'INSAR_NEW_OUTPUT': str(output)})
        extra = [] if input_name is None else ['--inputs', str(ROOT/'data/additional'/input_name), '--output', str(output)]
        if mode == 'additional-gnss': extra += ['--sample-source', 'windows']
    elif mode in CORE:
        name, extra = CORE[mode]
        script = ROOT / 'core' / name
    elif mode in FOLLOWUP:
        script = ROOT / 'followup' / FOLLOWUP[mode]
        extra = []
    else:
        script = ROOT / 'paired_support' / PAIRED[mode]
        extra = []
    start = time.perf_counter()
    print('REPRODUCE ' + mode, flush=True)
    subprocess.run([sys.executable, '-X', 'utf8', str(script), *extra], env=env,
                   cwd=analysis, check=True)
    return {'mode': mode, 'seconds': round(time.perf_counter() - start, 3)}


def compare(workspace, modes=None):
    from verify_results import verify_results
    return verify_results(ROOT, workspace, modes)


def main():
    choices = ['check-code', 'prepare', 'quick', 'standard', 'additional', 'compare', 'fetch-inputs',
               *CORE, *FOLLOWUP, *PAIRED, *ADDITIONAL]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=choices)
    parser.add_argument('--workspace', type=Path, default=ROOT / 'reproduction_output')
    parser.add_argument('--strengthening-root', type=Path)
    parser.add_argument('--dataset', choices=['acs', 's20', 'census-blocks', 'dwr-polygons', 'worldpop-thailand'])
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    if args.mode == 'check-code':
        check_code()
        return
    if args.mode == 'fetch-inputs':
        from fetch_reviewer_inputs import fetch_inputs
        prepare(workspace)
        if args.dataset is None:
            parser.error('--dataset is required for fetch-inputs')
        fetch_inputs(ROOT, workspace, args.dataset)
        return
    if args.mode == 'compare':
        compare(workspace)
        return
    check_code()
    if args.mode == 'prepare':
        strengthening, analysis = prepare(workspace)
        print(json.dumps({'strengthening': str(strengthening), 'analysis': str(analysis)}))
        return
    modes = QUICK if args.mode == 'quick' else STANDARD if args.mode == 'standard' else list(ADDITIONAL) if args.mode == 'additional' else [args.mode]
    timings = [run_mode(mode, workspace, args.strengthening_root) for mode in modes]
    if args.mode in {'annual-median','persistence'}:
        report={'status':'GENERATED','mode':args.mode,'timings':timings,
                'scope':'Optional upstream-processing step. No frozen output comparator is configured for this mode; inspect its documented method and outputs.'}
        print(json.dumps(report),flush=True)
        return
    result = compare(workspace, modes)
    result['timings'] = timings
    (workspace / 'reproduction_verification.json').write_text(
        json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'modes': modes, 'comparisons': len(result['tables'])}), flush=True)


if __name__ == '__main__':
    main()
