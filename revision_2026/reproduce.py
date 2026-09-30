"""Run retained scientific analyses with explicit input and output roots."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
FOLLOWUP = {
    'followup': 'followup_quality.py',
    'census-support': 'fresno_census_block_analysis.py',
    'annual-median': 'build_dwr_annual_rate_median.py',
    'census-rates': 'dwr_census_velocity_analysis.py',
    'temporal': 'fresno_temporal_class_stability.py',
    'persistence': 'fresno_six_year_persistence.py',
    'spatial': 'bangkok_spatial_effect_decomposition.py',
    'population-benchmark': 'fresno_census_population_model_benchmark.py',
    'download-population': 'download_fresno_population_rasters.py',
    'partial-identification': 'partial_identification.py',
    'published-method-benchmark': 'published_method_benchmark.py',
    'reporting-tolerance': 'reporting_tolerance_audit.py',
    'acs-population-benchmark': 'acs_population_benchmark.py',
}
PAIRED = {'paired': 'paired_support_pilot.py', 'direct-geometry': 'direct_geometry_audit.py',
          'observed-los-masking': 'actual_los_masking.py'}

def check_code():
    files = sorted(ROOT.rglob('*.py'))
    for path in files:
        ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path))
    manifest = ROOT / 'MANIFEST.sha256.json'
    checked = 0
    if manifest.is_file():
        for record in json.loads(manifest.read_text(encoding='utf-8')):
            path = ROOT / record['path']
            if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
                raise RuntimeError(f'Hash mismatch: {record["path"]}')
            checked += 1
    print(json.dumps({'python_syntax_files': len(files), 'manifest_hashes_verified': checked}))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['check-code', *FOLLOWUP, *PAIRED])
    parser.add_argument('--workspace', type=Path, help='External follow-up input/output directory.')
    parser.add_argument('--strengthening-root', type=Path,
                        help='Retained strengthening input root; read-only inputs are supported.')
    args = parser.parse_args()
    if args.mode == 'check-code':
        check_code()
        return
    if args.workspace is None:
        parser.error('--workspace is required for analysis modes.')
    workspace = args.workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update({'INSAR_FOLLOWUP_ROOT': str(workspace), 'INSAR_PAIRED_ROOT': str(workspace),
                        'INSAR_BASE': str(workspace), 'BASE': str(workspace)})
    if args.strengthening_root is not None:
        project = args.strengthening_root.resolve()
        if not project.is_dir():
            parser.error('The strengthening input root does not exist.')
        environment.update({'INSAR_STRENGTHENING_ROOT': str(project), 'INSAR_PROJECT': str(project),
                            'PROJROOT': str(project.parent)})
    required_project_modes = {'followup', 'census-support', 'census-rates', 'persistence', 'partial-identification', 'published-method-benchmark', *PAIRED}
    if args.mode in required_project_modes and args.strengthening_root is None:
        parser.error('--strengthening-root is required for this analysis.')
    if args.mode in FOLLOWUP:
        script = ROOT / 'followup' / FOLLOWUP[args.mode]
    else:
        script = ROOT / 'paired_support' / PAIRED[args.mode]
    print('Running', script.relative_to(ROOT), 'with workspace', workspace, flush=True)
    subprocess.run([sys.executable, '-X', 'utf8', str(script)], env=environment, cwd=workspace, check=True)

if __name__ == '__main__':
    from reviewer_reproduce import main as reviewer_main
    reviewer_main()
