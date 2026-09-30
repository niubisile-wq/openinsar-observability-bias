"""Verify manifest content hashes and relocation without downloading inputs."""
from pathlib import Path
import os,sys,json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;PROJECT=ROOT.parent;MASTER=PROJECT.parent.parent
os.environ.setdefault('INSAR_PROJECT_ROOT',str(PROJECT))
if 'INSAR_BACKUP_ROOT' not in os.environ:
    matches=[p.parent for p in (MASTER/'02_历史备份').glob('*/nature/revision_scientific_reports_v1') if (p/'inputs_needed/licsar_chao_28_unw').exists()]
    if len(matches)!=1:raise RuntimeError('Set INSAR_BACKUP_ROOT to the retained nature input root.')
    os.environ['INSAR_BACKUP_ROOT']=str(matches[0])
sys.path.insert(0,str(ROOT))
from common import resolve_recorded_path

def main():
    targets=[]
    for name in ['sampling_manifest.json','timeseries_manifest.json']:
        for rec in json.loads((ROOT/'external/licsar_expanded'/name).read_text(encoding='utf-8')):
            assert rec['status']=='ok'
            for kind,v in rec['files'].items():targets.append((v['path'],v['sha256'],name,kind))
    for rec in json.loads((ROOT/'results/fine/pair_manifest.json').read_text(encoding='utf-8')):
        for kind in ['cc','unw']:targets.append((rec[kind],rec[kind+'_sha256'],'fine/pair_manifest.json',kind))
    checked={};rows=[]
    for original,expected,manifest,kind in targets:
        p=resolve_recorded_path(original)
        if p not in checked:
            h=hashlib.sha256()
            with p.open('rb') as f:
                for b in iter(lambda:f.read(2**20),b''):h.update(b)
            checked[p]=h.hexdigest()
        assert checked[p]==expected,(str(p),manifest)
        rows.append(dict(recorded=original,resolved=str(p),sha256=expected,manifest=manifest,kind=kind))
    out=HERE/'evidence';out.mkdir(exist_ok=True)
    (out/'relocated_input_hashes.json').write_text(json.dumps(dict(unique_files=len(checked),manifest_references=len(rows),
        all_hashes_match=True,inputs=rows),indent=2,ensure_ascii=False),encoding='utf-8')
    print('Verified',len(checked),'unique original rasters;',len(rows),'manifest references; all recorded hashes match.')

if __name__=='__main__':main()
