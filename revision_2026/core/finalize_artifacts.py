"""Final local artifact integrity and completeness checks; does not publish anything."""
import ast,json,re
from pathlib import Path
from common import ROOT,OUT,dump,sha256

def main():
    paper=ROOT/'paper';build=json.loads((paper/'build_audit.json').read_text());render=json.loads((paper/'qa/render_inventory.json').read_text());visual=json.loads((paper/'qa/visual_review.json').read_text())
    assert len(build)==3 and all(not r['draft'] and not r['log_findings'] for r in build)
    for r in build:
        actual=sha256(paper/(r['document']+'.pdf'));assert actual==r['pdf_sha256']
        q=next(q for q in render if q['document']==r['document']);assert actual==q['pdf_sha256'] and not q['text_near_page_edge']
        assert visual['pdf_sha256'][r['document']]==actual
    assert visual['status']=='reviewed'
    for p in ROOT.glob('*.py'):ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    for p in [paper/'manuscript.tex',paper/'supplementary.tex',paper/'matched_methods.tex',paper/'matched_results.tex',paper/'matched_supplement.tex',paper/'response_to_reviewers.md']:
        t=p.read_text(encoding='utf-8').lower();assert not any(k in t for k in ['working-draft marker','outcome pending','still being processed','remains explicitly pending','@@'])
    manuscript=(paper/'manuscript.tex').read_text(encoding='utf-8')
    abstract=re.search(r'\\abstract\{(.*?)\}\n\\keywords',manuscript,re.S).group(1)
    title=re.search(r'\\title\[.*?\]\{(.*?)\}',manuscript).group(1)
    keywords=re.search(r'\\keywords\{(.*?)\}',manuscript).group(1).split(',')
    expanded=manuscript.replace(r'\input{matched_results}',(paper/'matched_results.tex').read_text(encoding='utf-8'))
    figures=expanded.count(r'\begin{figure}');tables=expanded.count(r'\begin{table}')
    assert len(abstract.split())<=200 and len(title.split())<=20 and len(keywords)<=6 and figures+tables<=8
    assert json.loads((OUT/'timeseries/audit.json').read_text())['all_twelve_stages_complete']
    checks=json.loads((OUT/'numerical_checks.json').read_text());assert checks['allocation_rows']==512 and checks['masking_replicate_rows']==21600 and checks['dwr/masking_replicate_rows']==21600
    records=[]
    files=list(ROOT.glob('*.py'))+list(ROOT.glob('*.md'))+list(ROOT.glob('requirements_*.txt'))+list(ROOT.glob('environment_*.json'))+[ROOT/'config.json']
    files+=list(paper.glob('*.tex'))+list(paper.glob('*.bib'))+list(paper.glob('*.cls'))+list(paper.glob('*.bst'))+list(paper.glob('*.pdf'))+list(paper.glob('*.md'))+list(paper.glob('*.json'))
    files+=list((paper/'figures').glob('*.pdf'))+list((paper/'figures').glob('*.png'))
    files+=[ROOT/'external/licsbas_windows_io/portability.diff',paper/'qa/visual_review.json',paper/'qa/render_inventory.json']
    for p in OUT.rglob('*'):
        if p.is_file() and p.suffix in ['.csv','.json'] and 'GEOC_' not in str(p) and 'TS_' not in str(p):files.append(p)
    for p in (ROOT/'external').rglob('*.json'):
        if not any(x.startswith('features_') for x in p.parts) and 'LiCSBAS2' not in p.parts:files.append(p)
    for p in sorted(set(files)):
        records.append(dict(path=str(p.relative_to(ROOT)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha256(p)))
    dump(ROOT/'delivery_manifest.json',dict(scope='Local scientific deliverables and concise provenance records; raw rasters/binaries have their own retained input manifests',files=records))
    dump(ROOT/'final_audit.json',dict(experiments='E0--E6 completed under the stated interpretation limits',geometry_tests=5,allocation_rows=512,masking_rows_per_region=21600,matched_stages=12,abstract_words=len(abstract.split()),title_words=len(title.split()),keywords=len(keywords),main_figures=figures,main_tables=tables,pdfs={r['document']:dict(pages=r['pages'],sha256=r['pdf_sha256']) for r in render},visual_review='completed and hash-matched',script_syntax='all top-level scripts parse',manifest_files=len(records),publication='No upload or external publication performed'))
    print('FINAL LOCAL AUDIT PASS;',len(records),'manifested files',flush=True)

if __name__=='__main__':main()
