"""Compile the local manuscript/SI/response; final mode rejects draft markers."""
import argparse,re,subprocess
from common import ROOT,dump,sha256

PAPER=ROOT/'paper'

def escape(s):
    table={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
    return ''.join(table.get(c,c) for c in s)

def inline(s):
    pieces=re.split(r'(\*\*.*?\*\*|`[^`]+`)',s)
    return ''.join(r'\textbf{'+escape(p[2:-2])+'}' if p.startswith('**') else r'\path{'+p[1:-1]+'}' if p.startswith('`') else escape(p) for p in pieces)

def response_tex():
    lines=(PAPER/'response_to_reviewers.md').read_text(encoding='utf-8').splitlines();body=[]
    for line in lines:
        if line.startswith('# '):body.append(r'\begin{center}\LARGE\bfseries '+escape(line[2:])+r'\end{center}\vspace{1em}')
        elif line.startswith('### '):body.append(r'\subsection*{'+escape(line[4:])+'}')
        elif line.startswith('## '):body.append(r'\section*{'+escape(line[3:])+'}')
        else:body.append(inline(line))
    pre=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=23mm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern,xurl,hyperref}
\hypersetup{hidelinks}
\DeclareUnicodeCharacter{00B2}{\ensuremath{^{2}}}
\emergencystretch=2em
\begin{document}
'''
    (PAPER/'response.tex').write_text(pre+'\n'.join(body)+'\n'+r'\end{document}'+'\n',encoding='utf-8')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--draft',action='store_true');args=ap.parse_args()
    if not args.draft:
        for name in ['matched_methods.tex','matched_results.tex','matched_supplement.tex','response_to_reviewers.md']:
            s=(PAPER/name).read_text(encoding='utf-8').lower()
            if any(x in s for x in ['working-draft marker','outcome pending','still being processed','remains explicitly pending']):raise ValueError('Unfinished scientific content: '+name)
    response_tex();records=[]
    for name in ['manuscript','supplementary','response']:
        commands=[['pdflatex','-interaction=nonstopmode','-halt-on-error',name+'.tex']]
        if name!='response':commands.append(['bibtex',name])
        commands += [['pdflatex','-interaction=nonstopmode','-halt-on-error',name+'.tex']]*2
        for n,command in enumerate(commands):
            result=subprocess.run(command,cwd=PAPER,capture_output=True,text=True,encoding='utf-8',errors='replace')
            (PAPER/(name+f'_build{n}.txt')).write_text(result.stdout+result.stderr,encoding='utf-8')
            if result.returncode:raise RuntimeError('Build failed: '+str(command))
        log=(PAPER/(name+'.log')).read_text(encoding='utf-8',errors='replace')
        problems=[line for line in log.splitlines() if any(k in line for k in ['Overfull','undefined','multiply defined'])]
        records.append(dict(document=name,draft=args.draft,pdf_sha256=sha256(PAPER/(name+'.pdf')),log_findings=problems))
        print(name,'compiled;',len(problems),'log findings',flush=True)
    dump(PAPER/'build_audit.json',records)

if __name__=='__main__':main()
