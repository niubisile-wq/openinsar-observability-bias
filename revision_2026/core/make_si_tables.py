import pandas as pd
from common import OUT,ROOT

PAPER=ROOT/'paper'

def table(path,columns,header,lines):
    path.write_text('\\begin{tabular}{'+columns+'}\\toprule\n'+header+r'\\\midrule'+'\n'+'\n'.join(lines)+'\n'+r'\bottomrule\end{tabular}'+'\n',encoding='utf-8')

def main():
    a=pd.read_csv(OUT/'lineage/encoding_registration_validity_factorial.csv');a=a[(a.domain=='inherited_strong')&(a.exposure=='population')&(a.statistic=='center_majority_binary')]
    lines=[]
    for r in a.itertuples():
        enc='Raw byte (incorrect)' if r.encoding=='raw_byte' else 'Byte / 255';reg='Legacy tiepoint' if r.registration=='legacy_tiepoint' else 'GDAL point'
        lines.append(f'{enc} & {reg} & {"Yes" if r.joint_unwrapped_validity else "No"} & {r.low_center_cells:,} & {r.unsupported/1e6:.6f} '+r'\\')
    table(PAPER/'lineage_table.tex','lllrr','Coherence encoding & Registration & Joint validity & Selected cells & Population',lines)
    s=pd.read_csv(OUT/'masking/summary.csv');s=s[(s.cell_size_fine_pixels==10)&(s.origin_fine_pixels==0)&(s.missing_pixel_fraction==.3)]
    lines=[f'{r.mechanism.replace("_"," ")} & {r.method.replace("_"," ")} & {100*r.mean_relative_error:.2f} & [{100*r.q025_relative_error:.2f}, {100*r.q975_relative_error:.2f}] '+r'\\' for r in s.itertuples()]
    table(PAPER/'masking_table.tex','llrr','Mechanism & Method & Mean error (\%) & 2.5th--97.5th (\%)',lines)

if __name__=='__main__':main()
