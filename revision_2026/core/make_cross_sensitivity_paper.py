"""Generate cross-experiment figure, numeric table and prose from retained outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent;D=R/'results/cross_sensitivity';P=R/'paper'

def main():
    audit=json.loads((D/'audit.json').read_text());df=pd.read_csv(D/'cross_sensitivity.csv')
    z={(r['model'],r['sample']):r for r in audit['summary']}
    gl=z['GHSL_2020_3ss','legacy28'];ge=z['GHSL_2020_3ss','expanded224']
    wl=z['WorldPop_2020_3ss','legacy28'];we=z['WorldPop_2020_3ss','expanded224']
    def span(r,a,b,d=2):return f'{r[a]:.{d}f}--{r[b]:.{d}f}'
    macros={'CrossGHLow':min(gl['pct_min'],ge['pct_min']),'CrossGHHigh':max(gl['pct_max'],ge['pct_max']),
            'CrossWPLow':min(wl['pct_min'],we['pct_min']),'CrossWPHigh':max(wl['pct_max'],we['pct_max'])}
    (P/'cross_numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+f'{v:.1f}'+'}' for k,v in macros.items())+'\n',encoding='utf-8')
    results=(r'\paragraph{Crossing population model and pair selection.} '+
        f'On the fixed GHSL--WorldPop common domain, both models are normalized to {audit["common_population"]/1e6:.3f} million population units. '+
        f'At 1 km, uniform allocation exceeds the native-model deficit by {span(gl,"pct_min","pct_max")}\\% for GHSL with the legacy 28 pairs and {span(ge,"pct_min","pct_max")}\\% with the expanded 224 pairs. '+
        f'The corresponding WorldPop differences are {span(wl,"pct_min","pct_max")}\\% and {span(we,"pct_min","pct_max")}\\%. '+
        f'All 16 model--sample--origin combinations have positive differences; normalized by the common total, these span {min(gl["pp_min"],ge["pp_min"]):.2f}--{max(gl["pp_max"],ge["pp_max"]):.2f} percentage points for GHSL and {min(wl["pp_min"],we["pp_min"]):.2f}--{max(wl["pp_max"],we["pp_max"]):.2f} for WorldPop. '+
        'The direction therefore persists across these inputs, while its magnitude is strongly population-model dependent. The original approximately 58\\% difference is specific to its full-domain GHSL case. Changing 28 to 224 pairs changes selection and composition as well as sample size; the comparison is not an isolated sample-size effect. Supplementary S12 reports every origin and the numerical convergence check.\n')
    (P/'cross_results.tex').write_text(results,encoding='utf-8')
    table=[r'\begin{table}[htbp]\centering\footnotesize',
      r'\caption{Crossed population-model and pair-selection comparison at 1 km and $q=0.3$. Both models use the same domain and total. Deficits and differences are in millions of allocation units (M). Ranges span four grid origins; they are not confidence intervals. $\Delta=U_{\rm uniform}-U_{\rm native}$, and pp denotes percentage points of the common total.}\label{tab:cross}',
      r'\setlength{\tabcolsep}{3.5pt}',r'\begin{tabular}{llrrrrr}\toprule',
      r'Model & Pairs & Native (M) & Uniform (M) & $\Delta$ (M) & $\Delta$/native (\%) & $\Delta$/total (pp)\\\midrule']
    order=[('GHSL_2020_3ss','legacy28'),('GHSL_2020_3ss','expanded224'),('WorldPop_2020_3ss','legacy28'),('WorldPop_2020_3ss','expanded224')]
    for key in order:
        r=z[key];model='GHSL' if key[0].startswith('GHSL') else 'WorldPop';sample='28' if key[1]=='legacy28' else '224'
        table.append(f'{model} & {sample} & {r["native_deficit"]/1e6:.3f} & {r["uniform_min"]/1e6:.3f}--{r["uniform_max"]/1e6:.3f} & {r["delta_min"]/1e6:.3f}--{r["delta_max"]/1e6:.3f} & '+span(r,'pct_min','pct_max')+' & '+span(r,'pp_min','pp_max')+r'\\')
    table += [r'\bottomrule\end{tabular}\end{table}']
    pd.DataFrame(audit['summary']).to_csv(D/'summary.csv',index=False)
    cons=pd.read_csv(D/'conservation.csv');drift=cons[cons.resolution_m.eq(25)].relative_drift.abs().max()
    supp=(r'\section*{S12. Crossed population-model and pair-selection sensitivity}'+'\n'+
      f'The common domain contains {audit["common_cells"]:,} native GHSL pixels and {audit["common_population"]:,.3f} GHSL population units; {audit["excluded_ghsl_population"]:,.3f} GHSL units from the full primary footprint are excluded. '+
      'It intersects the primary footprint with full coverage of GHSL 2020 and WorldPop 2020 only; the four-layer epoch comparison in S6 uses a separately defined intersection. Each model is normalized to the GHSL total in this two-model domain before any support calculation. We cross these two allocations with the original 28-pair and expanded 224-pair mean joint support fields, at $q=0.3$. Neither population product is treated as demographic truth.\n\n'+
      'All comparisons use 1-km cells in UTM zone 47N and offsets $(0,0)$, $(500,0)$, $(0,500)$ and $(500,500)$ m. Within each combination, native-model and uniform allocations have identical cell population totals and the same included-domain area. We project extensive area, supported area, population and supported population separately onto a 50-m integration grid. The population--support product is formed before projection. Retaining partial cells avoids full-cell extrapolation outside the domain.\n\n'+
      '\n'.join(table)+'\n\n'+
      f'Production first and cross moments conserve source totals to relative tolerance $10^{{-8}}$, and aggregation over each origin preserves those projected totals. Direct native-grid integration agrees with the coarse cross-moment reference. The 25-m diagnostic retains raw reprojection drift without renormalization (maximum absolute relative moment drift {drift:.3g}). Its maximum change in uniform-minus-native difference is {audit["max_integration_difference_pp_total"]:.6f} percentage points of the common total, below the prespecified 0.01-percentage-point numerical tolerance. These are numerical sensitivity checks, not uncertainty intervals.\n\n'+
      'All 16 production rows have positive uniform-minus-native differences. GHSL gives a larger effect than WorldPop under both pair selections. The result supports a persistent direction for these selected inputs, while showing that a single percentage cannot represent population-model sensitivity. Because pair number and pair composition change together, no causal effect of sample size is inferred. The four origins are deterministic sensitivity cases, not independent replicates.\n\n'+
      r'\begin{figure}[htbp]\centering\includegraphics[width=\linewidth]{figures/S5_cross_sensitivity.pdf}'+'\n'+
      r'\caption{Crossed aggregation sensitivity. Dots show four grid origins and horizontal lines span them. Left: difference relative to the native-model deficit. Right: the same difference as percentage points of the common population total. Colors identify population model; the 28 and 224 labels identify distinct pair selections.}\label{fig:cross}'+'\n'+r'\end{figure}'+'\n\n'+
      r'Complete 16-row results, the 25-m comparison, conservation ledger and input SHA-256 values are retained in \path{results/cross_sensitivity/}. The analysis and figure can be regenerated with \path{experiment_cross_sensitivity.py} and \path{make_cross_sensitivity_paper.py}.'+'\n')
    (P/'cross_supplement.tex').write_text(supp,encoding='utf-8')
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(7,2.75),sharey=True,layout='constrained')
    for ax,column,title in zip(axes,['difference_pct_native','difference_pp_total'],['a  Relative to native deficit','b  Relative to common total']):
        for i,key in enumerate(order):
            vals=df[df.model.eq(key[0])&df['sample'].eq(key[1])][column].to_numpy()
            color='#276A9D' if key[0].startswith('GHSL') else '#C36A26'
            ax.plot([vals.min(),vals.max()],[i,i],color=color,lw=2)
            ax.scatter(vals,i+np.array([-.075,-.025,.025,.075]),color=color,s=15,zorder=3)
        ax.set_title(title,loc='left',fontsize=10);ax.axvline(0,color='.5',lw=.6);ax.grid(axis='x',alpha=.2)
        ax.set_xlabel('Difference (%)' if column.endswith('native') else 'Difference (percentage points)')
        ax.set_xlim(left=0);ax.set_ylim(3.45,-.45)
    axes[0].set_yticks(range(4),['GHSL / 28','GHSL / 224','WorldPop / 28','WorldPop / 224'])
    fig.savefig(P/'figures/S5_cross_sensitivity.pdf');fig.savefig(P/'figures/S5_cross_sensitivity.png',dpi=220);plt.close(fig)
    print('Cross-experiment paper files generated',flush=True)

if __name__=='__main__':main()
