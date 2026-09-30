"""Report the second-region masking benchmark with both error denominators."""
import numpy as np,pandas as pd
from make_figures import plt,save,COLORS,LABELS
from common import OUT

def main():
    d=pd.read_csv(OUT/'dwr/masking/replicates.csv');d['whole_population_error']=100*d.signed_error/d.total_population
    mechanisms=['random','spatial','density_positive','density_negative'];labels=['Random','Spatial','Dense\nmissing','Sparse\nmissing']
    fig,axs=plt.subplots(1,3,figsize=(9,3.6),sharey=True,constrained_layout=True)
    rows=[]
    for ax,fraction in zip(axs,[.1,.3,.5]):
        for j,method in enumerate(COLORS):
            v=d[(d.cell_size_fine_pixels==10)&(d.origin_fine_pixels==0)&(d.missing_pixel_fraction==fraction)&(d.method==method)]
            q=v.groupby('mechanism').whole_population_error.agg(['mean',lambda v:v.quantile(.025),lambda v:v.quantile(.975)]).loc[mechanisms]
            yy=q.iloc[:,0].values;lo=q.iloc[:,1].values;hi=q.iloc[:,2].values
            ax.errorbar(np.arange(4)+(j-1)*.19,yy,yerr=np.maximum(0,np.vstack([yy-lo,hi-yy])),fmt='o',capsize=2,color=COLORS[method],label=LABELS[method])
            for mechanism,g in v.groupby('mechanism'):rows.append(dict(fraction=fraction,method=method,mechanism=mechanism,mean_pct_total=g.whole_population_error.mean(),mean_pct_removed=100*g.relative_error.mean()))
        ax.set_xticks(range(4),labels);ax.axhline(0,color='.4',lw=.8);ax.set_title(f'{int(fraction*100)}% observed pixels removed');ax.tick_params(axis='x',labelsize=10)
    axs[0].set_ylabel('Signed error (% of footprint population)');axs[0].legend(loc='upper left');save(fig,'S4_dwr_masking')
    pd.DataFrame(rows).to_csv(OUT/'dwr/masking/figure_summary.csv',index=False)

if __name__=='__main__':main()
