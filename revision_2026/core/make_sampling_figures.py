"""Actual acquisition graphs, baseline composition and descriptive spatial stability."""
from datetime import datetime
import json
import numpy as np,pandas as pd
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection
from matplotlib.patches import Arc
from make_figures import plt,save
from common import ROOT,OUT

def main():
    old=pd.read_csv(OUT/'fine/pair_exposure.csv')
    legacy=sorted(old.pair.unique())
    expanded=sorted(pd.read_csv(OUT/'sampling/expanded_pair_exposure.csv').pair.unique())
    hold=json.loads((OUT/'timeseries/holdout_design.json').read_text(encoding='utf-8'))
    full=sorted(hold['training']+hold['held_out'])
    fig=plt.figure(figsize=(9,6.5),constrained_layout=True);gs=fig.add_gridspec(3,2)
    ax=fig.add_subplot(gs[:2,:])
    for j,(name,pairs,color) in enumerate([('Legacy 28',legacy,'#0072b2'),('Expanded 224',expanded,'#009e73'),('Matched 437',full,'#d55e00')]):
        level=2-j;dates=set()
        for pair in pairs:
            da,db=[mdates.date2num(datetime.strptime(v,'%Y%m%d')) for v in pair.split('_')];dates.update([da,db]);span=db-da
            height=.08+.56*min(span/276,1)
            ax.add_patch(Arc(((da+db)/2,level),width=span,height=height,theta1=0,theta2=180,lw=.5,color=color,alpha=.45))
        ax.scatter(sorted(dates),np.full(len(dates),level),s=4,color=color)
    ax.set_yticks([2,1,0],['Legacy 28','Expanded 224','Matched 437']);ax.set_ylim(-.12,2.65);ax.set_xlim(mdates.date2num(datetime(2015,10,1)),mdates.date2num(datetime(2023,2,1)));ax.xaxis.set_major_locator(mdates.YearLocator());ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'));ax.set_title('a  Acquisition-date graphs (arcs are pair edges)',loc='left')
    bx=fig.add_subplot(gs[2,0]);bins=[0,12,24,48,96,192,300]
    for pairs,label,color in [(legacy,'Legacy 28','#0072b2'),(expanded,'Expanded 224','#009e73')]:
        days=np.array([(datetime.strptime(p[9:],'%Y%m%d')-datetime.strptime(p[:8],'%Y%m%d')).days for p in pairs]);counts=np.histogram(days,bins=bins)[0]/len(days)*100
        bx.plot(np.arange(6),counts,'o-',label=label,color=color)
    bx.set_xticks(range(6),['0–11','12–23','24–47','48–95','96–191','192–299'],rotation=25,ha='right');bx.set(title='b  Temporal baseline composition',xlabel='Baseline bin (days)',ylabel='Pairs (%)');bx.legend()
    cx=fig.add_subplot(gs[2,1]);rank=pd.read_csv(OUT/'sampling/spatial_rank_stability.csv')
    for metric,label in [('deficit','Deficit amount'),('deficit_share','Deficit share')]:
        z=rank[(rank.metric==metric)&rank['sample'].str.startswith('nested')].copy();z['n']=z['sample'].str.replace('nested','').astype(int);z=z.sort_values('n');cx.plot(z.n,z.spearman,'o-',label=label)
    cx.set(title='c  Block ranks versus the 224-pair pool',xlabel='Pairs in nested sample',ylabel='Spearman correlation',ylim=(.99,1.0005));cx.set_xticks([28,56,112,224]);cx.legend(loc='lower right');save(fig,'S2_network_ranks')

if __name__=='__main__':main()
