"""Publication figures drawn exclusively from saved experiment outputs."""
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullLocator,MaxNLocator
from common import OUT,ROOT,dump,sha256

FIG=ROOT/'paper/figures';FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':11,'axes.labelsize':11,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':220,'axes.spines.top':False,'axes.spines.right':False})
COLORS={'area_uniform':'#d55e00','cell_center':'#0072b2','building_weighted':'#009e73'}
EXPOSURE_COLORS={'area_m2':'#0072b2','population':'#d55e00','builtup_m2':'#009e73'}
LABELS={'area_uniform':'Uniform area','cell_center':'Cell center','building_weighted':'Built-area weights'}

def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');fig.savefig(FIG/(name+'.png'),bbox_inches='tight');plt.close(fig)
    print('FIGURE',name,flush=True)

def map_panel(ax,a,x,y,title,cmap='viridis',norm=None):
    im=ax.imshow(a,origin='lower',extent=[x[0],x[-1],y[0],y[-1]],cmap=cmap,norm=norm,interpolation='nearest',aspect=1/np.cos(np.deg2rad(np.mean(y))))
    ax.set_title(title,loc='left');ax.set_xlabel('Longitude (degrees)');ax.set_ylabel('Latitude (degrees)')
    ax.xaxis.set_major_locator(MaxNLocator(3))
    return im

def main():
    g=np.load(OUT/'fine/grid.npz');s=np.load(OUT/'fine/mean_support.npz');x,y=g['x_edges'],g['y_edges'];d=g['domain']
    fig,axs=plt.subplots(1,3,figsize=(9,4.5),constrained_layout=True)
    im=map_panel(axs[0],np.where(d,g['population'],np.nan),x,y,'a  GHSL 2020\npopulation',norm=LogNorm(1,10000));fig.colorbar(im,ax=axs[0],shrink=.35,label='Population per native pixel')
    im=map_panel(axs[1],np.where(d,s['q30'],np.nan),x,y,'b  Mean pair support\n(q = 0.3)');im.set_clim(0,1);fig.colorbar(im,ax=axs[1],shrink=.35,label='Area-pair fraction')
    im=map_panel(axs[2],np.where(d,g['population']*(1-s['q30']),np.nan),x,y,'c  Population-weighted\nsupport deficit',norm=LogNorm(.01,1000));fig.colorbar(im,ax=axs[2],shrink=.35,label='Mean allocation units per pixel')
    save(fig,'F1_chao_spatial')

    curves=pd.read_csv(OUT/'fine/support_curves.csv');fig,axs=plt.subplots(1,2,figsize=(9,3.5),constrained_layout=True)
    for e,label in [('area_m2','Area'),('builtup_m2','Built surface'),('population','Population')]:
        z=curves[(curves.domain=='geographic')&(curves.support=='q30')&(curves.exposure==e)]
        axs[0].plot(z.threshold,z.share_below,label=label,lw=2,color=EXPOSURE_COLORS[e])
    axs[0].set(title='a  Weighted distribution of fine-pixel support',xlabel='Mean area-pair support threshold',ylabel='Share below threshold',xlim=(0,1),ylim=(0,1));axs[0].legend()
    base=pd.read_csv(OUT/'fine/allocation_summary.csv')
    for e,label in [('population','Population'),('builtup_m2','Built surface'),('area_m2','Area')]:
        z=base[(base.domain=='geographic')&(base.exposure==e)&base.support.isin(['q20','q30','q40'])].sort_values('support')
        axs[1].plot([.2,.3,.4],1-z.weighted_support.values,'o-',label=label,color=EXPOSURE_COLORS[e])
    axs[1].set(title='b  Sensitivity to the coherence cutoff',xlabel='Coherence cutoff',ylabel='Mean weighted support deficit',ylim=(0,1));axs[1].legend();save(fig,'F2_support_curves')

    scale=pd.read_csv(OUT/'allocation/scale_origin_comparison.csv');fig,axs=plt.subplots(1,3,figsize=(9,3.4),constrained_layout=True)
    for ax,e,title in [(axs[0],'population','a  Population allocation'),(axs[1],'builtup_m2','b  Built-surface allocation')]:
        for method in ['area_uniform','cell_center']:
            z=scale[(scale.domain=='geographic')&(scale.support=='q30')&(scale.exposure==e)&(scale.method==method)]
            q=z.groupby('cell_size_m').error_fraction_total.agg(['min','max','mean'])*100
            ax.plot(q.index,q['mean'],'o-',color=COLORS[method],label=LABELS[method]);ax.fill_between(q.index,q['min'],q['max'],alpha=.2,color=COLORS[method])
        ax.axhline(0,color='.4',lw=.8);ax.set(title=title,xlabel='Cell size (m)',ylabel='Signed error (% of AOI total)',xscale='log');ax.set_xticks([250,500,1000,2000],['250','500','1000','2000']);ax.xaxis.set_minor_locator(NullLocator());ax.legend()
    cells=pd.read_csv(OUT/'allocation/primary_1000m_cells.csv');z=cells[cells.population>0]
    hb=axs[2].hexbin(z.support,(z.uniform-z.reference)/z.population,gridsize=40,mincnt=1,bins='log',cmap='cividis')
    axs[2].axhline(0,color='.3',lw=.8);axs[2].set(title='c  Local allocation error (1 km)',xlabel='Cell mean area-pair support',ylabel='Error / cell population');fig.colorbar(hb,ax=axs[2],label='Cells per hexagon')
    save(fig,'F3_allocation_scale')

    masks=pd.read_csv(OUT/'masking/replicates.csv');masks['error_pct_total']=100*masks.signed_error/masks.total_population;fig,axs=plt.subplots(1,3,figsize=(9,3.6),sharey=True,constrained_layout=True)
    mechanisms=['random','spatial','density_positive','density_negative'];names=['Random','Spatial','Dense\nmissing','Sparse\nmissing']
    for ax,fraction in zip(axs,[.1,.3,.5]):
        for j,method in enumerate(COLORS):
            z=masks[(masks.cell_size_fine_pixels==10)&(masks.origin_fine_pixels==0)&(masks.missing_pixel_fraction==fraction)&(masks.method==method)].groupby('mechanism').error_pct_total.agg(['mean',lambda v:v.quantile(.025),lambda v:v.quantile(.975)]).loc[mechanisms]
            yy=z.iloc[:,0].values;low=z.iloc[:,1].values;high=z.iloc[:,2].values
            ax.errorbar(np.arange(4)+(j-1)*.19,yy,yerr=np.maximum(0,np.vstack([yy-low,high-yy])),fmt='o',capsize=2,color=COLORS[method],label=LABELS[method])
        ax.set_xticks(range(4),names);ax.axhline(0,color='.4',lw=.8);ax.set_title(f'{int(fraction*100)}% pixels removed');ax.tick_params(axis='x',labelsize=10)
    axs[0].set_ylabel('Signed error (% of test-domain population)');axs[0].legend(loc='upper left');save(fig,'F4_masking')

    model=pd.read_csv(OUT/'population/model_epoch_comparison.csv');land=pd.read_csv(OUT/'population/landcover_strata.csv');fig,axs=plt.subplots(1,3,figsize=(9,3.5),constrained_layout=True)
    for name,label in [('GHSL_2020_3ss','GHSL 2020'),('WorldPop_2020_3ss','WorldPop 2020')]:
        z=model[(model.domain=='geographic')&(model.comparison=='joint_valid')&(model.model==name)&model.support.isin(['q20','q30','q40'])].sort_values('support')
        axs[0].plot([.2,.3,.4],z.unsupported_share*100,'o-',label=label)
    axs[0].set(title='a  Alternative population models',xlabel='Coherence cutoff',ylabel='Mean deficit (% of model total)');axs[0].legend()
    z=model[(model.domain=='geographic')&(model.comparison=='joint_valid')&(model.support=='q30')&model.model.isin(['GHSL_2015_30ss','GHSL_2020_30ss'])].sort_values('model')
    axs[1].bar(['2015','2020'],z.unsupported_share*100,color=['#7570b3','#1b9e77']);axs[1].set(title='b  GHSL epochs at fixed 30 arcsec',ylabel='Mean deficit (% of model total)',ylim=(0,30))
    classes=[50,10,40,80];names=['Built-up','Tree cover','Cropland','Water']
    for j,(e,label) in enumerate([('area_m2','Area'),('population','Population')]):
        z=land[(land.support=='q30')&(land.exposure==e)].set_index('class_code').loc[classes]
        axs[2].bar(np.arange(4)+(j-.5)*.36,(1-z.weighted_support)*100,width=.36,label=label,color=EXPOSURE_COLORS[e])
    axs[2].set_xticks(range(4),names,rotation=25,ha='right');axs[2].set(title='c  Fractional land-cover strata',ylabel='Mean weighted deficit (%)');axs[2].legend();save(fig,'F5_population_models')

    if (OUT/'dwr/grid.npz').exists():
        z=np.load(OUT/'dwr/grid.npz');fig,axs=plt.subplots(1,3,figsize=(9,3.7),constrained_layout=True)
        im=map_panel(axs[0],z['population'],z['x_edges'],z['y_edges'],'a  Fresno population',norm=LogNorm(.1,3000));fig.colorbar(im,ax=axs[0],shrink=.7,label='Population per GHSL pixel')
        im=map_panel(axs[1],z['support'],z['x_edges'],z['y_edges'],'b  Published DWR\nobservation coverage');im.set_clim(0,1);fig.colorbar(im,ax=axs[1],shrink=.7,label='Covered area fraction')
        f=pd.read_csv(OUT/'dwr/scale_origin_comparison.csv')
        for method in ['area_uniform','cell_center']:
            q=f[(f.exposure=='population')&(f.method==method)].copy();q['error_pct']=100*q.signed_error/q.total;q=q.groupby('cell_size_m').error_pct.agg(['min','max','mean'])
            axs[2].plot(q.index,q['mean'],'o-',color=COLORS[method],label=LABELS[method]);axs[2].fill_between(q.index,q['min'],q['max'],alpha=.2,color=COLORS[method])
        axs[2].set(title='c  Same-total allocation comparison',xlabel='Cell size (m)',ylabel='Signed error (% of AOI population)',xscale='log');axs[2].set_xticks([250,500,1000,2000],['250','500','1000','2000']);axs[2].xaxis.set_minor_locator(NullLocator());axs[2].legend();save(fig,'F6_dwr')

    if (OUT/'sampling/nested_summary.csv').exists():
        z=pd.read_csv(OUT/'sampling/nested_summary.csv');reps=pd.read_csv(OUT/'sampling/conditional_sampling_replicates.csv');rob=pd.read_csv(OUT/'sampling/temporal_robustness.csv');fig,axs=plt.subplots(1,3,figsize=(9,3.5),constrained_layout=True)
        q=reps.groupby('pairs').unsupported.agg(['mean',lambda v:v.quantile(.025),lambda v:v.quantile(.975)]);axs[0].plot(q.index,q.iloc[:,0]/1e6,'o-');axs[0].fill_between(q.index,q.iloc[:,1]/1e6,q.iloc[:,2]/1e6,alpha=.2);axs[0].set(title='a  Conditional sample variation',xlabel='Pairs sampled from the fixed 224-pair pool',ylabel='Mean deficit (million allocation units)')
        for sample in ['legacy28','expanded224']:
            q=rob[(rob['sample']==sample)&(rob.support=='q30')&(rob.exposure=='population')&rob.subset.str.startswith('leave_')]
            axs[1].plot(q.subset.str[6:].astype(int),q.mean_unsupported/1e6,'o-',label=sample)
        axs[1].set(title='b  Leave-one-year-out sensitivity',xlabel='Year omitted',ylabel='Mean deficit (million allocation units)');axs[1].legend();axs[1].tick_params(axis='x',rotation=30)
        for n in [28,56,112,224]:
            curves=pd.read_csv(OUT/'sampling/nested_curves.csv');q=curves[(curves.pairs==n)&(curves.support=='q30')&(curves.exposure=='population')];axs[2].plot(q.threshold,q.below_share,label=str(n))
        axs[2].set(title='c  Predeclared nested samples',xlabel='Mean support threshold',ylabel='Population share below threshold');axs[2].legend(title='Pairs');save(fig,'S1_sampling')
    dump(FIG/'figure_manifest.json',[dict(file=p.name,sha256=sha256(p)) for p in sorted(FIG.glob('*.pdf'))])

if __name__=='__main__':main()
