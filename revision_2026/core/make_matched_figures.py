"""Same-chain support/quality association and held-out prediction figures."""
import numpy as np,pandas as pd
from matplotlib.colors import ListedColormap,BoundaryNorm,TwoSlopeNorm
from matplotlib.patches import Patch
from make_figures import plt,save,map_panel
from common import OUT

def main():
    z=np.load(OUT/'timeseries/matched_quality.npz');x,y=z['x_edges'],z['y_edges']
    bins=pd.read_csv(OUT/'timeseries/support_quality_bins.csv');held=pd.read_csv(OUT/'timeseries/heldout_support_bins.csv')
    fig,axs=plt.subplots(2,2,figsize=(9,7),constrained_layout=True)
    im=map_panel(axs[0,0],z['legacy_support'][::-1],x,y,'a  Legacy 28-pair support');im.set_clim(0,1);fig.colorbar(im,ax=axs[0,0],shrink=.7,label='Mean area-pair fraction')
    categories=2*(z['legacy_support']>=.5).astype(int)+z['full_valid'].astype(int)
    colors=['#eeeeee','#56b4e9','#e69f00','#009e73'];names=['Low support / rejected','Low support / accepted','High support / rejected','High support / accepted']
    map_panel(axs[0,1],categories[::-1],x,y,'b  Legacy support versus quality mask',cmap=ListedColormap(colors),norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4))
    axs[0,1].legend(handles=[Patch(facecolor=c,label=n) for c,n in zip(colors,names)],fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.23),ncol=2)
    for name,label in [('legacy28','Legacy 28'),('matched437','Matched input 437')]:
        b=bins[bins.support==name];axs[1,0].plot((b.bin_lower+b.bin_upper)/2,b.final_valid_pixel_fraction,'o-',label=label)
    axs[1,0].set(title='c  Descriptive acceptance by support bin',xlabel='Mean pair support (bin midpoint)',ylabel='Fraction passing full quality mask',xlim=(0,1),ylim=(0,1));axs[1,0].legend()
    comparison=pd.read_csv(OUT/'timeseries/heldout_prediction_comparison.csv')
    for model,label in [('time_series','Time-series differences'),('linear_velocity','Training linear velocity'),('zero_displacement','Zero displacement')]:
        c=comparison[(comparison.model==model)&comparison.group.str.startswith('support_')].sort_values('group')
        axs[1,1].plot(np.arange(.1,1,.2),c.pooled_rmse_mm,'o-',label=label)
    axs[1,1].legend(fontsize=9)
    axs[1,1].set(title='d  Excluded-pair prediction residuals',xlabel='Training-pair support (bin midpoint)',ylabel='Pooled LOS displacement RMSE (mm)',xlim=(0,1))
    save(fig,'F7_matched_processing')
    fig,axs=plt.subplots(2,2,figsize=(9,7),constrained_layout=True)
    fields=[('full_vel_los','a  Unfiltered fitted LOS velocity','mm/year',True),('full_vstd','b  Bootstrap velocity spread','mm/year',False),('full_residual','c  Full-stack inversion residual','mm',False),('heldout_rmse_mm','d  Held-out pixel residual RMSE','mm',False)]
    for ax,(key,title,unit,masked) in zip(axs.ravel(),fields):
        a=z[key].copy()
        if masked:a[~z['full_valid']]=np.nan
        finite=a[np.isfinite(a)]
        im=map_panel(ax,a[::-1],x,y,title,cmap='RdBu_r' if masked else 'magma')
        if masked:
            span=max(float(np.quantile(np.abs(finite),.98)),1);im.set_clim(-span,span)
        else:im.set_clim(0,max(float(np.quantile(finite,.98)),.01))
        fig.colorbar(im,ax=ax,shrink=.7,label=unit,extend='both' if masked else 'max')
    save(fig,'S3_matched_quality')

if __name__=='__main__':main()
