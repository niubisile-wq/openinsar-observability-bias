"""Scientific plots from the retained six-experiment output tables."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'savefig.dpi':200,'font.family':'DejaVu Sans'})
COLORS=['#2368A0','#D87928','#339579']
REGION_LABELS={'Los_Angeles':'Los Angeles','Davis_Sacramento':'Davis-Sacramento'}


def save(fig,out,name):
    fig.savefig(out/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(out/(name+'.png'),bbox_inches='tight')
    plt.close(fig)


def main(results,out):
    out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(results/'real_block_masking/real_block_masking_summary.csv')
    d=d[d.removed_valid_pixel_fraction.eq(.3)&d.threshold_mm_per_year.eq(-5)]
    patterns=['random_cells','one_contiguous_hole','four_contiguous_holes']
    methods=['native_observed_allocation_fraction','area_observed_class_fraction','published_pixel_center_median_class']
    available=d.method.unique()
    methods=[m if m in available else next(x for x in available if 'median' in x) for m in methods]
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    x=np.arange(3)
    for j,(method,label) in enumerate(zip(methods,['Native prevalence','Area prevalence','Unit median rate'])):
        s=d[d.method.eq(method)].set_index('mask_geometry').loc[patterns]
        axs[0].bar(x+(j-1)*.24,s.mean_unestimated_allocation_pct,.23,color=COLORS[j],label=label)
        axs[1].bar(x+(j-1)*.24,s.mean_common_local_L1_pp,.23,color=COLORS[j])
    for ax in axs:
        ax.set_xticks(x,['Random pixels','One large hole','Four large holes']);ax.tick_params(axis='x',labelsize=9)
    axs[0].set_ylabel('Unestimated population allocation (%)');axs[0].set_ylim(0,36)
    axs[1].set_ylabel('Local L1 error on common-computable blocks (pp)')
    axs[0].set_title('(a) Availability');axs[1].set_title('(b) Conditional estimation error')
    axs[0].legend(fontsize=8,loc='upper left')
    save(fig,out,'additional_S25_real_blocks')

    d=pd.read_csv(results/'spatial_acquisition/spatial_calibration_summary.csv')
    fig,axs=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for j,model in enumerate(['GHSL','WorldPop']):
        for axis,marker,ls in [('west_east','o','-'),('south_north','s','--')]:
            s=d[d.model.eq(model)&d.axis.eq(axis)].sort_values('buffer_m')
            axs[0].plot(s.buffer_m/1000,s.WAPE_pct,marker=marker,ls=ls,color=COLORS[j],
                        label=f'{model}, '+('W-E bands' if axis=='west_east' else 'S-N bands'))
    axs[0].set_xlabel('Training exclusion buffer (km)');axs[0].set_ylabel('Held-out WAPE against ACS (%)')
    axs[0].set_title('(a) Spatial recalibration');axs[0].legend(fontsize=8)
    a=pd.read_csv(results/'spatial_acquisition/acquisition_date_influence.csv')
    a=a[a.support.eq('q30')&a.exposure.eq('population')].copy()
    # Rows in source are sorted by deleted date, a fixed chronological index.
    names=[c for c in a.columns if 'quarter' in c and 'change' in c]
    if len(names)!=1:raise RuntimeError(list(a.columns))
    axs[1].scatter(np.arange(len(a)),a[names[0]],s=12,color=COLORS[0],alpha=.8)
    axs[1].axhline(0,color='#888888',lw=.8);axs[1].set_xlabel('Deleted acquisition date (chronological index)')
    axs[1].set_ylabel('Equal-quarter population deficit change (pp)')
    axs[1].set_title('(b) Shared-date influence; no re-inversion')
    save(fig,out,'additional_S26_dependence')

    q=pd.read_csv(results/'training_quality/quality_coverage_error.csv')
    s=q[q.training_metric.eq('coh_avg_min')].sort_values('threshold')
    fig,axs=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for ax,metric,label in zip(axs,['pooled_rmse_mm','population_weighted_pooled_rmse_mm'],
                             ['Pooled excluded-edge RMSE (mm)','Population-weighted excluded-edge RMSE (mm)']):
        ax.plot(s.accepted_domain_population_pct,s[metric],'-o',color=COLORS[0],label='Training coherence rule')
        for _,r in s.iterrows():
            if r.threshold in [.2,.3,.4,.7]:ax.annotate(f'q={r.threshold:g}',(r.accepted_domain_population_pct,r[metric]),
                                                       xytext=(7,-13) if r.threshold==.4 else (5,6),textcoords='offset points',fontsize=8)
        for rule,marker,col,txt in [('finite_training_velocity','x','#666666','Finite velocity'),
                                    ('fixed_training_default_mask','s',COLORS[1],'Default training mask')]:
            r=q[q.rule.eq(rule)].iloc[0];ax.scatter(r.accepted_domain_population_pct,r[metric],marker=marker,color=col,s=45,label=txt)
        ax.set_xlabel('Accepted full-domain population allocation (%)');ax.set_ylabel(label);ax.set_xlim(37,105)
    axs[0].set_title('(a) Pixel/edge pooled residual');axs[1].set_title('(b) Population-weighted residual')
    axs[0].legend(fontsize=8)
    save(fig,out,'additional_S27_quality_curve')

    d=pd.read_csv(results/'new_regions/new_region_allocation_summary.csv')
    fig,axs=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for ax,region in zip(axs,['Los_Angeles','Davis_Sacramento']):
        for j,threshold in enumerate([.3,.4]):
            s=d[d.region.eq(region)&d.coherence_threshold.eq(threshold)].groupby('scale_m').uniform_minus_native_pp.agg(['min','mean','max'])
            xx=s.index.to_numpy();ax.plot(xx,s['mean'],marker='o',color=COLORS[j],label=f'Coherence >= {threshold:g}')
            ax.fill_between(xx,s['min'],s['max'],color=COLORS[j],alpha=.18)
        ax.set_xlabel('Reporting scale (m)');ax.set_ylabel('Uniform minus native deficit (pp)')
        ax.set_xticks([250,500,1000,2000]);ax.set_title(REGION_LABELS[region]);ax.legend(fontsize=8)
    save(fig,out,'additional_S28_new_regions')

    e=pd.read_csv(results/'residential_allocation/residential_landuse_eligibility.csv')
    e=e[e.residential_definition.eq('strict')].set_index('population_allocation_model').loc[['Area_uniform','GHSL_100m','CA_POP_100m']]
    fig,axs=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    bottom=np.zeros(3)
    for col,color,label in [('residential_share_of_all_units_pct',COLORS[0],'Residential eligibility'),
                            ('mapped_nonresidential_share_of_all_units_pct',COLORS[1],'Mapped other land use'),
                            ('unmapped_share_of_all_units_pct','#BBBBBB','Unmapped classification')]:
        axs[0].bar(np.arange(3),e[col],bottom=bottom,color=color,label=label);bottom+=e[col].to_numpy()
    axs[0].set_xticks(np.arange(3),['Area uniform','GHSL','CA-POP']);axs[0].set_ylabel('Complete allocation denominator (%)');axs[0].set_ylim(0,118)
    axs[0].legend(fontsize=8,loc='upper left');axs[0].set_title('(a) Full-domain eligibility diagnostic')
    d=pd.read_csv(results/'residential_allocation/residential_allocation_summary.csv')
    for j,definition in enumerate(['strict','broader_mixed_use']):
        s=d[d.residential_definition.eq(definition)]
        axs[1].bar(np.arange(6)+(j-.5)*.34,s.support_deficit_share_pct,.33,color=COLORS[j],label=definition.replace('_',' '))
    axs[1].set_xticks(np.arange(6),['Area','GHSL','CA-POP','Res. area','GHSL + res.','CA-POP + res.'],rotation=30,ha='right',fontsize=8)
    axs[1].set_ylabel('Common-domain support deficit (%)');axs[1].set_title('(b) Restricted domain: 89.13-89.22% of Census total')
    axs[1].legend(fontsize=8)
    save(fig,out,'additional_S29_residential')

    d=pd.read_csv(results/'gnss_transfer/gnss_transfer_summary.csv')
    fig,axs=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for ax,region,refs,n in zip(axs,['Los_Angeles','Davis_Sacramento'],[['P597','MHMS'],['P266','UCD1']],[11,4]):
        for j,reference in enumerate(refs):
            s=d[d.region.eq(region)&d.frame_control.eq(reference)].sort_values('pixel_radius')
            ax.plot(s.pixel_radius,s.RMSE_mm_yr,marker='o',color=COLORS[j],label=f'Frame control {reference}')
        ax.set_xticks([0,1,3]);ax.set_xlabel('Sampling window radius (native pixels)');ax.set_ylabel('Relative LOS RMSE against archived GNSS (mm/yr)')
        ax.set_title(REGION_LABELS[region]+f' ({n} validation stations)');ax.set_ylim(0,16);ax.legend(fontsize=8)
    save(fig,out,'additional_S30_GNSS')
    print('Saved six scientific figures as PDF and PNG',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=ROOT/'results');p.add_argument('--output',type=Path,default=ROOT/'02_Editable_Sources/figures');a=p.parse_args();main(a.results,a.output)
