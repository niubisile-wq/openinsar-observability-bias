"""Hidden-class estimation with explicit undefined units and common domains."""
from dataclasses import dataclass
import numpy as np

@dataclass
class Geometry:
    labels: np.ndarray
    count: int
    distance2: np.ndarray
    population: np.ndarray
    known_population: np.ndarray
    valid: np.ndarray

    def aggregate(self, field):
        return np.bincount(self.labels, weights=np.asarray(field).ravel(), minlength=self.count)

def geometry(population, valid, block_pixels, dy=0, dx=0):
    h,w=population.shape
    yy,xx=np.indices((h,w)); gy,gx=(yy+dy)//block_pixels,(xx+dx)//block_pixels
    labels=(gy*int(gx.max()+1)+gx).ravel();n=int(labels.max()+1)
    rmin,rmax=np.full(n,h),np.full(n,-1);cmin,cmax=np.full(n,w),np.full(n,-1)
    np.minimum.at(rmin,labels,yy.ravel());np.maximum.at(rmax,labels,yy.ravel())
    np.minimum.at(cmin,labels,xx.ravel());np.maximum.at(cmax,labels,xx.ravel())
    distance2=(yy.ravel()-(rmin+rmax)[labels]/2)**2+(xx.ravel()-(cmin+cmax)[labels]/2)**2
    known=np.bincount(labels,weights=(population*valid).ravel(),minlength=n)
    return Geometry(labels,n,distance2,population,known,valid)

def retained_statistics(g,rate,missing,include_median=False):
    observed=g.valid&~missing
    obs_pop=g.aggregate(g.population*observed);obs_count=g.aggregate(observed.astype(float))
    nearest_distance=np.full(g.count,np.inf)
    np.minimum.at(nearest_distance,g.labels,np.where(observed.ravel(),g.distance2,np.inf))
    ties=observed.ravel()&(g.distance2==nearest_distance[g.labels])
    result=dict(observed=observed,observed_population=obs_pop,observed_count=obs_count,
        missing_population=g.aggregate(g.population*missing),ties=ties,tie_count=g.aggregate(ties.astype(float)))
    if include_median:
        ids=np.flatnonzero(observed);order=np.lexsort((rate.ravel()[ids],g.labels[ids]));ordered_ids=ids[order]
        units,first,counts=np.unique(g.labels[ordered_ids],return_index=True,return_counts=True)
        sorted_rates=rate.ravel()[ordered_ids];medians=np.full(g.count,np.nan)
        medians[units]=.5*(sorted_rates[first+(counts-1)//2]+sorted_rates[first+counts//2])
        result['medians']=medians
    return result

def evaluate(g,rate,missing,threshold,stats,include_median=False):
    target=rate<=threshold;hidden=stats['missing_population']
    truth=g.aggregate(g.population*missing*target)
    observed_class=g.aggregate(g.population*stats['observed']*target)
    prevalence=np.full(g.count,np.nan)
    np.divide(observed_class,stats['observed_population'],out=prevalence,where=stats['observed_population']>0)
    class_ties=g.aggregate((stats['ties']&target.ravel()).astype(float));nearest=np.full(g.count,np.nan)
    np.divide(class_ties,stats['tie_count'],out=nearest,where=stats['tie_count']>0)
    methods=[('within_block_observed_population_prevalence',prevalence),('nearest_observed_pixel_to_unit_center',nearest)]
    if include_median:
        med=np.where(np.isfinite(stats['medians']),(stats['medians']<=threshold).astype(float),np.nan)
        methods.insert(1,('census_unit_observed_median_rate',med))
    predictions={method:np.where(hidden==0,0.0,fraction*hidden) for method,fraction in methods}
    common=np.logical_and.reduce([np.isfinite(x) for x in predictions.values()])
    total=float(g.known_population.sum());common_total=float(g.known_population[common].sum())
    full_class=float((g.population*g.valid*target).sum());retained_class=float(observed_class.sum())
    lower,upper=retained_class/total*100,(retained_class+hidden.sum())/total*100
    assert lower-1e-10<=full_class/total*100<=upper+1e-10
    def percent(numerator,denominator):
        return float(100*numerator/denominator) if denominator>0 else np.nan
    rows=[]
    for method,pred in predictions.items():
        available=np.isfinite(pred);own_total=float(g.known_population[available].sum())
        own_true=float(truth[available].sum());own_error=pred[available]-truth[available]
        common_error=pred[common]-truth[common];common_true=float(truth[common].sum())
        unknown=float(hidden[~available].sum())
        rows.append(dict(method=method,
            observed_valid_pixels=int(g.valid.sum()),domain_population_allocation_units=total,
            hidden_population_allocation_units=float(hidden.sum()),true_hidden_class_population_allocation_units=float(truth.sum()),
            reference_full_class_population_allocation_units=full_class,positive_population_units=int((g.known_population>0).sum()),
            computable_units=int((available&(g.known_population>0)).sum()),unestimated_units=int((~available&(hidden>0)).sum()),
            unestimated_allocation_units=unknown,unestimated_allocation_share_pct=percent(unknown,total),
            conditional_computable_allocation_units=own_total,conditional_true_hidden_class_units=own_true,
            conditional_estimated_hidden_class_units=float(pred[available].sum()),
            conditional_signed_error_pp=percent(own_error.sum(),own_total),conditional_local_L1_pp=percent(np.abs(own_error).sum(),own_total),
            common_computable_units=int((common&(g.known_population>0)).sum()),common_computable_allocation_units=common_total,
            common_true_hidden_class_units=common_true,common_estimated_hidden_class_units=float(pred[common].sum()),
            common_signed_error_units=float(common_error.sum()),common_signed_error_pp=percent(common_error.sum(),common_total),
            common_local_L1_pp=percent(np.abs(common_error).sum(),common_total),
            common_signed_error_normalized_to_full_pp=percent(common_error.sum(),total),
            common_relative_error_pct=percent(common_error.sum(),common_true),common_relative_error_defined=bool(common_true>0),
            hidden_allocation_share_pct=percent(hidden.sum(),total),whole_domain_lower_bound_pct=float(lower),
            whole_domain_upper_bound_pct=float(upper),whole_domain_bound_width_pp=float(upper-lower),whole_domain_bound_contains_reference=True))
    return rows

def summarize(df,keys):
    return df.groupby(keys,dropna=False).agg(
        cases=('replicate','size'),replicates=('replicate','nunique'),
        mean_unestimated_allocation_pct=('unestimated_allocation_share_pct','mean'),
        min_unestimated_allocation_pct=('unestimated_allocation_share_pct','min'),max_unestimated_allocation_pct=('unestimated_allocation_share_pct','max'),
        mean_common_population=('common_computable_allocation_units','mean'),mean_common_reference=('common_true_hidden_class_units','mean'),
        mean_common_signed_error_units=('common_signed_error_units','mean'),mean_common_signed_error_pp=('common_signed_error_pp','mean'),
        q025_common_signed_error_pp=('common_signed_error_pp',lambda x:x.quantile(.025)),q975_common_signed_error_pp=('common_signed_error_pp',lambda x:x.quantile(.975)),
        mean_common_local_L1_pp=('common_local_L1_pp','mean'),mean_common_relative_error_pct=('common_relative_error_pct','mean'),
        relative_error_defined_cases=('common_relative_error_defined','sum'),mean_lower_bound_pct=('whole_domain_lower_bound_pct','mean'),
        mean_upper_bound_pct=('whole_domain_upper_bound_pct','mean'),mean_bound_width_pp=('whole_domain_bound_width_pp','mean')).reset_index()
