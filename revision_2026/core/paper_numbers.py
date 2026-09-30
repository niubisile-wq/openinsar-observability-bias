"""Generate manuscript numerical macros and tables from saved result CSVs."""
import pandas as pd
from common import ROOT,OUT,dump

PAPER=ROOT/'paper'

def main():
    values={};records={}
    def add(key,value,precision=3,source=''):
        values[key]=f'{value:.{precision}f}';records[key]=dict(value=float(value),printed=values[key],source=source)
    base=pd.read_csv(OUT/'fine/allocation_summary.csv')
    def select(domain,support,exposure):return base[(base.domain==domain)&(base.support==support)&(base.exposure==exposure)].iloc[0]
    for q in ['q20','q30','q40']:
        a=select('geographic',q,'population');b=select('geographic',q,'builtup_m2')
        name={'q20':'Twenty','q30':'Thirty','q40':'Forty'}[q]
        add('Pop'+name,a.unsupported/1e6,source='fine/allocation_summary.csv')
        add('Built'+name,b.unsupported/1e6,1,source='fine/allocation_summary.csv')
        add('PopPct'+name,(1-a.weighted_support)*100,2,source='fine/allocation_summary.csv')
    add('PopTotal',a.total/1e6,source='fine/allocation_summary.csv');add('BuiltTotal',b.total/1e6,1,source='fine/allocation_summary.csv')
    scale=pd.read_csv(OUT/'allocation/scale_origin_comparison.csv');u=scale[(scale.domain=='geographic')&(scale.support=='q30')&(scale.exposure=='population')&(scale.method=='area_uniform')&(scale.cell_size_m==1000)]
    for stat in ['min','max']:
        add('UniformOne'+stat.title(),getattr(u.estimated_unsupported,stat)()/1e6,source='allocation/scale_origin_comparison.csv')
        add('UniformBias'+stat.title(),getattr(u.error_fraction_reference,stat)()*100,1,source='allocation/scale_origin_comparison.csv')
    model=pd.read_csv(OUT/'population/model_epoch_comparison.csv');j=model[(model.domain=='geographic')&(model.comparison=='joint_valid')&(model.support=='q30')]
    for name,prefix in [('GHSL_2020_3ss','GHSLJoint'),('WorldPop_2020_3ss','WPJoint'),('GHSL_2015_30ss','EpochOld'),('GHSL_2020_30ss','EpochNew')]:
        r=j[j.model==name].iloc[0];add(prefix+'Pct',r.unsupported_share*100,2,'population/model_epoch_comparison.csv');add(prefix+'Common',r.unsupported_common_total/1e6,3,'population/model_epoch_comparison.csv')
    add('JointExcluded',j.iloc[0].excluded_ghsl_population/1e3,1,'population/model_epoch_comparison.csv')
    sample=pd.read_csv(OUT/'sampling/nested_summary.csv')
    for n,name in [(28,'TwentyEight'),(56,'FiftySix'),(112,'OneTwelve'),(224,'TwoTwentyFour')]:
        r=sample[(sample.pairs==n)&(sample.support=='q30')&(sample.exposure=='population')].iloc[0];add('Sample'+name,r.unsupported/1e6,3,'sampling/nested_summary.csv')
    dwr=pd.read_csv(OUT/'dwr/summary.csv');r=dwr[dwr.exposure=='population'].iloc[0]
    add('DWRTotal',r.total,0,'dwr/summary.csv');add('DWRUnsupported',r.unsupported,0,'dwr/summary.csv');add('DWRSupportPct',r.support_share*100,2,'dwr/summary.csv')
    du=pd.read_csv(OUT/'dwr/scale_origin_comparison.csv');du=du[(du.exposure=='population')&(du.method=='area_uniform')&(du.cell_size_m==1000)]
    add('DWRUniformMin',du.estimated_unsupported.min(),0,'dwr/scale_origin_comparison.csv');add('DWRUniformMax',du.estimated_unsupported.max(),0,'dwr/scale_origin_comparison.csv')
    (PAPER/'numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in values.items())+'\n',encoding='utf-8')
    dump(PAPER/'number_ledger.json',records)
    rows=[]
    for q,label in [('unw','Valid unwrapped'),('q20',r'Joint, $q=0.2$'),('q30',r'Joint, $q=0.3$'),('q40',r'Joint, $q=0.4$')]:
        a=select('geographic',q,'population');b=select('geographic',q,'builtup_m2');c=select('geographic',q,'area_m2')
        rows.append(f'{label} & {a.unsupported/1e6:.3f} & {100*(1-a.weighted_support):.2f} & {b.unsupported/1e6:.1f} & {100*(1-c.weighted_support):.2f} '+r'\\')
    (PAPER/'baseline_table.tex').write_text(r'\begin{tabular}{lrrrr}\toprule'+'\n'+r'Support rule & Population & Pop. (\%) & Built surface & Area (\%)\\\midrule'+'\n'+'\n'.join(rows)+'\n'+r'\bottomrule\end{tabular}'+'\n',encoding='utf-8')
    print('PAPER NUMBERS',len(values),flush=True)

if __name__=='__main__':main()
