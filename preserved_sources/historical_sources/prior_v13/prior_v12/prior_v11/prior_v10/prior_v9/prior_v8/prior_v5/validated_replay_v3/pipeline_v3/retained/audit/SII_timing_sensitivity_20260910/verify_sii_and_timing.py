import statsmodels.api as sm
from scipy import stats
from threadpoolctl import threadpool_limits
person=pd.read_pickle(A/'SII_person.private.pkl')
timing=pd.read_pickle(A/'timing.private.pkl');timing.person_id=timing.person_id.astype(str)
bio=pd.read_pickle(R/'biomarker_inputs.private.pkl');bio.person_id=bio.person_id.astype(str)
res['male_scale']=res.male_scale.map({1:'01',2:'02'})
M=pd.read_csv(R/'variant_map.tsv',sep='	')
def rn(n):return 'suppressed' if 0<int(n)<20 else int(round(int(n)/100)*100)
summary={'counts_rounded_to_100':{},'validation':{},'timing':{},'selected_results':[]}
for label in ['original_units','strict_units']:
 dd=pd.read_pickle(A/(label+'_timestamps.private.pkl')).reset_index()
 pp=dd.groupby('person_id').SII.count()
 summary['counts_rounded_to_100'][label]={'participants':rn(len(pp)),'timestamps':rn(len(dd)),'timestamps_per_person_median_IQR':pp.quantile([.25,.5,.75]).tolist()}
 ca=dd.merge(bio[['person_id','psoriasis_case','unrelated']],on='person_id',how='inner',validate='many_to_one');ca=ca[ca.psoriasis_case.eq(1)&ca.unrelated].merge(timing,on='person_id',how='left',validate='many_to_one')
 day=ca.dt.dt.tz_localize(None).dt.normalize()
 ca['before_dx']=day<pd.to_datetime(ca.first_dx);ca['on_after_dx']=day>=pd.to_datetime(ca.first_dx)
 ca['before_tx']=day<pd.to_datetime(ca.first_tx);ca['on_after_tx']=day>=pd.to_datetime(ca.first_tx)
 ca['known_dx']=ca.first_dx.notna();ca['known_tx']=ca.first_tx.notna()
 tmp=ca.groupby('person_id')[['before_dx','on_after_dx','before_tx','on_after_tx','known_dx','known_tx']].any()
 summary['timing'][label]={'pruned_cases_with_concurrent_counts':rn(len(tmp)),**{k:rn(tmp[k].sum()) for k in tmp.columns},'percent_cases_any_before_dx':round(100*tmp.before_dx.mean(),1),'percent_cases_any_on_after_dx':round(100*tmp.on_after_dx.mean(),1),'percent_cases_documented_tx':round(100*tmp.known_tx.mean(),1),'percent_tx_cases_any_before_tx':round(100*tmp.loc[tmp.known_tx,'before_tx'].mean(),1),'percent_tx_cases_any_on_after_tx':round(100*tmp.loc[tmp.known_tx,'on_after_tx'].mean(),1)}
for (unit,subset,status),d in res.groupby(['unit_rules','subset','status']):
 summary['counts_rounded_to_100'][unit+'_'+subset+'_'+status]={'model_N_min':rn(d.n.min()),'model_N_max':rn(d.n.max())}
key=['unit_rules','subset','status','male_scale','variant']
assert res.groupby(key).n.nunique().eq(1).all()
summary['validation']['matched_N_identical_all_400_pairs']=True
checks=[]
for scale,variant,outcome in [('02','chr19:10352442:G:C','strict_units'),('02','chr19:10352442:G:C','original_rebuilt'),('02','chrX:5998724:A:G','strict_units'),('01','chrX:5998724:A:G','strict_units')]:
 d=bio[bio.unrelated].merge(person.reset_index(),on='person_id',how='left',validate='one_to_one');d=d[d.strict_units.notna()&d.original_rebuilt.notna()].copy()
 y=np.log(d[outcome]);y=(y-y.mean())/y.std(ddof=1)
 v=M[M.variant_id.eq(variant)].iloc[0];g=d[v.aligned_column].copy()
 if variant.startswith('chrX:') and scale=='01':g.loc[d.sex.eq(1)]/=2
 C=d[['age','sex']+['PC'+str(i) for i in range(1,11)]+['psoriasis_case']].copy();C['g']=g;C=sm.add_constant(C)
 good=np.isfinite(y)&np.isfinite(C).all(axis=1)
 with threadpool_limits(limits=2):fit=sm.OLS(y[good],C.loc[good]).fit()
 rr=res[res.unit_rules.eq('strict_units')&res.subset.eq('pruned')&res.status.eq('full')&res.male_scale.eq(scale)&res.variant.eq(variant)&res.outcome.eq(outcome)].iloc[0]
 checks.append({'N_match':bool(fit.nobs==rr.n),'max_abs_difference':float(max(abs(fit.params['g']-rr.beta),abs(fit.bse['g']-rr.se),abs(fit.pvalues['g']-rr.p),abs(fit.conf_int().loc['g',0]-rr.CI_low),abs(fit.conf_int().loc['g',1]-rr.CI_high)))})
summary['validation']['independent_statsmodels']=checks
sel=res[res.subset.eq('pruned')&res.status.eq('full')&((res.variant.isin(['chr19:10352442:G:C','chr5:132395634:A:G','chr5:132497907:C:T'])&res.male_scale.eq('02'))|(res.variant.str.startswith('chrX:')&res.unit_rules.eq('strict_units')))]
summary['selected_results']=sel.drop(columns=['n','df','rank']).round(10).to_dict('records')
(A/'manuscript_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))