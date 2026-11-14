import re, traceback
from scipy import stats
from threadpoolctl import threadpool_limits
state={'stage':'starting'}
def run_sii():
 global sii_person, sii_draws, sii_results, audit, strict_map
 try:
  state['stage']='cleaning counts'
  c=labs[labs.biomarker.isin(['ALC','ANC','platelet_count'])].copy()
  c['person_id']=c.person_id.astype(str)
  c['value_as_number']=pd.to_numeric(c.value_as_number,errors='coerce')
  u=c.unit_name.astype(str).str.lower()+' | '+c.unit_source_value.astype(str).str.lower()
  one=u.str.contains(r'/ul|/µl|/mcl|per microliter|cells/ul|cells per microliter',case=False,na=False,regex=True)
  thousand=u.str.contains('|'.join(map(re.escape,['10*3','10^3','thousand','k/ul','10*9/l','billion per liter'])),case=False,na=False,regex=True)
  c['original_std']=c.value_as_number/np.where(one&~thousand,1000.,1.)
  standard={'billion per liter':1.,'thousand per microliter':1.,'thousand per cubic millimeter':1.,'x10(3)/mcl':1.,'cells per microliter':.001,'cells/ul':.001,'per microliter':.001,'/mm3':.001,'per cubic millimeter':.001,'cells per cubic millimeter':.001,'million per microliter':1000.,'ten thousand per microliter':10.,'billion per microliter':1000000.}
  source={'10*9/l':1.,'10*3/ul':1.,'10*3/mm3':1.,'k/ul':1.,'k/mm3':1.,'{cells}/ul':.001,'/ul':.001,'/mm3':.001,'{cells}/mm3':.001,'10*6/ul':1000.,'10*4/ul':10.,'10*9/ul':1000000.}
  a=c.unit_name.str.lower().map(standard);b=c.unit_source_value.astype(str).str.lower().map(source)
  fac=a.where(a.notna()&(b.isna()|np.isclose(a,b,equal_nan=False)))
  c['strict_std']=c.value_as_number*fac
  strict_map=c[['unit_name','unit_source_value']].copy();strict_map['factor']=fac
  strict_map=strict_map.drop_duplicates();strict_map.to_csv(A/'unit_mapping.tsv',sep='	',index=False)
  c['dt']=pd.to_datetime(c.measurement_datetime,errors='coerce',utc=True,format='mixed')
  limits={'ALC':(.01,100),'ANC':(.01,100),'platelet_count':(1,2000)}
  def cleaned(col):
   good=c.value_as_number.gt(0)&np.isfinite(c[col])
   for k,(lo,hi) in limits.items():good &= (~c.biomarker.eq(k))|c[col].between(lo,hi)
   return c.loc[good,['person_id','biomarker','dt',col]].rename(columns={col:'v'})
  co=cleaned('original_std');cs=cleaned('strict_std')
  med=co.groupby(['person_id','biomarker']).v.median().unstack()
  med['original_rebuilt']=med.platelet_count*med.ANC/med.ALC
  oo=orig.copy();oo.person_id=oo.person_id.astype(str);oo=oo.set_index('person_id')
  check=med.join(oo[['ALC','ANC','platelet_count','inflammatory_index']],rsuffix='_saved')
  audit={'original_components_match':{k:bool(np.allclose(check[k],check[k+'_saved'],equal_nan=True,rtol=1e-10,atol=1e-10)) for k in ['ALC','ANC','platelet_count']},'original_formula_match':bool(np.allclose(check.original_rebuilt,check.inflammatory_index,equal_nan=True,rtol=1e-10,atol=1e-10))}
  state['stage']='matching exact timestamps'
  sii_draws={};sii_person=med[['original_rebuilt']].copy()
  for label,cc in [('original_units',co),('strict_units',cs)]:
   cc=cc[cc.dt.notna()&cc.dt.ne(cc.dt.dt.normalize())].drop_duplicates()
   gg=cc.groupby(['person_id','dt','biomarker']).v.agg(['nunique','first'])
   wide=gg['first'].unstack('biomarker');counts=gg['nunique'].unstack('biomarker')
   wide=wide.loc[counts.eq(1).all(axis=1)].dropna(subset=['ALC','ANC','platelet_count'])
   wide['SII']=wide.platelet_count*wide.ANC/wide.ALC
   wide=wide[np.isfinite(wide.SII)&wide.SII.gt(0)]
   sii_draws[label]=wide
   per=wide.SII.groupby(level='person_id').median()
   sii_person[label]=per
   audit[label]={'participants_private':len(per),'timestamps_private':len(wide)}
  sii_person.to_pickle(A/'SII_person.private.pkl')
  for k,v in sii_draws.items():v.to_pickle(A/(k+'_timestamps.private.pkl'))
  state['stage']='fitting matched regressions'
  bd=bio.copy();bd.person_id=bd.person_id.astype(str)
  bd=bd.merge(sii_person.reset_index(),on='person_id',how='left',validate='one_to_one')
  M=pd.read_csv(R/'variant_map.tsv',sep='	');pcs=['PC'+str(i) for i in range(1,11)]
  runner=(R/'run_biomarkers.py').read_text();fun=runner[runner.index('def ols('):runner.index('br=[]')]
  exec(fun,globals())
  out=[]
  with threadpool_limits(limits=2):
   for unit in ['original_units','strict_units']:
    for subset in ['all','pruned']:
     for status in ['full','cases']:
      d=bd.loc[(bd.unrelated if subset=='pruned' else np.ones(len(bd),bool))&(bd.psoriasis_case.eq(1) if status=='cases' else True)].copy()
      d=d[d[unit].notna()&d.original_rebuilt.notna()]
      cols=['age','sex']+pcs+(['psoriasis_case'] if status=='full' else [])
      for outcome in ['original_rebuilt',unit]:
       vals=np.log(d[outcome]);d['Y_'+outcome]=(vals-vals.mean())/vals.std(ddof=1)
      for scale in ['02','01']:
       for outcome in ['original_rebuilt',unit]:
        family=[]
        for _,v in M.iterrows():
         g=d[v.aligned_column].to_numpy(float).copy()
         if str(v.variant_id).startswith('chrX:') and scale=='01':g[d.sex.eq(1).to_numpy()]/=2
         fit=ols(d['Y_'+outcome].to_numpy(float),g,d[cols].to_numpy(float))
         family.append({'unit_rules':unit,'subset':subset,'status':status,'male_scale':scale,'outcome':outcome,'variant':v.variant_id,'EA':v.effect_allele,**fit})
        p=np.array([r['p'] for r in family]);order=np.argsort(p);q=np.empty(25);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*25/np.arange(1,26))[::-1])[::-1])
        for r,qq in zip(family,q):r['q25']=float(qq);r['bonf25']=r['p']<.002;r['historical_threshold']=r['p']<.0005
        out+=family
  sii_results=pd.DataFrame(out);sii_results.to_csv(A/'matched_results.private.tsv',sep='	',index=False)
  audit['models']=len(out);audit['correlations']={k:float(np.log(sii_person[['original_rebuilt',k]]).corr().iloc[0,1]) for k in ['original_units','strict_units']}
  (A/'audit.private.json').write_text(json.dumps(audit,indent=2));state['stage']='complete'
 except Exception:
  state['stage']='failed';state['error']=traceback.format_exc();(A/'error.txt').write_text(state['error'])
worker=threading.Thread(target=run_sii);worker.start()
print('SII_WORKER_STARTED')
