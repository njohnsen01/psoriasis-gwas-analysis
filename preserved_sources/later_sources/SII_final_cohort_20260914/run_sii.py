from pathlib import Path
import os,json,sys,hashlib,datetime,subprocess,traceback
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
R=Path('/home/jupyter/psoriasis_reproducibility_20260910')
O=B/'reconciliation_20260912/SII_final_cohort_20260914';O.mkdir(exist_ok=True)
def save(n,x): (O/n).write_text(json.dumps(x,indent=2,default=str))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
try:
 L=pd.read_csv(B/'submission_assets_20260911/pooled_autosomal_reporting_groups.private.tsv',sep='	')
 L=L.loc[~((L.CHROM==6)&L.GENPOS.between(25000000,34000000))].copy();assert L.ID.is_unique
 L.to_csv(O/'final_variant_list.tsv',sep='	',index=False)
 protocol={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'fixed_before_new_fits','scope':'post-discovery sensitivity analysis, not independent replication','variants':L.ID.tolist(),'definitions':['original_rebuilt separate component medians','strict concurrent timestamp SII followed by participant median'],'sample':'final keep, provider-pruned, all five outcomes and covariates complete; variant-specific genotype completeness identical across outcomes','scaling':'natural logs; both SII definitions use mean and sample SD of concurrent SII in common baseline sample; each component uses own baseline log mean/SD; fixed for all ancestry and outlier models','components':'participant medians from same valid concurrent timestamp draws; component models are not an exact algebraic decomposition of median SII','covariates':['age','age2','SEX','age_sex','PC1-PC10','psoriasis_case'],'relatedness':'existing fixed provider-pruned set, require no retained provider pair kinship>0.1; distant or unavailable relationships not excluded by this guarantee','estimation':'OLS with HC3 covariance and asymptotic two-sided normal tests; rank checks, drop constant nuisance columns only','families':{'main':5*len(L),'robustness':20*len(L),'ancestry_heterogeneity':5*len(L)},'multiplicity':'Bonferroni and BH with all planned hypotheses, unavailable tests count as P=1','outliers':'winsorize log outcomes at fixed pooled baseline 1st and 99th percentiles, then apply original common scaling','ancestry':['EUR','AFR','AMR'],'estimability':'N>=500, nonmissing diploid minor allele count>=20, full design rank; not a power claim','excluded':'no treatment or diagnosis timing analyses'}
 save('PROTOCOL.json',protocol)
 sys.path.insert(0,str(B/'scripts'));from run_autosomes import download
 import zstandard
 tasks=json.loads((B/'autosomal_tasks.json').read_text());need=set(L.ID);found={};G=[];receipts=[]
 for t in tasks:
  if int(t['chromosome']) not in set(L.CHROM):continue
  key='chr'+t['chromosome']+'_'+t['stem'];p=B/'results/autosomes'/key/'passed_variants.txt'
  if not p.exists():continue
  hits=need.intersection(p.read_text().split())
  if not hits:continue
  print('EXTRACT',key,'variants',len(hits),flush=True)
  w=O/'source_work'/key;w.mkdir(parents=True,exist_ok=True)
  for ext in ['pgen','pvar','psam']:
   target=w/('input.pvar.zst' if ext=='pvar' else 'input.'+ext);receipts.append(download(t[ext],target))
   if ext=='pvar':
    with target.open('rb') as fi,(w/'input.pvar').open('wb') as fo:zstandard.ZstdDecompressor().copy_stream(fi,fo)
  v=pd.read_csv(w/'input.pvar',sep='	',comment='#',header=None)
  # Verify with explicit header parsing instead of assuming metadata lines.
  with (w/'input.pvar').open() as f:
   header=next(x for x in f if x.startswith('#CHROM')).strip().lstrip('#').split('	')
  v.columns=header;z=v.set_index('ID').loc[sorted(hits)]
  for vid,row in z.iterrows():
   expected=L.set_index('ID').loc[vid];assert str(row.ALT)==str(expected.ALLELE1) and str(row.REF)==str(expected.ALLELE0)
  (w/'extract.txt').write_text('\n'.join(sorted(hits))+'\n')
  (w/'alleles.txt').write_text('\n'.join(vid+'	'+str(z.loc[vid,'ALT']) for vid in sorted(hits))+'\n')
  cmd=['plink2','--pfile',str(w/'input'),'--keep',str(B/'inputs/keep.tsv'),'--extract',str(w/'extract.txt'),'--export','A','--export-allele',str(w/'alleles.txt'),'--threads','2','--memory','2048','--out',str(w/'dosage')]
  with (w/'export.console.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
  d=pd.read_csv(w/'dosage.raw',sep='\s+',dtype={'IID':str}).set_index('IID');out=pd.DataFrame(index=d.index)
  for vid in sorted(hits):
   col=vid+'_'+str(z.loc[vid,'ALT']);assert col in d;out[vid]=pd.to_numeric(d[col],errors='raise');assert out[vid].dropna().between(0,2).all();found[vid]=key
  G.append(out)
 assert set(found)==need, 'Missing exact finalized genotypes'
 g=pd.concat(G,axis=1);assert g.index.is_unique and g.columns.is_unique;g.to_pickle(O/'final_genotypes.private.pkl');save('source_receipts.private.json',receipts)
 print('GENOTYPES_READY',len(found),flush=True)
 p=R/'reconstructed_sii_v2';bio=pd.read_pickle(p/'biomarker_inputs.private.pkl');bio.index=bio.index.astype(str)
 anc=pd.read_pickle(R/'reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');anc.index=anc.index.astype(str)
 keep=pd.read_csv(B/'inputs/keep.tsv',sep='\s+',dtype=str);iid='IID' if 'IID' in keep else keep.columns[-1];ids=set(keep[iid])
 cov=pd.read_csv(B/'inputs/covariates.tsv',sep='\s+',dtype={'IID':str}).set_index('IID');ph=pd.read_csv(B/'inputs/phenotype.tsv',sep='\s+',dtype={'IID':str}).set_index('IID')
 cov.columns=[('SEX' if x=='sex' else x) for x in cov.columns]
 cols=['age','age2','SEX','age_sex']+['PC'+str(i) for i in range(1,11)]
 base=cov.loc[cov.index.isin(ids),cols].join(ph[['psoriasis_case']]).join(bio[['person_id','unrelated']]).join(anc[['ancestry_pred']]);assert base.index.is_unique
 base=base.loc[base.unrelated.eq(True)].copy()
 pairs=pd.read_csv(B.parent/'relatedness_sensitivity_20260910/samples_relatedness.tsv',sep='\s+',dtype={'i.s':str,'j.s':str})
 assert not ((pairs.kin>.1)&pairs['i.s'].isin(base.index)&pairs['j.s'].isin(base.index)).any()
 sii=pd.read_pickle(p/'SII_person.private.pkl');sii.index=sii.index.astype(str)
 draws=pd.read_pickle(p/'strict_units_timestamps.private.pkl');comp=draws.groupby(level='person_id')[['platelet_count','ANC','ALC']].median();comp.index=comp.index.astype(str)
 outcomes=['original_rebuilt','strict_units','platelet_count','ANC','ALC'];Y=sii[['original_rebuilt','strict_units']].join(comp)
 base['person_key']=base.person_id.astype('Int64').astype(str);base=base.join(Y,on='person_key');base=base.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+['psoriasis_case']+outcomes);base=base.loc[base[outcomes].gt(0).all(axis=1)].copy();base=base.join(g)
 assert len(base)>=500
 scaling={};winsor={}
 for y in outcomes:
  ly=np.log(base[y]);ref=np.log(base.strict_units) if y in outcomes[:2] else ly
  scaling[y]={'mean':float(ref.mean()),'sd':float(ref.std(ddof=1))};winsor[y]=ly.quantile([.01,.99]).tolist()
 save('scaling.private.json',{'scaling':scaling,'winsor_bounds':winsor});base.to_pickle(O/'matched_inputs.private.pkl')
 from scipy.stats import norm,chi2
 from statsmodels.stats.multitest import multipletests
 rows=[]
 for mode in ['baseline','winsorized','EUR','AFR','AMR']:
  a=base if mode in ['baseline','winsorized'] else base.loc[base.ancestry_pred.eq(mode)]
  for vid in L.ID:
   d=a.dropna(subset=[vid]);n=len(d);mac=min(d[vid].sum(),2*n-d[vid].sum());nuis=cols+['psoriasis_case'];nuis=[c for c in nuis if d[c].nunique()>1]
   X=np.column_stack([np.ones(n),d[nuis].to_numpy(float),d[vid].to_numpy(float)])
   ok=n>=500 and mac>=20 and np.linalg.matrix_rank(X)==X.shape[1]
   if ok:
    inv=np.linalg.inv(X.T@X);hat=np.einsum('ij,jk,ik->i',X,inv,X);assert hat.max()<1
   for y in outcomes:
    row={'mode':mode,'variant':vid,'effect_allele':str(L.set_index('ID').loc[vid,'ALLELE1']),'outcome':y,'N':n,'MAC':mac,'covariates':','.join(nuis),'status':'estimated' if ok else 'not_estimable'}
    if ok:
     raw=np.log(d[y].to_numpy(float));raw=np.clip(raw,*winsor[y]) if mode=='winsorized' else raw
     yy=(raw-scaling[y]['mean'])/scaling[y]['sd'];beta=inv@X.T@yy;res=yy-X@beta;meat=X.T@((res/(1-hat))[:,None]**2*X);vc=inv@meat@inv;se=np.sqrt(vc[-1,-1]);b=beta[-1];pv=2*norm.sf(abs(b/se));row.update(beta=b,se=se,CI_low=b-1.95996398454*se,CI_high=b+1.95996398454*se,p=pv)
    rows.append(row)
  print('MODE_COMPLETE',mode,flush=True)
 res=pd.DataFrame(rows)
 for label,mask,m in [('main',res['mode'].eq('baseline'),5*len(L)),('robustness',res['mode'].ne('baseline'),20*len(L))]:
  pv=res.loc[mask,'p'].fillna(1).to_numpy();assert len(pv)==m;res.loc[mask,'family']=label;res.loc[mask,'p_bonferroni']=np.minimum(1,pv*m);res.loc[mask,'q_BH']=multipletests(pv,method='fdr_bh')[1]
 res.to_csv(O/'SII_component_results.private.tsv',sep='	',index=False)
 het=[]
 for vid in L.ID:
  for y in outcomes:
   z=res.loc[res['mode'].isin(['EUR','AFR','AMR'])&res.variant.eq(vid)&res.outcome.eq(y)&res.status.eq('estimated')];h={'variant':vid,'outcome':y,'groups':','.join(z['mode']),'status':'estimated' if len(z)>=2 else 'not_estimable'}
   if len(z)>=2:
    w=1/z.se**2;mu=(w*z.beta).sum()/w.sum();Q=(w*(z.beta-mu)**2).sum();h.update(Q=Q,df=len(z)-1,p=chi2.sf(Q,len(z)-1))
   het.append(h)
 h=pd.DataFrame(het);pv=h.p.fillna(1);h['p_bonferroni']=np.minimum(1,pv*len(h));h['q_BH']=multipletests(pv,method='fdr_bh')[1];h.to_csv(O/'ancestry_heterogeneity.private.tsv',sep='	',index=False)
 save('COMPLETE.json',{'status':'complete','variants':len(L),'planned_model_tests':len(res),'estimated_models':int(res.status.eq('estimated').sum()),'main_bonferroni_significant':int((res.family.eq('main')&res.p_bonferroni.lt(.05)).sum()),'robustness_bonferroni_significant':int((res.family.eq('robustness')&res.p_bonferroni.lt(.05)).sum()),'protocol_sha256':sha(O/'PROTOCOL.json'),'outputs_remain_private':True,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 print('SII_ANALYSIS_COMPLETE',flush=True)
except Exception:
 save('FAILED.json',{'traceback':traceback.format_exc()});traceback.print_exc();raise
