from pathlib import Path
import sys,json,hashlib,datetime,traceback,warnings,importlib.metadata
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.special import expit
from scipy.stats import norm
from threadpoolctl import threadpool_limits
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');R=Path('/home/jupyter/psoriasis_reproducibility_20260910');I=N/'reconciliation_20260912';O=I/'EHR_final_cohort_sensitivity';O.mkdir(exist_ok=True);W=N.parent/'phenotype_validation_20260910';X=N.parent/'X_source_reanalysis_20260910';sources={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def record(p):sources[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size};return p
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def main():
 d=pd.read_pickle(record(W/'private_phenotype_analysis.pkl'));d.index=d.person_id.astype(str);assert d.index.is_unique
 met=pd.read_csv(record(X/'private_provider_qc/genomics_metrics_Dec142023_1859_02_tz0000.tsv'),sep=chr(9),dtype={'research_id':str}).set_index('research_id');assert met.index.is_unique
 covfail=pd.read_csv(record(X/'private_provider/coverage_failures.tsv'),sep=chr(9),dtype=str).research_id
 flags=pd.read_csv(record(X/'private_provider_qc/flagged_samples.tsv'),sep=chr(9),dtype=str).s
 prune=pd.read_csv(record(N.parent/'relatedness_sensitivity_20260910/samples_relatedness_flagged_samples.tsv'),sep=chr(9),dtype=str).sample_id
 ploidy=met.dragen_sex_ploidy.reindex(d.index);canonical=(d.sex.eq(1)&ploidy.eq('XY'))|(d.sex.eq(0)&ploidy.eq('XX'))
 good=~d.index.isin(covfail)&~d.index.isin(flags)&canonical;d=d.loc[good].copy();d['unrelated']=~d.index.isin(prune)
 ph=pd.read_csv(record(N/'inputs/phenotype.tsv'),sep=chr(9),dtype={'IID':str}).set_index('IID');cv=pd.read_csv(record(N/'inputs/covariates.tsv'),sep=chr(9),dtype={'IID':str}).set_index('IID')
 aa=d.loc[d.A_original];assert set(aa.index)==set(ph.index);assert aa.y.reindex(ph.index).eq(ph.psoriasis_case).all()
 for c in cv.columns.drop('FID'):assert np.allclose(aa[c].reindex(cv.index),cv[c],rtol=1e-12,atol=1e-12),c
 cases=set(ph.index[ph.psoriasis_case.eq(1)])
 for definition in ['A_original','B_observed','C_no_drug']:assert set(d.index[d[definition]&d.y.eq(1)])==cases
 G=np.load(record(W/'private_all_sample_lead_genotypes.npy'));H=np.load(record(W/'private_all_sample_lead_heterozygotes.npy'));ps=pd.read_csv(record(W/'private_all_sample_leads.psam'),sep=chr(9),dtype=str)
 sid='IID' if 'IID' in ps else '#IID';assert G.shape==(4,len(ps)) and H.shape==G.shape;assert np.array_equal(ps[sid].iloc[d.source_index.to_numpy(int)].to_numpy(),d.index.to_numpy())
 D=pd.read_pickle(record(R/'reconstructed_x_v2/X_inputs.private.pkl'));D.index=D.index.astype(str)
 leads=[('rs183018930',5998724,'A','G'),('rs5970952',22917619,'C','T'),('rs5964607',67736589,'T','C'),('rs5919457',67842269,'G','T')]
 for j,(rs,_,_,_) in enumerate(leads):
  g=G[j,d.source_index.to_numpy(int)].astype(float);g[(g<0)|(d.sex.eq(1).to_numpy()&H[j,d.source_index.to_numpy(int)])]=np.nan;d[rs]=g
  assert np.allclose(d[rs].reindex(D.index),D[rs],rtol=0,atol=0,equal_nan=True),(rs,'genotype mismatch')
 pairs=pd.read_csv(record(N.parent/'relatedness_sensitivity_20260910/samples_relatedness.tsv'),sep=chr(9),dtype={'i.s':str,'j.s':str});kept=set(d.index[d.unrelated]);assert not (pairs['i.s'].isin(kept)&pairs['j.s'].isin(kept)&pairs.kin.gt(.1)).any()
 qc=pd.read_csv(record(I/'estimate_audit/target_QC.private.tsv.gz'),sep=chr(9));qc=qc.loc[qc.domain.eq('X')&qc.group.eq('pooled')].set_index('ID');assert qc.index.is_unique
 protocol={'scope':'Post-hoc targeted sensitivity; original frozen primary unchanged','definitions':json.loads((W/'phenotype_sensitivity_protocol.json').read_text()),'new_sample_QC':'Apply same provider coverage, flagged-sample, canonical XX/XY concordance exclusions to original and restored controls','relatedness':'all vs fixed provider-pruned; verify no retained released pair kin>0.1','covariates':['intercept','sex','age','age2','age_sex']+['PC'+str(i) for i in range(1,11)]+['sex_PC'+str(i) for i in range(1,11)],'QC':'cases>=200 controls>=500 biological MAC>=20 per case/control and >=100 overall; call rate>=0.98 each','dosage':'female0/1/2 male0/2 or0/1; all effects per coded ALT unit','inference':'nominal Wald tests; four historical targets; all HWE-failed targets diagnostic only; no discovery or replication claim','release':'HELD'};save('PROTOCOL.private.json',protocol);save('input_files.private.json',sources);d.to_pickle(O/'inputs.private.pkl')
 rows=[];cohorts=[]
 with threadpool_limits(limits=2):
  for definition in ['A_original','B_observed','C_no_drug']:
   for subset in ['all','provider_pruned']:
    z=d.loc[d[definition]&(d.unrelated if subset=='provider_pruned' else True)].copy();y0=z.y.to_numpy(float);male=z.sex.eq(1).to_numpy();cohorts.append({'definition':definition,'subset':subset,'n':len(z),'cases':int(y0.sum()),'controls':int((y0==0).sum())})
    pc=['PC'+str(i) for i in range(1,11)];cov=z[['sex','age','age2','age_sex']+pc].to_numpy(float);cov=np.c_[cov,z[pc].to_numpy(float)*male[:,None]]
    for rs,pos,ref,ea in leads:
     g0=z[rs].to_numpy(float);valid=np.isfinite(g0)&np.isfinite(cov).all(1);y=y0[valid];sx=male[valid];g=g0[valid];c=cov[valid];biological=g/np.where(sx,2.,1.);pl=np.where(sx,1.,2.);counts={};why=[]
     for label,v in [('case',1),('control',0)]:
      use=y==v;ac=float(biological[use].sum());an=int(pl[use].sum());counts.update({label+'_AC':ac,label+'_AN':an,label+'_MAC':min(ac,an-ac),label+'_called':int(use.sum()),label+'_callrate':float(valid[y0==v].mean())})
     mac=min(counts['case_AC']+counts['control_AC'],counts['case_AN']+counts['control_AN']-counts['case_AC']-counts['control_AC']);counts['total_MAC']=mac
     if counts['case_called']<200 or counts['control_called']<500:why.append('sample_size')
     if min(counts['case_MAC'],counts['control_MAC'])<20 or mac<100:why.append('MAC')
     if min(counts['case_callrate'],counts['control_callrate'])<.98:why.append('call_rate')
     vid=f'chrX:{pos}:{ref}:{ea}';base={'definition':definition,'subset':subset,'variant':vid,'rsID':rs,'EA':ea,'n':len(y),'cases':int(y.sum()),'controls':int((y==0).sum()),'primary_QC_PASS':bool(qc.loc[vid,'PASS']),'interpretation':'exploratory sensitivity' if bool(qc.loc[vid,'PASS']) else 'diagnostic only; excluded from primary by HWE',**counts}
     for coding in ['male02','male01']:
      rec={**base,'coding':coding}
      if why:rows.append({**rec,'status':'not_estimable','reason':','.join(why)});continue
      sd=c.std(0);assert (sd>0).all();xg=g/np.where(sx,2.,1.) if coding=='male01' else g;xx=np.c_[np.ones(len(y)),(c-c.mean(0))/sd,xg];rank=int(np.linalg.matrix_rank(xx.T@xx));assert rank==xx.shape[1]
      with warnings.catch_warnings(record=True) as caught:
       warnings.simplefilter('always');fit=sm.GLM(y,xx,family=sm.families.Binomial()).fit(maxiter=100,tol=1e-10)
      beta=float(fit.params[-1]);se=float(fit.bse[-1]);mu=expit(xx@fit.params);score=float(np.max(np.abs(xx.T@(y-mu))));se2=float(np.sqrt(np.linalg.inv(xx.T@((mu*(1-mu))[:,None]*xx))[-1,-1]));assert fit.converged and np.isfinite([beta,se]).all() and se>0 and score<1e-5 and abs(se2-se)<1e-7
      rows.append({**rec,'status':'estimated','beta':beta,'se':se,'OR':float(np.exp(beta)),'CI_low':float(np.exp(beta-1.95996398454*se)),'CI_high':float(np.exp(beta+1.95996398454*se)),'p':float(2*norm.sf(abs(beta/se))),'converged':True,'rank':rank,'score_max_abs':score,'SE_independent_check':True,'warnings':[str(w.message) for w in caught]});pd.DataFrame(rows).to_csv(O/'results.private.tsv',sep=chr(9),index=False)
      print('EHR_DONE',definition,subset,rs,coding,flush=True)
 pd.DataFrame(cohorts).to_csv(O/'cohort_counts.private.tsv',sep=chr(9),index=False);save('COMPLETE.private.json',{'status':'validated','models':len(rows),'same_original_strict_participants_and_labels':True,'same_case_set_in_each_definition':True,'source_genotypes_equal_reconstructed_X':True,'pruning_verified':True,'release':'HELD','versions':{k:importlib.metadata.version(k) for k in ['numpy','pandas','statsmodels','scipy']},'script_sha256':sha(Path(__file__))});print('EHR_SENSITIVITY_COMPLETE',flush=True)
try:main()
except Exception as e:save('FAILED.private.json',{'error':repr(e),'traceback':traceback.format_exc(),'release':'HELD'});raise
