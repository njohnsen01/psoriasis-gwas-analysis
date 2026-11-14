from pathlib import Path
import argparse,json,hashlib
import pandas as pd,numpy as np
N=Path(__file__).resolve().parents[1];B=N.parent;R=Path('/home/jupyter/psoriasis_reproducibility_20260910');X=B/'X_source_reanalysis_20260910'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(out,reference):
 out.mkdir(parents=True,exist_ok=False)
 inputs=[R/'replay_cohort_v2/psoriasis_pheno.tsv',R/'replay_cohort_v2/psoriasis_covariates.tsv',X/'private_provider_qc/genomics_metrics_Dec142023_1859_02_tz0000.tsv',X/'private_provider/coverage_failures.tsv',X/'private_provider_qc/flagged_samples.tsv',B/'phenotype_validation_20260910/private_full_ancestry.tsv']
 ph=pd.read_csv(inputs[0],sep='\t',dtype={'FID':str,'IID':str});cv=pd.read_csv(inputs[1],sep='\t',dtype={'FID':str,'IID':str})
 metrics=pd.read_csv(inputs[2],sep='\t',dtype={'research_id':str}).set_index('research_id')
 coverage=pd.read_csv(inputs[3],sep='\t',dtype=str);flags=pd.read_csv(inputs[4],sep='\t',dtype=str)
 ancestry=pd.read_csv(inputs[5],sep='\t',dtype={'research_id':str}).set_index('research_id')
 assert ph.IID.is_unique and cv.IID.is_unique and metrics.index.is_unique and ancestry.index.is_unique
 assert len(ph)==392155 and ph.psoriasis_case.sum()==5805
 q=cv.set_index('IID').reindex(ph.IID);m=metrics.reindex(ph.IID)
 canonical=m.dragen_sex_ploidy.isin(['XX','XY']);concordant=(q.sex.eq(0)&m.dragen_sex_ploidy.eq('XX'))|(q.sex.eq(1)&m.dragen_sex_ploidy.eq('XY'))
 valid=canonical&concordant&~ph.IID.isin(coverage.research_id).to_numpy()&~ph.IID.isin(flags.s).to_numpy()
 p=ph.loc[valid.to_numpy(),['FID','IID','psoriasis_case']].copy();c=cv.set_index('IID').loc[p.IID].reset_index()
 cols=['FID','IID','sex','age','age2','age_sex']+['PC'+str(i) for i in range(1,11)];c=c[cols]
 assert len(p)==386069 and p.psoriasis_case.sum()==5705 and np.isfinite(c.drop(columns=['FID','IID']).to_numpy(float)).all()
 p.to_csv(out/'phenotype.tsv',sep='\t',index=False);c.to_csv(out/'covariates.tsv',sep='\t',index=False)
 p[['FID','IID']].to_csv(out/'keep.tsv',sep='\t',header=False,index=False)
 for name,subset in [('cases',p.loc[p.psoriasis_case.eq(1)]),('controls',p.loc[p.psoriasis_case.eq(0)])]:subset[['FID','IID']].to_csv(out/(name+'.keep.tsv'),sep='\t',header=False,index=False)
 labels=ancestry.ancestry_pred.reindex(p.IID).str.upper();assert labels.notna().all()
 for group in ['EUR','AFR','AMR']:
  subset=p.loc[labels.eq(group).to_numpy()&p.psoriasis_case.eq(0).to_numpy()];subset[['FID','IID']].to_csv(out/(group+'.controls.keep.tsv'),sep='\t',header=False,index=False)
 checks={}
 if reference:
  for name in ['phenotype.tsv','covariates.tsv']:
   got=pd.read_csv(out/name,sep='\t',dtype={'FID':str,'IID':str}).set_index('IID').sort_index();want=pd.read_csv(reference/name,sep='\t',dtype={'FID':str,'IID':str}).set_index('IID').sort_index()
   pd.testing.assert_frame_equal(got,want[got.columns],check_dtype=False,check_exact=False,rtol=1e-12,atol=1e-12);checks[name]='all participant keys and values match retained primary inputs'
  for name in ['keep.tsv','cases.keep.tsv','controls.keep.tsv','EUR.controls.keep.tsv','AFR.controls.keep.tsv','AMR.controls.keep.tsv']:
   got=pd.read_csv(out/name,sep='\t',header=None,dtype=str);want=pd.read_csv(reference/name,sep='\t',header=None,dtype=str)
   assert set(map(tuple,got.to_numpy()))==set(map(tuple,want.to_numpy()));checks[name]='exact participant key sets match'
 receipt={'status':'validated' if reference else 'constructed','inputs':{str(x):sha(x) for x in inputs},'outputs':{x.name:sha(x) for x in out.glob('*.tsv')},'comparison':checks,'script_sha256':sha(Path(__file__)),'release_status':'participant files must remain in Workbench'}
 (out/'VALIDATION.private.json').write_text(json.dumps(receipt,indent=2));print('PRIMARY_INPUT_BOOTSTRAP_MATCHED',len(checks),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--reference',type=Path);a=p.parse_args();main(a.output,a.reference)
