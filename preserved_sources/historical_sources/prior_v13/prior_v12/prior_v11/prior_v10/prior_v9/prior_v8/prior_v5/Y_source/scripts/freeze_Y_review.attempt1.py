from pathlib import Path
import sys,json,hashlib,re,datetime,importlib.metadata as im,shutil,traceback
import numpy as np,pandas as pd
from PIL import Image
Y=Path(__file__).resolve().parents[1];N=Y.parent;F=Y/'freeze_review';F.mkdir(exist_ok=True)
sys.path.insert(0,str(N/'scripts'))
from reconcile_X_site_filters import sha,save
from finalize_pooled_calibration import quantiles
sys.path.insert(0,str(N.parent/'X_source_reanalysis_20260910/python_dependencies'))
import pgenlib
def main():
 result=json.loads((Y/'results/COMPLETE.json').read_text());assert result['status']=='validated'
 assert 'Y_PIPELINE_STAGE_RECEIPTS_VERIFIED' in (Y/'pipeline_verify.console.log').read_text()
 q=pd.read_csv(Y/'projected/final_variant_QC.private.tsv.gz',sep=chr(9)).set_index('ID')
 qc=json.loads((Y/'projected/COMPLETE.json').read_text());assert int(q.PASS.sum())==qc['passed']
 for f,h in qc['inputs'].items():assert sha(Y/'projected'/f)==h
 assert sha(Y/'scripts/project_Y_source.py')==qc['script_sha256'] and sha(Y/'scripts/run_Y_models.py')==result['script_sha256']
 checks=[]
 for scale in ['02','01']:
  name='male'+scale;r=result['reports'][name];p=Y/'results'/(name+'_psoriasis_case.regenie')
  assert sha(p)==r['source_output_sha256'];log=(Y/'results'/(name+'.log')).read_text()
  assert sha(Y/'results'/(name+'.log'))==r['model_log_sha256']
  warnings=[l for l in log.splitlines() if re.search('warning|error|converg',l,re.I)];assert not warnings,warnings
  d=pd.read_csv(p,sep=r'\s+').set_index('ID');derived=pd.read_csv(Y/'results'/(name+'.derived.private.tsv.gz'),sep=chr(9)).set_index('ID');d=d.reindex(derived.index)
  assert np.isfinite(d[['BETA','SE','LOG10P']]).all().all() and d.SE.gt(0).all() and d.LOG10P.ge(0).all()
  for col in ['N','N_CASES','N_CONTROLS','BETA','SE','LOG10P']:assert d[col].equals(derived[col])
  assert d.N.eq(d.N_CASES+d.N_CONTROLS).all()
  lam,qq=quantiles(d.LOG10P.to_numpy());stored=pd.read_csv(Y/'results'/(name+'_QQ_full_ranks.private.tsv.gz'),sep=chr(9))
  pd.testing.assert_frame_equal(qq,stored,check_exact=False,rtol=1e-14,atol=1e-14);assert lam==r['lambda_GC']
  assert (np.diff(stored.full_output_rank)>0).all() and stored.full_output_rank.iloc[0]==1 and stored.full_output_rank.iloc[-1]==len(d)
  assert np.allclose(derived.OR_per_coded_unit,np.exp(d.BETA)) and np.allclose(derived.CI95_Wald_low,np.exp(d.BETA-1.959963984540054*d.SE)) and np.allclose(derived.CI95_Wald_high,np.exp(d.BETA+1.959963984540054*d.SE))
  variant_q=q.reindex(d.index);assert variant_q.PASS.all()
  checks.append({'model':name,'source_and_derived_rows_match':True,'warnings':warnings,'full_rank_QQ_verified':True,'lambda_GC':lam,'tests':len(d),'N_matches_nonmissing_genotypes_all':bool(d.N.eq(variant_q.case_called+variant_q.control_called).all()),'N_cases_matches_nonmissing_cases_all':bool(d.N_CASES.eq(variant_q.case_called).all()),'N_controls_matches_nonmissing_controls_all':bool(d.N_CONTROLS.eq(variant_q.control_called).all())})
  with Image.open(Y/'results'/('Chromosome_Y_'+name+'_600dpi.tiff')) as img:assert all(abs(v-600)<.1 for v in img.info['dpi'])
 environment={'python':sys.version,'python_executable':sys.executable,'packages':{k:im.version(k) for k in ['numpy','pandas','scipy','matplotlib','hail','google-cloud-storage','google-crc32c']},'binaries':{b:{'path':shutil.which(b),'sha256':sha(Path(shutil.which(b)))} for b in ['plink2','regenie']},'pgenlib_binary':{'path':pgenlib.__file__,'sha256':sha(Path(pgenlib.__file__))},'version_lines':{'plink2':(Y/'counts/cases.log').read_text().splitlines()[0],'regenie':[l.strip() for l in (Y/'results/male02.log').read_text().splitlines() if 'REGENIE v' in l]}}
 save(F/'ENVIRONMENT.private.json',environment)
 save(F/'QUANTITATIVE_AUDIT.private.json',{'status':'validated','checks':checks,'dosage_comparison':result['dosage_comparison'],'scope':'Technical integrity and coding consistency; not validation of all scientific assumptions','release':'HELD'})
 save(F/'VISUAL_QA.json',{'status':'reviewed','figures':['Chromosome_Y_male02','Chromosome_Y_male01'],'method':'Rendered both previews in the secure console and visually inspected via browser screenshots','findings':'Titles, axis labels, test counts and calibration labels are visible; no clipping or overlap observed. TIFF metadata is 600 dpi. Final manuscript placement and numbering remain pending.'})
 files=[]
 for p in sorted(Y.rglob('*')):
  if not p.is_file() or F in p.parents or '__pycache__' in p.parts:continue
  files.append({'path':str(p.relative_to(Y)),'bytes':p.stat().st_size,'sha256':sha(p)})
 for p in sorted(F.glob('*.json')):files.append({'path':str(p.relative_to(Y)),'bytes':p.stat().st_size,'sha256':sha(p)})
 save(F/'SNAPSHOT_INDEX.private.json',{'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'validated_review_snapshot','files':files,'parent_autosomal_X_snapshot_sha256':sha(N/'freeze_review_20260912/SNAPSHOT_INDEX.private.json'),'manuscript_integration':'pending checklist item 2','release':'HELD'})
 print('Y_FREEZE_AUDIT_COMPLETE',len(files),sha(F/'SNAPSHOT_INDEX.private.json'),flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(F/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
