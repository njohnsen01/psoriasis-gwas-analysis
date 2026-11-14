import os
from pathlib import Path
import json,hashlib,re,numpy as np,pandas as pd
mroot=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09')
ma=mroot/'model_threshold_audit_20260910';ma.mkdir(exist_ok=True)
primary=Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}')
print('PRIMARY_FILES',[(p.name,p.is_dir()) for p in primary.iterdir()])
print('STEP1_FILES',[(str(p.relative_to('/home/jupyter')),p.stat().st_size) for p in Path('/home/jupyter/psoriasis_regenie_loco_indexed_range_step1').rglob('*') if p.is_file() and ('pred' in p.name or 'plan' in p.name or 'valid' in p.name)])
print('STEP2_CODE', '\n'.join(Path('/home/jupyter/psoriasis_step2_parallel_runner_v2.py').read_text().splitlines()[239:267]))
print('AUDIT_READY')

step1=Path(f'/home/jupyter/psoriasis_regenie_loco_indexed_range_step1/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}')
predlist=step1/'regenie_step1/psoriasis_primary_loco_pruned_pred.list'
print('PRED_LIST',predlist.read_text())
print('STEP1_VALIDATION',(step1/'status/regenie_step1_pruned_validation.json').read_text())
print('STEP2_STATUS_FILES',[(p.name,p.stat().st_size) for p in (primary/'status').iterdir() if p.is_file()][:25])
print('DOWNSTREAM',[(str(p.relative_to('/home/jupyter')),p.stat().st_size) for base in ['psoriasis_corrected_downstream_2026-07-31','psoriasis_corrected_downstream_2026-08-01','psoriasis_corrected_downstream_2026-08-07'] for p in Path('/home/jupyter',base).rglob('*') if p.is_file() and (p.suffix=='.py' or p.suffix=='.json' or ('biomarker' in p.name and 'private' not in str(p)))][:45])
print('AUDIT_INVENTORY_DONE')

pplan=json.loads((primary/'status/step2_plan.json').read_text());print('PLAN_KEYS',list(pplan))
print('PLAN_RELEVANT',{k:v for k,v in pplan.items() if k!='tasks' and any(s in k for s in ['prediction','phenotype','covariate','version','sha','attempt'])})
print('TASK_KEYS',list(pplan['tasks'][0]) if isinstance(pplan.get('tasks'),list) else type(pplan.get('tasks')).__name__)
print('SUMMARY',(primary/'status/step2_summary.json').read_text())
fa=primary/'final_aggregate'
for fn in ['psoriasis_corrected_regenie_all_results.tsv.gz','psoriasis_corrected_non_mhc_lead_groups.tsv','psoriasis_corrected_lead_centered_100kb_groups.tsv']:
 import gzip
 p=fa/fn
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:print('HEADER',fn,next(f).strip())
bio=Path('/home/jupyter/psoriasis_corrected_downstream_2026-07-31/aggregate_outputs')
print('BIO_SUMMARY',(bio/'corrected_25_lead_biomarker_analysis_summary.json').read_text())
print('BIO_CODE_CANDIDATES',(bio/'corrected_downstream_biomarker_code_candidates.tsv').read_text())
print('READY_NEXT')

print('CHUNK_SAMPLE',[(str(p.relative_to(primary)),p.stat().st_size) for p in (primary/'aggregate_results/chrX').iterdir()][:8])
print('PRIMARY_LOGS',len(list(primary.rglob('*.log'))))
print('PRIVATE_BIO_NAMES',[p.name for p in (bio.parent/'private_participant_level').iterdir()])
print('BIO_VALIDATION',(bio/'corrected_25_lead_biomarker_model_validation.tsv').read_text())
print('BIO_RESULT_HEADER',pd.read_csv(bio/'corrected_25_lead_biomarker_results_full_cohort.tsv',sep='\t',nrows=0).columns.tolist())
print('PRIMARY_GROUPS_N',len(pd.read_csv(fa/'psoriasis_corrected_non_mhc_lead_groups.tsv',sep='\t')))
print('STEP1_PLAN_KEYS',list(json.loads((step1/'status/regenie_step1_plan_pruned.json').read_text())))
print('INVENTORY_2_DONE')

audit_code='import os\nfrom pathlib import Path\nimport json,hashlib,gzip,re\nimport numpy as np,pandas as pd\nfrom scipy.stats import chi2,t\nR=Path(__file__).resolve().parent\nP=Path(f\'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}\');F=P/\'final_aggregate\'\nplan=json.loads((P/\'status/step2_plan.json\').read_text());S=Path(plan[\'prediction_list\']).parent\nreport={\'scope\':\'original primary output, not source-corrected reanalysis\'}\n# Validate all LOCO rows without exporting IDs or predictions.\nloco=Path(plan[\'loco_file\']);h=hashlib.sha256();rows=[];header=None\nwith loco.open(\'rb\') as f:\n for i,line in enumerate(f):\n  h.update(line);v=line.decode().split()\n  if i==0:header=v[1:];assert v[0]==\'FID_IID\';continue\n  a=np.array(v[1:],dtype=float);assert len(a)==len(header) and np.isfinite(a).all()\n  rows.append(v[0])\n  if v[0]==\'23\':xhash=hashlib.sha256(line).hexdigest();xn=len(a);xvar=float(a.var())\nassert rows==[str(i) for i in range(1,24)] and len(set(header))==len(header)\nph=pd.read_csv(plan[\'phenotype_file\'],sep=r\'\\s+\',dtype=str);ids=(ph.iloc[:,0]+\'_\'+ph.iloc[:,1]).tolist();assert set(ids)==set(header)\nreport[\'prediction\']={\'rows\':rows,\'n\':xn,\'row23_all_finite\':True,\'row23_nonconstant\':xvar>0,\'sample_ids_match_phenotype\':True,\'loco_sha256\':h.hexdigest(),\'row23_sha256\':xhash,\'prediction_list_sha256\':hashlib.sha256(Path(plan[\'prediction_list\']).read_bytes()).hexdigest()}\n# Reconcile selected reported variants and source output conversions.\nB=Path(\'/home/jupyter/psoriasis_corrected_downstream_2026-07-31/aggregate_outputs\')\nbio=pd.read_csv(B/\'corrected_25_lead_biomarker_results_full_cohort.tsv\',sep=\'\\t\')\nselected=set(bio.variant_id)|{\'chr6:31299200:G:A\',\'chr6:163803171:CTTATT:C\'}\nthresholds=[5e-8,3e-8,5e-9,1e-9,0.05/13255796]\ncounts={str(q):{\'all\':0,\'autosome\':0,\'X\':0,\'MHC\':0,\'non_MHC\':0} for q in thresholds}\nerrors={k:0.0 for k in [\'OR_relative\',\'CI_lower_relative\',\'CI_upper_relative\',\'P_relative\',\'CHISQ_LOG10P_abs\']}\nnt=0;extra={};found=[];invalid=0\nfor df in pd.read_csv(F/\'psoriasis_corrected_regenie_all_results.tsv.gz\',sep=\'\\t\',chunksize=250000,dtype={\'CHROM\':str}):\n nt+=len(df);v=df[[\'BETA\',\'SE\',\'LOG10P\',\'CHISQ\',\'OR\',\'CI_LOWER_OR\',\'CI_UPPER_OR\',\'P\']].to_numpy(float);invalid+=int((~np.isfinite(v).all(axis=1)|(df.SE<=0)|(df.P<=0)|(df.P>1)).sum())\n for k,o,e in [(\'OR_relative\',df.OR,np.exp(df.BETA)),(\'CI_lower_relative\',df.CI_LOWER_OR,np.exp(df.BETA-1.96*df.SE)),(\'CI_upper_relative\',df.CI_UPPER_OR,np.exp(df.BETA+1.96*df.SE)),(\'P_relative\',df.P,10.0**(-df.LOG10P))]:errors[k]=max(errors[k],float(np.max(np.abs(o/e-1))))\n errors[\'CHISQ_LOG10P_abs\']=max(errors[\'CHISQ_LOG10P_abs\'],float(np.max(np.abs(-chi2.logsf(df.CHISQ,1)/np.log(10)-df.LOG10P))))\n for k,n in df.EXTRA.fillna(\'NA\').value_counts().items():extra[str(k)]=extra.get(str(k),0)+int(n)\n isx=df.CHROM.isin([\'23\',\'X\',\'chrX\']);mhc=(df.CHROM.isin([\'6\',\'chr6\']))&df.POS.between(25000000,34000000)\n for q in thresholds:\n  hit=df.LOG10P>-np.log10(q);d=counts[str(q)]\n  for label,mask in [(\'all\',np.ones(len(df),bool)),(\'X\',isx),(\'autosome\',~isx),(\'MHC\',mhc),(\'non_MHC\',~mhc)]:d[label]+=int((hit&mask).sum())\n found.append(df[df.ID.isin(selected)])\n if nt%1000000==0:print(\'ROWS\',nt,flush=True)\nleads=pd.concat(found);assert len(leads)==len(selected) and leads.ID.is_unique\nleadmap=leads.set_index(\'ID\');raw=[];logs=[]\nfor f in sorted((P/\'aggregate_results\').glob(\'chr*/*.regenie\')):\n d=pd.read_csv(f,sep=r\'\\s+\',dtype={\'CHROM\':str});q=d[d.ID.isin(selected)]\n if len(q):q=q.assign(source_file=str(f.relative_to(P)));raw.append(q)\n log=Path(str(f).replace(\'_psoriasis_case.regenie\',\'.log\'));assert log.exists();s=log.read_text();fails=re.findall(r\'Number of failed tests\\s*:\\s*(\\d+)\',s,re.I)\n logs.append({\'file\':str(log.relative_to(P)),\'sha256\':hashlib.sha256(log.read_bytes()).hexdigest(),\'failed\':int(fails[-1]) if fails else None,\'completed\':\'End time:\' in s})\nraw=pd.concat(raw);assert len(raw)==len(selected) and raw.ID.is_unique\nmerged=raw.merge(leads,on=\'ID\',suffixes=(\'_raw\',\'_merged\'),validate=\'one_to_one\')\nmatch={c:bool(np.allclose(merged[c+\'_raw\'],merged[c+\'_merged\'],rtol=1e-12,atol=1e-12)) for c in [\'BETA\',\'SE\',\'LOG10P\',\'CHISQ\',\'N\']}\nmatch[\'alleles\']=bool((merged.ALLELE0==merged.REF).all() and (merged.ALLELE1==merged.ALT).all());assert all(match.values())\nraw.to_csv(R/\'primary_leads_source_rows.tsv\',sep=\'\\t\',index=False)\nleads.to_csv(R/\'primary_leads_derived_rows.tsv\',sep=\'\\t\',index=False)\nreport.update({\'n_tests\':nt,\'invalid_numeric_rows\':invalid,\'conversion_max_errors\':errors,\'EXTRA_counts\':extra,\'threshold_counts\':counts,\'selected_leads_n\':len(selected),\'raw_merge_checks\':match,\'raw_chunk_logs_n\':len(logs),\'logs_missing_failed_counter\':sum(x[\'failed\'] is None for x in logs),\'failed_tests_sum\':sum(x[\'failed\'] or 0 for x in logs),\'logs_without_end_marker\':sum(not x[\'completed\'] for x in logs)})\n# Check both original biomarker model sets using their own estimates, SEs and df.\nbiocheck=[]\nfor suffix in [\'full_cohort\',\'cases_only\']:\n d=pd.read_csv(B/(\'corrected_25_lead_biomarker_results_\'+suffix+\'.tsv\'),sep=\'\\t\');m=d.merge(leads[[\'ID\',\'ALT\',\'BETA\',\'P\']],left_on=\'variant_id\',right_on=\'ID\',validate=\'many_to_one\')\n assert (m.effect_allele==m.ALT).all() and np.allclose(m.corrected_gwas_beta,m.BETA) and np.allclose(m.corrected_gwas_p,m.P)\n biocheck.append({\'analysis\':suffix,\'n\':len(d),\'finite\':bool(np.isfinite(d[[\'beta\',\'se\',\'p_value\',\'df_resid\']]).all().all()),\'positive_se_df\':bool((d.se>0).all() and (d.df_resid>0).all()),\'effect_alleles_match\':True,\'t_p_error\':float(np.max(abs(2*t.sf(abs(d.beta/d.se),d.df_resid)-d.p_value))),\'ci_1_96_error\':float(max(abs(d.ci_lower-(d.beta-1.96*d.se)).max(),abs(d.ci_upper-(d.beta+1.96*d.se)).max()))})\nreport[\'biomarker_checks\']=biocheck\n(R/\'primary_log_checks.json\').write_text(json.dumps(logs,indent=2));(R/\'audit_summary.json\').write_text(json.dumps(report,indent=2));print(\'AUDIT_COMPLETE\',json.dumps(report),flush=True)\n'
(ma/'audit_primary.py').write_text(audit_code)
import subprocess
auditproc=subprocess.Popen(['python',str(ma/'audit_primary.py')],stdout=(ma/'audit_primary.console.log').open('w'),stderr=subprocess.STDOUT)
print('AUDIT_PID',auditproc.pid)
print('CURRENT_X_LOG','\n'.join((primary/'aggregate_results/chrX/chrX_pso_chunk_0020.log').read_text().splitlines()[-18:]))

print('AUDIT_STATUS',auditproc.poll())
print((ma/'audit_primary.console.log').read_text()[-9000:])
print('MODEL_FILES', '\n'.join(str(p) for p in mroot.rglob('*') if p.is_file() and any(w in p.name.lower() for w in ['retained','biomarker','dosage']) and p.suffix in ['.py','.pgen','.json','.tsv']))
print('FOLLOWUP_LOGS','\n'.join(str(p) for p in Path('/home/jupyter/psoriasis_corrected_downstream_2026-08-07/x_followup_consolidated').rglob('*.log'))[:4000])

exec(In[222].replace("'.pgen']))[:7000])","'.pgen'])[:7000])"))

print('FOLLOWUP_ROOTS','\n'.join(str(p) for p in Path('/home/jupyter').iterdir() if any(w in p.name.lower() for w in ['phase3','phase4','sex','biomarker'])))
print('BIO_LOG',(Path('/home/jupyter/psoriasis_corrected_downstream_2026-07-31/private_participant_level/corrected_missing_lead_exports/chrX_chunk_0022_corrected_leads.log').read_text())[:4500])
bdpath=Path('/home/jupyter/psoriasis_corrected_downstream_2026-07-31/private_participant_level/corrected_25_lead_dosages_effect_allele.tsv.gz')
print('BIO_HEADER',pd.read_csv(bdpath,sep='\t',nrows=0).columns.tolist())
print('RETAINED_FILES','\n'.join(str(p) for p in (mroot/'Xq12_conditional/private_genotypes').iterdir()))
print('FOURLEAD_FILES','\n'.join(str(p) for p in mroot.rglob('*.pgen'))[:6500])
print('SEX_MANIFESTS','\n'.join(str(p) for p in Path('/home/jupyter/psoriasis_corrected_downstream_2026-08-07/x_followup_consolidated').iterdir()))

print('PHASE4_FILES','\n'.join(str(p) for p in Path('/home/jupyter/phase4_chrX_sex_aware_2026-07-23').rglob('*') if p.is_file() and p.suffix in ['.py','.log','.json','.tsv'])[:9500])
print('DOSAGE_MANIFESTS','\n'.join(str(p) for p in bio.iterdir() if any(w in p.name for w in ['dosage','manifest','align'])))
print('RETAINED_STAGE','\n'.join(str(p) for p in Path('/home/jupyter/psoriasis_corrected_downstream_2026-07-31/private_participant_level/targeted_chunk_staging').rglob('*.pgen')))

ph4=Path('/home/jupyter/phase4_chrX_sex_aware_2026-07-23')
for name in ['phase4_model_specification.json','phase4_genotype_source.json']:
 print(name,(ph4/'aggregate_outputs'/name).read_text())
print('MALE_COMMAND',(ph4/'private_participant_level/plink2_chrX_sex_stratified/chrX_three_leads_male.log').read_text()[:3300])
for name in ['corrected_25_lead_dosage_availability.tsv','corrected_25_lead_dosage_qc.tsv','corrected_8_missing_lead_extraction_manifest.tsv']:
 d=pd.read_csv(bio/name,sep='\t');print(name,d.columns.tolist());print(d[d.astype(str).apply(lambda z:z.str.contains('chrX')).any(axis=1)].to_string(index=False))
print('PRIMARY_RAW_X',pd.read_csv(ma/'primary_leads_source_rows.tsv',sep='\t').query('CHROM == 23').to_json(orient='records'))

bq=pd.read_csv(bio/'corrected_25_lead_dosage_qc.tsv',sep='\t');bx=bq[bq.variant_id.str.startswith('chrX:')]
bd=pd.read_csv(bdpath,sep='\t',dtype={'research_id':str}).set_index('research_id')
cv=pd.read_csv(pplan['covariate_file'],sep=r'\s+',dtype={'FID':str,'IID':str});cv=cv.set_index('IID');sx=cv.sex.reindex(bd.index);assert sx.notna().all()
oldraw=Path('/home/jupyter/psoriasis_phase2_reproducibility/secure_participant_level/lead_variant_dosages/chrX_lead_dosages.raw')
raws=[pd.read_csv(oldraw,sep=r'\s+',dtype={'IID':str}).set_index('IID'),pd.read_csv(Path('/home/jupyter/psoriasis_corrected_downstream_2026-07-31/private_participant_level/corrected_missing_lead_exports/chrX_chunk_0022_corrected_leads.raw'),sep=r'\s+',dtype={'IID':str}).set_index('IID')]
dchecks=[]
for z in bx.itertuples():
 ds=bd[z.aligned_column];raw=next(r for r in raws if z.variant_id+'_'+z.ref in r.columns);orig=raw[z.variant_id+'_'+z.ref].reindex(bd.index);v=2-orig;dchecks.append({'variant':z.variant_id,'effect_allele':z.effect_allele,'column':z.aligned_column,'raw_REF_to_ALT_2_minus_dosage_exact':bool(np.allclose(ds,v,equal_nan=True)),'male_stored_codes':sorted(ds[sx==1].dropna().unique().tolist()),'female_stored_codes':sorted(ds[sx==0].dropna().unique().tolist())})
print('DOSAGE_VECTOR_CHECKS',json.dumps(dchecks))
print('PLINK_RAW','\n'.join(p.read_text() for p in (ph4/'private_participant_level/plink2_chrX_sex_stratified').glob('*.glm.logistic.hybrid'))[:6000])
print('CONSOLIDATED_FILES','\n'.join(str(p) for p in Path('/home/jupyter/psoriasis_corrected_downstream_2026-08-07').rglob('*') if p.is_file())[:6500])
(ma/'biomarker_dosage_checks.json').write_text(json.dumps(dchecks,indent=2))

exec(In[228].replace("'.gz']))[:2000])","'.gz'])[:2000])"))

xf=Path(f'/home/jupyter/psoriasis_corrected_downstream_2026-08-01/x_chromosome_followup_aggregate_outputs/attempt_{os.environ["PSO_X_FOLLOWUP_RUN_ID"]}');print('X_CURRENT',[(p.name,p.is_file(),p.is_symlink()) for p in xf.iterdir()])
xrec=next(x for x in pplan['chunk_records'] if str(x['chromosome'])=='X');print('X_CHUNK_KEYS',list(xrec));print('X_PGEN_URIS',[u for u in xrec.get('pgen_uris',[]) if any(t in u for t in ['0007','0011','0022'])])
print('LOCAL_X_SCRIPT',[(p.name) for p in Path('/home/jupyter').glob('*.py') if any(t in p.name for t in ['followup','biomarker','sex'])])

from google.cloud import storage
mc=storage.Client();mb=mc.bucket(f'{os.environ["WORKSPACE_BUCKET"].removeprefix("gs://").rstrip("/")}')
print('ARCHIVED_X_FILES',[b.name for b in mc.list_blobs(mb,prefix='psoriasis_corrected_downstream_2026-08-01/x_chromosome_followup_aggregate_outputs/') if not b.name.endswith('/')][:30])
(ma/'private_primary_inputs').mkdir(exist_ok=True)
import concurrent.futures
mxpool=concurrent.futures.ThreadPoolExecutor(max_workers=3)
def get_primary_input(uri):
 name=uri.split('/',3)[3];dest=ma/'private_primary_inputs'/name.split('/')[-1];blob=mb.blob(name);blob.reload();blob.download_to_filename(str(dest));return {'name':dest.name,'bytes':dest.stat().st_size,'md5':blob.md5_hash}
uris=[u for k in ['pgen_uris','pvar_uris'] for u in xrec[k] if any(t in u for t in ['0007','0011'])]+[xrec['psam_uri']]
input_jobs=[mxpool.submit(get_primary_input,u) for u in uris]
print('INPUT_DOWNLOAD_JOBS',len(input_jobs))

print('PRIMARY_DOWNLOADS',[j.result() for j in input_jobs if j.done()]);print('PENDING',sum(not j.done() for j in input_jobs))
import pgenlib
primary_dosage_checks=[]
for chunk,subset in [('0007',['chrX:5998724:A:G']),('0011',['chrX:22917619:C:T']),('0022',['chrX:67736589:T:C','chrX:67842269:G:T'])]:
 pr=ma/'private_primary_inputs'/('chrX_pso_chunk_'+chunk) if chunk!='0022' else mroot/'Xq12_conditional/private_genotypes/chrX_pso_chunk_0022';psfile=ma/'private_primary_inputs/chrX_pso_step2.psam';ps=pd.read_csv(psfile,sep=r'\s+',dtype={'IID':str});pv=pd.read_csv(str(pr)+'.pvar.zst',sep='\t');gv=np.empty(len(ps),dtype=np.int8)
 with pgenlib.PgenReader(str(pr.with_suffix('.pgen')).encode()) as rr:
  for vid in subset:
   j=pv.index[pv.ID==vid].item();rr.read(j,gv);g=pd.Series(gv.copy(),index=ps.IID).reindex(bd.index).astype(float);g[g<0]=np.nan;col=bx.set_index('variant_id').loc[vid,'aligned_column'];primary_dosage_checks.append({'variant':vid,'exact_match_biomarker_including_missing':bool(np.allclose(g,bd[col],equal_nan=True)),'male_codes':sorted(g[sx==1].dropna().unique().tolist()),'female_codes':sorted(g[sx==0].dropna().unique().tolist()),'pgen_sha256':hashlib.sha256(pr.with_suffix('.pgen').read_bytes()).hexdigest()})
print('PRIMARY_VECTOR_MATCH',json.dumps(primary_dosage_checks));(ma/'primary_dosage_checks.json').write_text(json.dumps(primary_dosage_checks,indent=2))

d=pd.read_csv(ma/'primary_leads_derived_rows.tsv',sep='\t').to_dict('records');ar=[[x['ID'],x['REF'],x['ALT']]+[format(x[c],'.12g') for c in ['N','BETA','SE','LOG10P']] for x in sorted(d,key=lambda x:x['ID'])];print('PRIMARY27_CANONICAL_SHA',hashlib.sha256(json.dumps(ar,separators=(',',':')).encode()).hexdigest())
for suffix in ['full_cohort','cases_only']:
 d=pd.read_csv(bio/('corrected_25_lead_biomarker_results_'+suffix+'.tsv'),sep='\t').to_dict('records');ar=[[x['variant_id'],x['trait'],x['effect_allele']]+[format(x[c],'.12g') for c in ['n','beta','se','ci_lower','ci_upper','p_value']] for x in sorted(d,key=lambda x:(x['variant_id'],x['trait']))];print('BIO_CANONICAL_SHA',suffix,hashlib.sha256(json.dumps(ar,separators=(',',':')).encode()).hexdigest())
import statsmodels.api as sm
from threadpoolctl import threadpool_limits
ph=pd.read_csv(pplan['phenotype_file'],sep=r'\s+',dtype={'FID':str,'IID':str}).set_index('IID');sdata=cv.join(ph[['psoriasis_case']]).join(bd);sexcols=['age','age2']+['PC'+str(k) for k in range(1,11)]
sexrefit=[]
with threadpool_limits(limits=2):
 for sexval in [0,1]:
  for vid in ['chrX:22917619:C:T','chrX:5998724:A:G','chrX:67842269:G:T']:
   col=bx.set_index('variant_id').loc[vid,'aligned_column'];dd=sdata.loc[sdata.sex==sexval,sexcols+[col,'psoriasis_case']].dropna();xx=dd[sexcols].astype(float);xx=(xx-xx.mean())/xx.std(ddof=0);xx[col]=dd[col];xx=sm.add_constant(xx);fit=sm.GLM(dd.psoriasis_case,xx,family=sm.families.Binomial()).fit(maxiter=100,tol=1e-8);sexrefit.append({'variant':vid,'sex':'male' if sexval else 'female','beta':float(fit.params[col]),'se':float(fit.bse[col]),'p':float(fit.pvalues[col]),'n':len(dd),'converged':bool(fit.converged),'rank':int(np.linalg.matrix_rank(xx))})
print('SEX_REFIT',json.dumps(sexrefit));(ma/'sex_refit_checks.json').write_text(json.dumps(sexrefit,indent=2))

for suffix in ['full_cohort','cases_only']:
 d=pd.read_csv(bio/('corrected_25_lead_biomarker_results_'+suffix+'.tsv'),sep='\t').to_dict('records');ar=[[x['variant_id'],x['trait'],x['effect_allele']]+[format(x[c],'.8g') for c in ['n','beta','se','ci_lower','ci_upper','p_value']] for x in sorted(d,key=lambda x:(x['variant_id'],x['trait']))];print('BIO_CANONICAL_SHA_8G',suffix,hashlib.sha256(json.dumps(ar,separators=(',',':')).encode()).hexdigest(),ar[0])
print('PRIMARY_STEP1_AUTOSOMAL',json.loads((step1/'status/step1_plan_pruned.json').read_text()).get('chromosomes')) if (step1/'status/step1_plan_pruned.json').exists() else None