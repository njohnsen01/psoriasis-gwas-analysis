import os
from pathlib import Path
import json,hashlib,gzip,re
import numpy as np,pandas as pd
from scipy.stats import chi2,t
R=Path(__file__).resolve().parent
P=Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}');F=P/'final_aggregate'
plan=json.loads((P/'status/step2_plan.json').read_text());S=Path(plan['prediction_list']).parent
report={'scope':'original primary output, not source-corrected reanalysis'}
# Validate all LOCO rows without exporting IDs or predictions.
loco=Path(plan['loco_file']);h=hashlib.sha256();rows=[];header=None
with loco.open('rb') as f:
 for i,line in enumerate(f):
  h.update(line);v=line.decode().split()
  if i==0:header=v[1:];assert v[0]=='FID_IID';continue
  a=np.array(v[1:],dtype=float);assert len(a)==len(header) and np.isfinite(a).all()
  rows.append(v[0])
  if v[0]=='23':xhash=hashlib.sha256(line).hexdigest();xn=len(a);xvar=float(a.var())
assert rows==[str(i) for i in range(1,24)] and len(set(header))==len(header)
ph=pd.read_csv(plan['phenotype_file'],sep=r'\s+',dtype=str);ids=(ph.iloc[:,0]+'_'+ph.iloc[:,1]).tolist();assert set(ids)==set(header)
report['prediction']={'rows':rows,'n':xn,'row23_all_finite':True,'row23_nonconstant':xvar>0,'sample_ids_match_phenotype':True,'loco_sha256':h.hexdigest(),'row23_sha256':xhash,'prediction_list_sha256':hashlib.sha256(Path(plan['prediction_list']).read_bytes()).hexdigest()}
# Reconcile selected reported variants and source output conversions.
B=Path('/home/jupyter/psoriasis_corrected_downstream_2026-07-31/aggregate_outputs')
bio=pd.read_csv(B/'corrected_25_lead_biomarker_results_full_cohort.tsv',sep='\t')
selected=set(bio.variant_id)|{'chr6:31299200:G:A','chr6:163803171:CTTATT:C'}
thresholds=[5e-8,3e-8,5e-9,1e-9,0.05/13255796]
counts={str(q):{'all':0,'autosome':0,'X':0,'MHC':0,'non_MHC':0} for q in thresholds}
errors={k:0.0 for k in ['OR_relative','CI_lower_relative','CI_upper_relative','P_relative','CHISQ_LOG10P_abs']}
nt=0;extra={};found=[];invalid=0
for df in pd.read_csv(F/'psoriasis_corrected_regenie_all_results.tsv.gz',sep='\t',chunksize=250000,dtype={'CHROM':str}):
 nt+=len(df);v=df[['BETA','SE','LOG10P','CHISQ','OR','CI_LOWER_OR','CI_UPPER_OR','P']].to_numpy(float);invalid+=int((~np.isfinite(v).all(axis=1)|(df.SE<=0)|(df.P<=0)|(df.P>1)).sum())
 for k,o,e in [('OR_relative',df.OR,np.exp(df.BETA)),('CI_lower_relative',df.CI_LOWER_OR,np.exp(df.BETA-1.96*df.SE)),('CI_upper_relative',df.CI_UPPER_OR,np.exp(df.BETA+1.96*df.SE)),('P_relative',df.P,10.0**(-df.LOG10P))]:errors[k]=max(errors[k],float(np.max(np.abs(o/e-1))))
 errors['CHISQ_LOG10P_abs']=max(errors['CHISQ_LOG10P_abs'],float(np.max(np.abs(-chi2.logsf(df.CHISQ,1)/np.log(10)-df.LOG10P))))
 for k,n in df.EXTRA.fillna('NA').value_counts().items():extra[str(k)]=extra.get(str(k),0)+int(n)
 isx=df.CHROM.isin(['23','X','chrX']);mhc=(df.CHROM.isin(['6','chr6']))&df.POS.between(25000000,34000000)
 for q in thresholds:
  hit=df.LOG10P>-np.log10(q);d=counts[str(q)]
  for label,mask in [('all',np.ones(len(df),bool)),('X',isx),('autosome',~isx),('MHC',mhc),('non_MHC',~mhc)]:d[label]+=int((hit&mask).sum())
 found.append(df[df.ID.isin(selected)])
 if nt%1000000==0:print('ROWS',nt,flush=True)
leads=pd.concat(found);assert len(leads)==len(selected) and leads.ID.is_unique
leadmap=leads.set_index('ID');raw=[];logs=[]
for f in sorted((P/'aggregate_results').glob('chr*/*.regenie')):
 d=pd.read_csv(f,sep=r'\s+',dtype={'CHROM':str});q=d[d.ID.isin(selected)]
 if len(q):q=q.assign(source_file=str(f.relative_to(P)));raw.append(q)
 log=Path(str(f).replace('_psoriasis_case.regenie','.log'));assert log.exists();s=log.read_text();fails=re.findall(r'Number of failed tests\s*:\s*\((\d+)/(\d+)\)',s,re.I)
 logs.append({'file':str(log.relative_to(P)),'sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'failed':int(fails[-1]) if fails else None,'completed':'End time:' in s})
raw=pd.concat(raw);assert len(raw)==len(selected) and raw.ID.is_unique
merged=raw.merge(leads,on='ID',suffixes=('_raw','_merged'),validate='one_to_one')
match={c:bool(np.allclose(merged[c+'_raw'],merged[c+'_merged'],rtol=1e-12,atol=1e-12)) for c in ['BETA','SE','LOG10P','CHISQ','N']}
match['alleles']=bool((merged.ALLELE0==merged.REF).all() and (merged.ALLELE1==merged.ALT).all());assert all(match.values())
raw.to_csv(R/'primary_leads_source_rows.tsv',sep='\t',index=False)
leads.to_csv(R/'primary_leads_derived_rows.tsv',sep='\t',index=False)
report.update({'n_tests':nt,'invalid_numeric_rows':invalid,'conversion_max_errors':errors,'EXTRA_counts':extra,'threshold_counts':counts,'selected_leads_n':len(selected),'raw_merge_checks':match,'raw_chunk_logs_n':len(logs),'logs_missing_failed_counter':sum(x['failed'] is None for x in logs),'failed_tests_sum':sum(x['failed'] or 0 for x in logs),'logs_without_end_marker':sum(not x['completed'] for x in logs)})
# Check both original biomarker model sets using their own estimates, SEs and df.
biocheck=[]
for suffix in ['full_cohort','cases_only']:
 d=pd.read_csv(B/('corrected_25_lead_biomarker_results_'+suffix+'.tsv'),sep='\t');m=d.merge(leads[['ID','ALT','BETA','P']],left_on='variant_id',right_on='ID',validate='many_to_one')
 assert (m.effect_allele==m.ALT).all() and np.allclose(m.corrected_gwas_beta,m.BETA) and np.allclose(m.corrected_gwas_p,m.P)
 biocheck.append({'analysis':suffix,'n':len(d),'finite':bool(np.isfinite(d[['beta','se','p_value','df_resid']]).all().all()),'positive_se_df':bool((d.se>0).all() and (d.df_resid>0).all()),'effect_alleles_match':True,'t_p_error':float(np.max(abs(2*t.sf(abs(d.beta/d.se),d.df_resid)-d.p_value))),'ci_1_96_error':float(max(abs(d.ci_lower-(d.beta-1.96*d.se)).max(),abs(d.ci_upper-(d.beta+1.96*d.se)).max()))})
report['biomarker_checks']=biocheck
(R/'primary_log_checks.json').write_text(json.dumps(logs,indent=2));(R/'audit_summary.json').write_text(json.dumps(report,indent=2));print('AUDIT_COMPLETE',json.dumps(report),flush=True)
