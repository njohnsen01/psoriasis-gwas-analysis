from pathlib import Path
import hashlib,importlib.metadata as im,json,os,sys,traceback
import numpy as np,pandas as pd
from google.cloud import bigquery
N=Path(__file__).resolve().parents[1];B=N.parent;R=Path('/home/jupyter/psoriasis_reproducibility_20260910');X=B/'X_source_reanalysis_20260910';O=N/'cohort_audit';O.mkdir(exist_ok=True)
sys.path.insert(0,str(N/'scripts'));from reconcile_X_site_filters import sha,save
sources={}
def csv(path,**kw):
 path=Path(path);sources[str(path)]={'sha256':sha(path),'bytes':path.stat().st_size};return pd.read_csv(path,sep=chr(9),**kw)
def tally(d,stage):
 cases=int(d.psoriasis_case.sum());return {'stage':stage,'N':len(d),'cases':cases,'controls':len(d)-cases}
def summary(d,col):
 rows=[]
 for case,z in d.groupby('psoriasis_case',sort=True):
  a=pd.to_numeric(z[col],errors='raise').dropna().astype(float)
  rows.append({'phenotype':'case' if case==1 else 'control','N_total':len(z),'N_available':len(a),'N_missing':len(z)-len(a),'median':a.median(),'q1':a.quantile(.25),'q3':a.quantile(.75),'mean':a.mean(),'sd':a.std(),'zero':int(a.eq(0).sum())})
 return pd.DataFrame(rows)
def main():
 ph=csv(R/'replay_cohort_v2/psoriasis_pheno.tsv',dtype={'FID':str,'IID':str}).set_index('IID')
 cv=csv(R/'replay_cohort_v2/psoriasis_covariates.tsv',dtype={'FID':str,'IID':str}).set_index('IID');d=ph.join(cv.drop(columns='FID'),validate='one_to_one');assert len(d)==392155 and d.psoriasis_case.sum()==5805
 metrics=csv(X/'private_provider_qc/genomics_metrics_Dec142023_1859_02_tz0000.tsv',dtype={'research_id':str}).set_index('research_id');assert d.index.isin(metrics.index).all();d=d.join(metrics,validate='one_to_one')
 coverage=csv(X/'private_provider/coverage_failures.tsv',dtype=str);flags=csv(X/'private_provider_qc/flagged_samples.tsv',dtype=str)
 d['SEX']=d.sex.map({1:1,0:2});assert d.SEX.notna().all()
 d['coverage_ok']=~d.index.isin(coverage.research_id);d['canonical']=d.dragen_sex_ploidy.isin(['XX','XY'])
 d['sex_concordant']=((d.SEX==1)&(d.dragen_sex_ploidy=='XY'))|((d.SEX==2)&(d.dragen_sex_ploidy=='XX'))
 d['provider_flag']=d.index.isin(flags.s);d['strict']=d.coverage_ok&d.sex_concordant&~d.provider_flag
 exclusions=[];flow=[tally(d,'original_analysis_cohort')];remaining=d.copy()
 for name,bad in [('provider_coverage_failure',~d.coverage_ok),('noncanonical_or_unavailable_ploidy',~d.canonical),('canonical_sex_discordance',d.canonical&~d.sex_concordant),('provider_QC_flag',d.provider_flag)]:
  removed=remaining.loc[bad.reindex(remaining.index)];exclusions.append(tally(removed,name));remaining=remaining.loc[~bad.reindex(remaining.index)];flow.append(tally(remaining,'after_'+name))
 keep=csv(N/'inputs/phenotype.tsv',dtype={'FID':str,'IID':str}).set_index('IID');assert set(remaining.index)==set(keep.index)==set(d.index[d.strict])
 assert len(remaining)==386069 and remaining.psoriasis_case.sum()==5705
 assert sum(x['N'] for x in exclusions)+len(remaining)==len(d)
 pd.DataFrame(flow).to_csv(O/'sample_flow.private.tsv',sep=chr(9),index=False);pd.DataFrame(exclusions).to_csv(O/'sequential_exclusions.private.tsv',sep=chr(9),index=False)
 oldpath=Path('/home/jupyter/psoriasis_phase2_reproducibility/secure_participant_level/ancestry_preds.tsv');newpath=B/'phenotype_validation_20260910/private_full_ancestry.tsv'
 old=csv(oldpath,usecols=['research_id','ancestry_pred'],dtype={'research_id':str}).set_index('research_id');new=csv(newpath,usecols=['research_id','ancestry_pred'],dtype={'research_id':str}).set_index('research_id');assert old.index.is_unique and new.index.is_unique
 d['legacy_ancestry']=old.ancestry_pred.reindex(d.index).str.upper().fillna('UNASSIGNED');d['complete_provider_ancestry']=new.ancestry_pred.reindex(d.index).str.upper();assert d.complete_provider_ancestry.notna().all()
 ar=[]
 for pop,rows,col in [('original_cohort_legacy_file',d,'legacy_ancestry'),('original_cohort_complete_file',d,'complete_provider_ancestry'),('strict_cohort_complete_file',d.loc[d.strict],'complete_provider_ancestry')]:
  for group,z in rows.groupby(col):ar.append({'population':pop,'group':group,**tally(z,'before_variant_specific_QC')})
 pd.DataFrame(ar).to_csv(O/'ancestry_group_counts.private.tsv',sep=chr(9),index=False)
 old_assigned=d.legacy_ancestry.ne('UNASSIGNED');disagree=old_assigned&d.legacy_ancestry.ne(d.complete_provider_ancestry)
 ancestry={'legacy_file_participants':len(old),'complete_file_participants':len(new),'legacy_assigned_original_cohort':int(old_assigned.sum()),'legacy_unassigned_original_cohort':int((~old_assigned).sum()),'complete_assigned_original_cohort':len(d),'labels_differ_among_previously_assigned':int(disagree.sum()),'legacy_source_version':'not inferred from file size; retained-file SHA identifies exact input','complete_source':'echo_v4_r2.ancestry_preds.tsv; provider generation 1759636181271952','historical_followup_and_corrected_source_groups_are_separate':True}
 assert ancestry['labels_differ_among_previously_assigned']==986;save(O/'ancestry_reconciliation.json',ancestry)
 clinicalpath=B/'phenotype_validation_20260910/private_phenotype_audit.pkl';clinical=pd.read_pickle(clinicalpath);clinical.person_id=clinical.person_id.astype(str);assert clinical.person_id.is_unique
 sources[str(clinicalpath)]={'sha256':sha(clinicalpath),'bytes':clinicalpath.stat().st_size,'source_SQL_sha256':sha(B/'phenotype_validation_20260910/phenotype_audit.sql')}
 d=d.join(clinical.set_index('person_id')[['has_ehr_data','visit_dates','visit_span_days']],validate='one_to_one');assert d.has_ehr_data.notna().all()
 query=(B/'ehr_observation_span.sql').read_text();assert os.environ['WORKSPACE_CDR']=='wb-silky-artichoke-2408.C2024Q3R8';(O/'ehr_observation_span.sql').write_text(query)
 client=bigquery.Client(project=os.environ['GOOGLE_PROJECT']);job=client.query(query,job_config=bigquery.QueryJobConfig(use_query_cache=True,maximum_bytes_billed=7000000000));ehr=job.to_dataframe();ehr.person_id=ehr.person_id.astype(str);assert ehr.person_id.is_unique
 save(O/'ehr_query_receipt.private.json',{'job_id':job.job_id,'location':job.location,'cache_hit':job.cache_hit,'bytes_processed':job.total_bytes_processed,'query_sha256':sha(O/'ehr_observation_span.sql'),'CDR':os.environ['WORKSPACE_CDR']})
 d=d.join(ehr.set_index('person_id'),validate='one_to_one');d['ehr_span_years']=(pd.to_datetime(d.last_ehr_date)-pd.to_datetime(d.first_ehr_date)).dt.days/365.25;d.loc[d.has_ehr_data.ne(1),'ehr_span_years']=np.nan
 for label,frame in [('original',d),('strict',d.loc[d.strict])]:
  summary(frame,'age').to_csv(O/(label+'_age_summary.private.tsv'),sep=chr(9),index=False);summary(frame,'ehr_span_years').to_csv(O/(label+'_EHR_summary.private.tsv'),sep=chr(9),index=False)
  flags=[]
  for case,z in frame.groupby('psoriasis_case'):
   flags.append({'phenotype':'case' if case else 'control','N':len(z),'has_EHR_flag':int(z.has_ehr_data.eq(1).sum()),'without_EHR_flag':int(z.has_ehr_data.ne(1).sum()),'EHR_and_two_visit_dates_and_365_day_span':int((z.has_ehr_data.eq(1)&z.visit_dates.ge(2)&z.visit_span_days.ge(365)).sum())})
  pd.DataFrame(flags).to_csv(O/(label+'_observation_flags.private.tsv'),sep=chr(9),index=False)
 d.to_pickle(O/'cohort_clinical.private.pkl')
 save(O/'input_receipts.private.json',sources)
 save(O/'COMPLETE.json',{'status':'complete','original':tally(d,'original'),'strict':tally(d.loc[d.strict],'strict'),'ancestry_reconciliation':ancestry,'exclusions_are_sequential_nonoverlapping':True,'EHR_span_definition':'earliest to latest recorded condition/procedure/drug/clinical-visit date; not continuous follow-up; missing when EHR flag absent','eligible_control_observation_definition':'EHR flag, at least two distinct clinical visit dates and visit span at least 365 days; diagnostic only, not silently imposed on primary','release_status':'private detailed counts; disclosure review required before public tables','code_sha256':sha(Path(__file__)),'versions':{p:im.version(p) for p in ['pandas','numpy','google-cloud-bigquery']}})
 print('CORRECTED_COHORT_FLOW_AND_OBSERVATION_AUDIT_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
