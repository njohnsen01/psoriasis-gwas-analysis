from pathlib import Path
import json,hashlib,datetime,traceback,importlib.metadata
import numpy as np,pandas as pd
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911'); B=N.parent; I=N/'reconciliation_20260912'; O=I/'cohort'; R=Path('/home/jupyter/psoriasis_reproducibility_20260910'); Q=B/'phenotype_validation_20260910'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def write(name,d):d.to_csv(O/(name+'.private.tsv'),sep='	',index=False)
def read_tsv(p):return pd.read_csv(p,sep='	',dtype={'FID':str,'IID':str}).set_index('IID')
def tally(d):
 c=int(d.psoriasis_case.sum());return {'N':len(d),'cases':c,'controls':len(d)-c}
def summ(d,col):
 rows=[]
 for case,z in d.groupby('psoriasis_case',sort=True):
  a=pd.to_numeric(z[col],errors='raise').dropna().astype(float)
  rows.append({'phenotype':'case' if case else 'control','N_total':len(z),'N_available':len(a),'N_missing':len(z)-len(a),'median':a.median(),'q1':a.quantile(.25),'q3':a.quantile(.75),'mean':a.mean(),'sd':a.std(),'zero':int(a.eq(0).sum())})
 return pd.DataFrame(rows)
def main():
 checks={};sources={}
 for receipt in [O/'INPUTS.private.json',N/'cohort_audit/input_receipts.private.json']:
  for name,meta in json.loads(receipt.read_text()).items():
   expected=meta if isinstance(meta,str) else meta['sha256']; actual=sha(name); assert actual==expected, 'changed input '+name
   sources[name]={'sha256':actual,'bytes':Path(name).stat().st_size}
 checks['preparation_and_prior_audit_source_hashes']=True
 up=pd.read_pickle(O/'upstream_eligibility.private.pkl').set_index('person_id');up.index=up.index.astype(str);assert up.index.is_unique
 cl=pd.read_pickle(Q/'private_phenotype_audit.pkl').set_index('person_id');cl.index=cl.index.astype(str);assert cl.index.is_unique
 d=pd.read_pickle(N/'cohort_audit/cohort_clinical.private.pkl');d.index=d.index.astype(str);assert d.index.is_unique
 ph=read_tsv(N/'inputs/phenotype.tsv'); cv=read_tsv(N/'inputs/covariates.tsv'); strict=d.loc[d.strict].copy()
 assert set(ph.index)==set(strict.index)==set(cv.index)
 pd.testing.assert_series_equal(ph.psoriasis_case.sort_index(),strict.psoriasis_case.sort_index(),check_names=False)
 pd.testing.assert_frame_equal(cv.drop(columns='FID').sort_index(),strict[cv.columns.drop('FID')].sort_index(),check_names=False,check_dtype=False)
 checks['strict_participants_phenotype_and_covariates']=True
 # Sequential flow; phenotype is not assigned to participants outside the eligible clinical universe.
 flow=[{'stage':'cb_search_person_source','N_retained':len(up),'N_excluded_at_stage':0,'cases_retained':pd.NA,'controls_retained':pd.NA,'definition':'Source inventory; not all source participants have eligible WGS'}]
 alive=up.index; exclusions=[]
 stages=[('WGS_available',up.has_whole_genome_variant.eq(1),'has_whole_genome_variant = 1'),('person_record_available',up.has_person_record,'Matching CDR person record'),('age_eligible',up.age_at_cdr.between(18,125),'Provider age_at_cdr from 18 through 125 inclusive'),('no_recorded_death',~up.has_death_record,'No participant record in the CDR death table'),('eligible_recorded_sex',up.sex_at_birth_concept_id.isin([45880669,45878463]),'Recorded sex_at_birth_concept_id in the two source-specified concepts')]
 for label,mask,definition in stages:
  keep=mask.reindex(alive).fillna(False).to_numpy(dtype=bool); removed=alive[~keep];alive=alive[keep]
  exclusions.append(set(removed));flow.append({'stage':label,'N_retained':len(alive),'N_excluded_at_stage':len(removed),'cases_retained':pd.NA,'controls_retained':pd.NA,'definition':definition})
 assert set(alive)==set(cl.index)
 for label,bad,definition in [('clinical_condition_exclusion',~cl.case_gate & cl.any_excluded_condition,'Without original case-arm eligibility, any source-specified excluded condition'),('control_drug_exclusion',~cl.case_gate & cl.any_excluded_drug,'After condition exclusion, without case-arm eligibility, any source-specified excluded drug')]:
  reject=bad.reindex(alive).fillna(False).to_numpy(dtype=bool);removed=alive[reject];alive=alive[~reject];exclusions.append(set(removed))
  flow.append({'stage':label,'N_retained':len(alive),'N_excluded_at_stage':len(removed),'cases_retained':pd.NA,'controls_retained':pd.NA,'definition':definition})
 assert set(alive)==set(up.index[up.selected_by_original_SQL.fillna(False)])
 clinical_selected=set(alive)
 complete=pd.read_csv(Q/'private_full_ancestry.tsv',sep='	',usecols=['research_id'],dtype=str); assert complete.research_id.is_unique
 keep=alive.isin(complete.research_id); removed=alive[~keep];alive=alive[keep];exclusions.append(set(removed)); assert set(alive)==set(d.index)
 flow.append({'stage':'original_analysis_after_PC_source_join','N_retained':len(alive),'N_excluded_at_stage':len(removed),'cases_retained':int(d.psoriasis_case.sum()),'controls_retained':int(d.psoriasis_case.eq(0).sum()),'definition':'Inner join to complete retained provider PC source; case >=2 distinct qualifying calendar dates, all remaining selected participants controls'})
 alive=d.index
 for label,bad,definition in [('provider_coverage_failure',~d.coverage_ok,'Remove provider coverage-failure list'),('noncanonical_or_unavailable_ploidy',~d.canonical,'Remove provider ploidy other than XX or XY, including unavailable'),('canonical_sex_discordance',d.canonical & ~d.sex_concordant,'Remove canonical XX/XY discordant with analysis sex covariate'),('provider_QC_flag',d.provider_flag,'Remove remaining provider flagged samples')]:
  reject=bad.reindex(alive).to_numpy(dtype=bool);removed=alive[reject];alive=alive[~reject];exclusions.append(set(removed)); t=tally(d.loc[alive]);flow.append({'stage':label,'N_retained':t['N'],'N_excluded_at_stage':len(removed),'cases_retained':t['cases'],'controls_retained':t['controls'],'definition':definition})
 assert set(alive)==set(strict.index)
 union=set()
 for s in exclusions:assert not union.intersection(s);union.update(s)
 assert len(union)+len(alive)==len(up) and not union.intersection(alive)
 write('full_sequential_flow',pd.DataFrame(flow));checks['upstream_partition_and_original_SQL_selection']=True
 # Phenotype edge-case lifecycle; counts refer to condition records and distinct dates separately.
 life=[]
 for stage,ids in [('eligible_before_clinical_selection',cl.index),('after_original_SQL_selection',sorted(clinical_selected)),('original_analysis',d.index),('strict_analysis',strict.index)]:
  z=cl.loc[ids]
  for label,mask in [('exactly_one_qualifying_record',z.qualifying_records.eq(1)),('exactly_one_qualifying_date',z.qualifying_dates.eq(1)),('no_qualifying_dates',z.qualifying_dates.eq(0)),('two_or_more_qualifying_dates',z.qualifying_dates.ge(2)),('any_undated_qualifying_record',z.undated_records.gt(0)),('qualifying_record_without_qualifying_date',z.qualifying_records.gt(0)&z.qualifying_dates.eq(0))]:
   row={'stage':stage,'criterion':label,'N':int(mask.sum()),'cases':pd.NA,'controls':pd.NA}
   if stage in ['original_analysis','strict_analysis']:
    p=d.loc[z.index[mask]].psoriasis_case; row.update(cases=int(p.sum()),controls=int(p.eq(0).sum()))
   life.append(row)
 write('phenotype_record_and_date_lifecycle',pd.DataFrame(life))
 assert cl.loc[d.index].qualifying_dates.ge(2).astype(int).equals(d.psoriasis_case)
 # Full private ancestry table + migration matrix. Missing legacy labels do not mean excluded from pooled scan.
 ars=[]
 for population,z,col in [('original_legacy',d,'legacy_ancestry'),('original_complete',d,'complete_provider_ancestry'),('strict_complete',strict,'complete_provider_ancestry')]:
  for group,part in z.groupby(col,sort=True):
   tt=tally(part); ars.append({'population':population,'group':group,**tt,'selected_for_three_ancestry_followups':group in ['EUR','AFR','AMR'],'passes_followup_sample_gate':tt['cases']>=200 and tt['controls']>=500,'gate':'at least 200 cases and 500 controls; variant gates remain separate'})
 write('ancestry_all_groups',pd.DataFrame(ars)); write('legacy_to_complete_migration',pd.crosstab(d.legacy_ancestry,d.complete_provider_ancestry).reset_index())
 gap=[]
 for label,mask in [('legacy_EUR_AFR_AMR',d.legacy_ancestry.isin(['EUR','AFR','AMR'])),('other_legacy_groups',~d.legacy_ancestry.isin(['EUR','AFR','AMR','UNASSIGNED'])),('absent_from_legacy_file',d.legacy_ancestry.eq('UNASSIGNED'))]:gap.append({'group':label,**tally(d.loc[mask])})
 assert sum(r['N'] for r in gap)==len(d) and gap[1]['N']+gap[2]['N']==174094
 write('historical_ancestry_gap',pd.DataFrame(gap))
 # Five subgroups: exact IDs, source hash, phenotype and covariate values.
 cfg=json.loads((N/'stratified_step1/strata.json').read_text());strata=[]
 for group,meta in cfg.items():
  mask=strict.SEX.eq(1 if group=='male' else 2) if group in ['female','male'] else strict.complete_provider_ancestry.eq(group)
  z=strict.loc[mask];p=N/'stratified_step1'/group; phs=read_tsv(p/'phenotype.tsv'); cvs=read_tsv(p/'covariates.tsv')
  assert set(z.index)==set(phs.index)==set(cvs.index)
  assert sha(p/'phenotype.tsv')==meta['phenotype_sha256'] and sha(p/'covariates.tsv')==meta['covariate_sha256']
  pd.testing.assert_series_equal(phs.psoriasis_case.sort_index(),z.psoriasis_case.sort_index(),check_names=False)
  pd.testing.assert_frame_equal(cvs[meta['covariates']].sort_index(),z[meta['covariates']].sort_index(),check_names=False,check_dtype=False)
  assert tally(z)=={'N':meta['n'],'cases':meta['cases'],'controls':meta['controls']}
  strata.append({'group':group,**tally(z),'covariates':', '.join(meta['covariates']),'assignment':meta['group_assignment'],'denominator':'cohort before variant-specific QC/missingness; not variant-specific analyzed N'})
  for f in ['phenotype.tsv','covariates.tsv','VALIDATED.json']:sources[str(p/f)]={'sha256':sha(p/f),'bytes':(p/f).stat().st_size}
 write('five_followup_groups',pd.DataFrame(strata));checks['all_five_strata_exact_IDs_values_hashes']=True
 # Reproduce age and EHR summaries; zero is retained, not treated as missing.
 for population,z in [('original',d),('strict',strict)]:
  for label,col in [('age','age'),('EHR','ehr_span_years')]:
   s=summ(z,col); old=pd.read_csv(N/'cohort_audit'/(population+'_'+label+'_summary.private.tsv'),sep='	');pd.testing.assert_frame_equal(s,old,check_dtype=False);write(population+'_'+label+'_summary',s)
 checks['age_and_EHR_summaries_reproduced']=True
 obs=[]
 for case,z in strict.groupby('psoriasis_case'):
  for label,mask in [('has_EHR_flag',z.has_ehr_data.eq(1)),('EHR_span_available',z.ehr_span_years.notna()),('EHR_span_zero',z.ehr_span_years.eq(0)),('EHR_flag_and_two_visit_dates_and_365_day_span',z.has_ehr_data.eq(1)&z.visit_dates.ge(2)&z.visit_span_days.ge(365))]:
   obs.append({'phenotype':'case' if case else 'control','criterion':label,'N_meets':int(mask.sum()),'N_does_not_meet':len(z)-int(mask.sum()),'N_total':len(z)})
 write('strict_observation_criteria',pd.DataFrame(obs))
 for p in [O/'upstream_eligibility.private.pkl',O/'QUERY_RECEIPT.private.json',O/'upstream_eligibility_query.sql',N/'inputs/phenotype.tsv',N/'inputs/covariates.tsv',N/'stratified_step1/strata.json',N/'scripts/refit_stratified_step1.py']:
  sources[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
 save(O/'RECONCILIATION_INPUTS.private.json',sources)
 save(O/'RECONCILED.private.json',{'status':'validated_private_reconciliation','created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'original':tally(d),'strict':tally(strict),'script_sha256':sha(__file__),'versions':{p:importlib.metadata.version(p) for p in ['numpy','pandas']},'release':'HELD: exact small counts; coordinated disclosure review and approved export required','limits':['No new chart validation','EHR-observed control sensitivity is not silently substituted for the frozen primary','Full S1/S2/manuscript integration not yet done','Analysis subgroup sample gates do not demonstrate statistical power']})
 print('COHORT_TABLE_RECONCILIATION_VALIDATED',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'RECONCILIATION_FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise
