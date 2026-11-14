from pathlib import Path
import ast,datetime,hashlib,json,os,sys,traceback
import pandas as pd
from google.cloud import bigquery
I=Path(__file__).resolve().parents[1];N=I.parent;B=N.parent;R=Path('/home/jupyter/psoriasis_reproducibility_20260910');P=R/'pipeline_v3';A=N/'cohort_audit'
O=I/'cohort';O.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+chr(10))
def main():
 assert os.environ['WORKSPACE_CDR']=='wb-silky-artichoke-2408.C2024Q3R8'
 manifest=json.loads((P/'manifest.json').read_text());sources={}
 for i in [3,4,6,11,12,13,14,15,16]:
  rel='cells/cell_'+str(i).zfill(3)+'.py';path=P/rel;assert sha(path)==manifest['files'][rel];sources[str(path)]=sha(path)
 tree=ast.parse((P/'cells/cell_003.py').read_text());node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='dataset_32891047_person_sql' for t in n.targets));ns={'os':os};exec(compile(ast.Module(body=[node],type_ignores=[]),'original_SQL_assignment','exec'),ns)
 original=ns['dataset_32891047_person_sql'];(O/'original_person_query.sql').write_text(original)
 cdr=os.environ['WORKSPACE_CDR']
 query='WITH original_selected AS ('+original+''')
SELECT s.person_id, s.has_whole_genome_variant, s.age_at_cdr,
 p.person_id IS NOT NULL AS has_person_record,
 p.sex_at_birth_concept_id,
 EXISTS(SELECT 1 FROM `'''+cdr+'''.death` de WHERE de.person_id=s.person_id) AS has_death_record,
 o.person_id IS NOT NULL AS selected_by_original_SQL
FROM `'''+cdr+'''.cb_search_person` s
LEFT JOIN `'''+cdr+'''.person` p USING(person_id)
LEFT JOIN original_selected o USING(person_id)
'''
 (O/'upstream_eligibility_query.sql').write_text(query)
 client=bigquery.Client(project=os.environ['GOOGLE_PROJECT'])
 dry=client.query(query,job_config=bigquery.QueryJobConfig(dry_run=True,use_query_cache=True));budget=100*1024**3
 save(O/'QUERY_ESTIMATE.private.json',{'estimated_bytes':dry.total_bytes_processed,'maximum_bytes_billed':budget,'query_sha256':sha(O/'upstream_eligibility_query.sql')})
 assert dry.total_bytes_processed<=budget,'Query exceeds documented ceiling; inspect before execution'
 print('UPSTREAM_QUERY_STARTED',flush=True)
 job=client.query(query,job_config=bigquery.QueryJobConfig(use_query_cache=True,maximum_bytes_billed=budget));base=job.to_dataframe();base.person_id=base.person_id.astype(str);assert base.person_id.is_unique
 base.to_pickle(O/'upstream_eligibility.private.pkl')
 save(O/'QUERY_RECEIPT.private.json',{'job_id':job.job_id,'location':job.location,'bytes_processed':job.total_bytes_processed,'bytes_billed':job.total_bytes_billed,'cache_hit':job.cache_hit,'CDR':cdr,'query_sha256':sha(O/'upstream_eligibility_query.sql')})
 path=B/'phenotype_validation_20260910/private_phenotype_audit.pkl';clinical=pd.read_pickle(path);clinical.person_id=clinical.person_id.astype(str);clinical=clinical.set_index('person_id');assert clinical.index.is_unique;sources[str(path)]=sha(path)
 eligible=base.has_person_record&base.has_whole_genome_variant.eq(1)&base.age_at_cdr.between(18,125)&~base.has_death_record&base.sex_at_birth_concept_id.isin([45880669,45878463])
 assert set(base.loc[eligible,'person_id'])==set(clinical.index),'Upstream eligibility differs from retained clinical audit'
 reconstructed=clinical.case_gate|(~clinical.any_excluded_condition&~clinical.any_excluded_drug)
 selected=set(base.loc[base.selected_by_original_SQL,'person_id']);assert set(clinical.index[reconstructed])==selected,'Original SQL and reconstructed phenotype gates disagree'
 ancestry_path=B/'phenotype_validation_20260910/private_full_ancestry.tsv';ancestry=pd.read_csv(ancestry_path,sep=chr(9),dtype={'research_id':str},usecols=['research_id','ancestry_pred']).set_index('research_id');assert ancestry.index.is_unique;sources[str(ancestry_path)]=sha(ancestry_path)
 retained=pd.read_csv(R/'replay_cohort_v2/psoriasis_pheno.tsv',sep=chr(9),dtype={'IID':str,'FID':str}).set_index('IID');assert retained.index.is_unique
 expected=selected&set(ancestry.index);assert expected==set(retained.index),'Original SQL plus PC-input membership differs from retained cohort'
 pd.testing.assert_series_equal(clinical.qualifying_dates.reindex(retained.index).ge(2).astype(int),retained.psoriasis_case.astype(int),check_names=False,check_dtype=False)
 assert len(retained)==392155 and int(retained.psoriasis_case.sum())==5805
 for path in [R/'replay_cohort_v2/psoriasis_pheno.tsv',R/'replay_cohort_v2/psoriasis_covariates.tsv',A/'cohort_clinical.private.pkl',A/'input_receipts.private.json',B/'phenotype_validation_20260910/phenotype_audit.sql']:sources[str(path)]=sha(path)
 save(O/'INPUTS.private.json',sources)
 save(O/'PREPARED.private.json',{'status':'validated','original_SQL_selection_matches_retained_clinical_rules':True,'clinical_eligibility_keys_match_new_source_query':True,'phenotype_and_PC_membership_match_retained_original':True,'source_notebook_sha256':manifest['notebook_sha256'],'script_sha256':sha(Path(__file__)),'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'release':'HELD'})
 print('UPSTREAM_COHORT_PREPARATION_VALIDATED',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise
