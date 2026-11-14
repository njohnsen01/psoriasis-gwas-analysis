from pathlib import Path
import ast,datetime,hashlib,json,os
from google.cloud import bigquery
P=Path('/home/jupyter/psoriasis_reproducibility_20260910/pipeline_v3');D=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912/supplement_staging_v5/methods_sources')
ids=set()
for n in ast.walk(ast.parse((P/'cells/cell_097.py').read_text())):
 if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id.endswith('_IDS'):ids.update(ast.literal_eval(n.value))
assert os.environ['WORKSPACE_CDR']=='wb-silky-artichoke-2408.C2024Q3R8'
q='SELECT concept_id,concept_name,concept_code,vocabulary_id,domain_id,standard_concept,invalid_reason FROM '+chr(96)+os.environ['WORKSPACE_CDR']+'.concept'+chr(96)+' WHERE concept_id IN ('+','.join(sorted(ids,key=int))+') ORDER BY concept_id'
(D/'lab_vocabulary_query.sql').write_text(q+chr(10));client=bigquery.Client(project=os.environ['GOOGLE_PROJECT']);dry=client.query(q,job_config=bigquery.QueryJobConfig(dry_run=True,use_query_cache=True));assert dry.total_bytes_processed<=1024**3
job=client.query(q,job_config=bigquery.QueryJobConfig(use_query_cache=True,maximum_bytes_billed=1024**3));df=job.result().to_dataframe();df.to_csv(D/'lab_vocabulary_current.private.tsv',sep=chr(9),index=False)
r={'queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'CDR':os.environ['WORKSPACE_CDR'],'query_sha256':hashlib.sha256(q.encode()).hexdigest(),'job_id':job.job_id,'estimated_bytes':dry.total_bytes_processed,'bytes_processed':job.total_bytes_processed,'bytes_billed':job.total_bytes_billed,'cache_hit':job.cache_hit,'requested_concepts':sorted(ids,key=int),'returned_concepts':list(map(str,df.concept_id)),'unreturned_concepts':sorted(ids-set(map(str,df.concept_id))),'release':'HELD'}
(D/'lab_vocabulary_query_receipt.private.json').write_text(json.dumps(r,indent=2)+chr(10));print('VOCABULARY_QUERY_VERIFIED',len(df),'returned; unmatched:',r['unreturned_concepts'])
