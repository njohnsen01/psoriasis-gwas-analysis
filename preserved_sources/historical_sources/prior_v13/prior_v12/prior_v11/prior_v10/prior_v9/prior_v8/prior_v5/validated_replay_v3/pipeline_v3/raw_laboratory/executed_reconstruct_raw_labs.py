from pathlib import Path
import ast,contextlib,hashlib,importlib.metadata as im,json,os,shutil,subprocess,sys,traceback
import numpy as np,pandas as pd
from google.cloud import bigquery,bigquery_storage
R=Path('/home/jupyter/psoriasis_reproducibility_20260910');V=R/'pipeline_v2';O=R/'raw_lab_reconstruction_v3';P=R/'fresh_lab_processing_v3'
def save(name,value):
 p=O/name;t=p.with_suffix(p.suffix+'.partial');t.write_text(json.dumps(value,indent=2));t.replace(p)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as handle:
  for b in iter(lambda:handle.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def row_fingerprints(path):
 arrays=[];counts={};cols=['person_id','measurement_concept_id','lab_name','lab_code','measurement_datetime','value_as_number','unit_name','unit_source_value','biomarker']
 for frame in pd.read_csv(path,sep=chr(9),dtype=str,chunksize=250000,usecols=cols):
  frame=frame[cols].copy()
  for c in cols:
   if c=='value_as_number':frame[c]=pd.to_numeric(frame[c],errors='raise').astype('float64')
   elif c=='measurement_datetime':frame[c]=pd.to_datetime(frame[c],utc=True,errors='raise',format='mixed').astype('int64')
   else:frame[c]=frame[c].fillna('').astype(str)
  for k,v in frame.biomarker.value_counts().items():counts[k]=counts.get(k,0)+int(v)
  a=np.empty(len(frame),dtype=[('h1','<u8'),('h2','<u8')])
  a['h1']=pd.util.hash_pandas_object(frame,index=False,hash_key='0123456789abcdef').to_numpy(np.uint64)
  a['h2']=pd.util.hash_pandas_object(frame,index=False,hash_key='fedcba9876543210').to_numpy(np.uint64)
  arrays.append(a)
 answer=np.concatenate(arrays);answer.sort(order=['h1','h2']);return answer,counts

def main():
 dry=json.loads((O/'dry_run.json').read_text());query=(O/'query.sql').read_text();assert sha(O/'query.sql')==dry['query_sha256']
 assert os.environ['WORKSPACE_CDR']==dry['source_CDR']=='wb-silky-artichoke-2408.C2024Q3R8'
 assert sha(V/'cells/cell_095.py')==dry['original_extraction_cell_sha256']
 ids=pd.read_csv(R/'replay_cohort_v2/psoriasis_pheno.tsv',sep=chr(9),dtype={'IID':str}).IID.drop_duplicates().tolist();assert len(ids)==392155
 assert shutil.disk_usage(O).free>128*1024**3
 save('versions.json',{'python':sys.version,'packages':{p:im.version(p) for p in ['numpy','pandas','google-cloud-bigquery','google-cloud-bigquery-storage','pyarrow']},'script_sha256':sha(Path(__file__))})
 client=bigquery.Client(project=os.environ['GOOGLE_PROJECT'])
 if (O/'query_job.private.json').exists():
  meta=json.loads((O/'query_job.private.json').read_text());job=client.get_job(meta['job_id'],location=meta['location'])
 else:
  cfg=bigquery.QueryJobConfig(use_query_cache=False,maximum_bytes_billed=45000000000,query_parameters=[bigquery.ArrayQueryParameter('cohort_ids','STRING',ids)])
  job=client.query(query,job_config=cfg)
  save('query_job.private.json',{'job_id':job.job_id,'location':job.location,'project':job.project,'maximum_bytes_billed':45000000000,'query_sha256':dry['query_sha256']})
 save('status.json',{'status':'executing_original_lab_query_with_single_cohort_binding'})
 raw=O/'psoriasis_biomarker_labs_raw.tsv'
 if not (O/'RAW_COMPLETE.json').exists():
  results=job.result();reader=bigquery_storage.BigQueryReadClient();n=0;part=0;tmp=O/'raw.partial.tsv'
  columns=['person_id','measurement_concept_id','lab_name','lab_code','lab_vocabulary','measurement_datetime','value_as_number','value_source_value','unit_concept_id','unit_name','unit_source_value','measurement_source_value','measurement_source_concept_id','source_lab_name','source_lab_code','source_vocabulary']
  buffer=[];buffer_n=0;idset=set(ids)
  def flush_frames(frames):
   nonlocal n,part
   frame=pd.concat(frames,ignore_index=True);assert frame.person_id.astype(str).isin(idset).all()
   frame.to_csv(tmp,sep=chr(9),index=False,mode='w' if part==0 else 'a',header=part==0)
   n+=len(frame);part+=1;save('status.json',{'status':'streaming_raw_labs','rows_written':n,'write_batches':part});print('RAW_LAB_STREAM',n,flush=True)
  for frame in results.to_dataframe_iterable(bqstorage_client=reader,max_queue_size=8):
   assert frame.columns.tolist()==columns;buffer.append(frame);buffer_n+=len(frame)
   if buffer_n>=100000:flush_frames(buffer);buffer=[];buffer_n=0
  if buffer:flush_frames(buffer)
  assert n==results.total_rows and n>0;tmp.replace(raw)
  save('RAW_COMPLETE.json',{'status':'raw_query_complete','rows':n,'bytes':raw.stat().st_size,'sha256':sha(raw),'total_bytes_processed':job.total_bytes_processed,'total_bytes_billed':job.total_bytes_billed,'cache_hit':job.cache_hit,'query_sha256':dry['query_sha256']})
 if not (O/'CLASSIFICATION_COMPLETE.json').exists():
  save('status.json',{'status':'classifying_with_original_cell_097'});os.chdir(O)
  with (O/'classification.private.log').open('w') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):exec(compile((V/'cells/cell_097.py').read_text(),str(V/'cells/cell_097.py'),'exec'),{'__name__':'__main__'})
  classified=O/'psoriasis_biomarker_labs_classified.tsv'
  save('CLASSIFICATION_COMPLETE.json',{'status':'complete','code_sha256':sha(V/'cells/cell_097.py'),'classified_sha256':sha(classified),'classified_bytes':classified.stat().st_size})
 classified=O/'psoriasis_biomarker_labs_classified.tsv'
 if not (P/'COMPLETE.json').exists():
  assert not P.exists() or not any(P.iterdir()),'Incomplete prior processing needs review'
  cmd=[sys.executable,str(V/'run.py'),'lab-processing','--cohort-dir',str(R/'replay_cohort_v2'),'--classified-labs',str(classified),'--out',str(P)]
  save('processing_command.json',cmd);save('status.json',{'status':'replaying_original_lab_processing'})
  with (O/'processing.private.log').open('w') as log:result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
  assert result.returncode==0,'Laboratory processing failed; inspect log'
 save('status.json',{'status':'validating_classified_measurements_and_person_level_endpoints'})
 oldroot=Path(f'{os.environ["PSO_WORKSPACE_MOUNT"].rstrip("/")}')
 before,bcounts=row_fingerprints(oldroot/'psoriasis_biomarker_labs_classified.tsv');after,acounts=row_fingerprints(classified)
 rowmatch=np.array_equal(before,after);del before,after
 x=pd.read_csv(oldroot/'psoriasis_biomarker_person_level.tsv',sep=chr(9),dtype={'person_id':str}).set_index('person_id').sort_index()
 y=pd.read_csv(P/'psoriasis_biomarker_person_level.tsv',sep=chr(9),dtype={'person_id':str}).set_index('person_id').sort_index()
 keys=x.index.is_unique and y.index.is_unique and x.index.equals(y.index) and x.columns.equals(y.columns)
 values=bool(keys and np.allclose(x.to_numpy(float),y.to_numpy(float),rtol=1e-12,atol=0,equal_nan=True))
 validation={'status':'validated' if rowmatch and keys and values and bcounts==acounts else 'requires_review','classified_row_multiset_128bit_fingerprints_match':rowmatch,'fingerprint_normalization':'all nine classified columns; numeric values as float64; datetimes normalized to UTC nanoseconds; duplicate rows retained; two independent 64-bit pandas hashes sorted lexicographically','biomarker_counts_match':bcounts==acounts,'biomarker_counts':acounts,'person_level_keys_match':bool(keys),'person_level_values_match':values,'numeric_tolerance':{'rtol':1e-12,'atol':0,'equal_nan':True},'rows_person_level':len(y),'columns_person_level':y.columns.tolist(),'new_person_level_sha256':sha(P/'psoriasis_biomarker_person_level.tsv'),'scope':'fresh CDR extraction plus exact original classification/processing; downstream model code was separately validated in v2; this run does not silently alter lab definition or units'}
 save('VALIDATION.json',validation);assert validation['status']=='validated','Fresh extraction differs; inspect secure validation and inputs'
 save('status.json',{'status':'complete','validation':'VALIDATION.json'});print('FRESH_RAW_LAB_RECONSTRUCTION_VALIDATED',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:
  save('FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});save('status.json',{'status':'failed','error':str(e)});raise
