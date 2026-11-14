from pathlib import Path
import os,json,traceback,multiprocessing,concurrent.futures
from run_autosomes import N,process,step1,save,sha
def worker(t,pred):
 return process(t,Path(pred))
def main():
 tasks=json.loads((N/'autosomal_tasks.json').read_text());assert len(tasks)==1104
 keys=['chr'+t['chromosome']+'_'+t['stem'] for t in tasks];assert len(set(keys))==1104
 pred=step1();done=[];pending=[]
 for t,key in zip(tasks,keys):
  p=N/'results/autosomes'/key/'COMPLETE.json'
  if p.exists():done.append(process(t,pred))
  else:pending.append(t)
 save(N/'status/autosomal_parallel_resume.json',{'validated_existing_chunks':len(done),'pending_chunks':len(pending),'workers':6,'threads_per_worker':4,'analysis_code_sha256':sha(N/'scripts/run_autosomes.py'),'driver_sha256':sha(Path(__file__))})
 batches=[(pending[:2],2,'parallel_pilot'),(pending[2:],6,'running')]
 for batch,workers,state in batches:
  save(N/'status/autosomal_status.json',{'status':state,'completed':len(done),'expected':1104,'pending_batch':len(batch),'pid':os.getpid()})
  with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
   futures={pool.submit(worker,t,str(pred)):t for t in batch}
   for f in concurrent.futures.as_completed(futures):
    result=f.result();assert result['status']=='validated'
    if state=='parallel_pilot':assert result.get('failed_or_unavailable',0)==0
    done.append(result)
    save(N/'status/autosomal_status.json',{'status':state,'completed':len(done),'expected':1104,'last_completed':result['task'],'pid':os.getpid()})
    print('CHUNK_COMPLETE',result['task'],result['passed'],flush=True)
  if state=='parallel_pilot':save(N/'status/autosomal_parallel_pilot_validated.json',{'status':'validated','chunks':len(batch),'unchanged_analysis_function':True,'failed_tests':0})
 assert len(done)==1104 and len({x['task'] for x in done})==1104
 save(N/'status/autosomal_status.json',{'status':'complete','completed':1104,'expected':1104,'source_variants':sum(r['source_variants'] for r in done),'passed':sum(r['passed'] for r in done),'failed_or_unavailable':sum(r.get('failed_or_unavailable',0) for r in done),'scope':'historical fixed autosomal panel re-QC in strict cohort; not maximal source-wide autosomal discovery'})
 print('AUTOSOMES_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:
  save(N/'status/autosomal_status.json',{'status':'failed','error':str(e),'traceback':traceback.format_exc()});raise
