from pathlib import Path
import json,os,subprocess,time,traceback
import pandas as pd
N=Path(__file__).resolve().parents[1];O=N/'X_projected_source';done=set();env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4')
try:
 while True:
  if (O/'FAILED.json').exists():raise RuntimeError('X projection failed; inspect its FAILED.json')
  for folder in sorted(O.glob('chunk_*')):
   if folder.name in done or not (folder/'COMPLETE.json').exists():continue
   try:meta=json.loads((folder/'COMPLETE.json').read_text())
   except json.JSONDecodeError:continue
   assert meta['status']=='projected_and_count_QC_complete'
   if (folder/'HWE_COMPLETE.json').exists():done.add(folder.name);continue
   q=pd.read_csv(folder/'variant_QC.private.tsv.gz',sep='\t').set_index('ID')
   for an in ['EUR','AFR','AMR']:
    cmd=['plink2','--pfile',str(folder/'input'),'--split-par','hg38','--keep',str(N/'inputs'/(an+'.female_controls.keep.tsv')),'--hardy','--threads','4','--memory','4096','--out',str(folder/(an+'.female_controls'))]
    (folder/(an+'.hwe.command.json')).write_text(json.dumps(cmd,indent=2))
    with (folder/(an+'.hwe.console.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
    if r.returncode:raise RuntimeError(folder.name+' '+an+' HWE failed')
    tables=[pd.read_csv(p,sep=r'\s+').set_index('ID') for p in [folder/(an+'.female_controls.hardy'),folder/(an+'.female_controls.hardy.x')] if p.exists()]
    ht=pd.concat(tables);assert ht.index.is_unique and set(ht.index)==set(q.index);q[an+'_female_control_HWE_P']=ht.P.reindex(q.index)
   q['HWE_PASS']=q[[an+'_female_control_HWE_P' for an in ['EUR','AFR','AMR']]].ge(1e-15).all(axis=1);q['PASS_except_site_FILTER']=q.count_call_maf_PASS&q.HWE_PASS
   q.to_csv(folder/'QC_with_HWE.private.tsv.gz',sep='\t')
   (folder/'HWE_COMPLETE.json').write_text(json.dumps({'variants':len(q),'passed_except_site_FILTER':int(q.PASS_except_site_FILTER.sum()),'site_FILTER_pending':True},indent=2));done.add(folder.name);print('X_HWE_COMPLETE',folder.name,int(q.PASS_except_site_FILTER.sum()),flush=True)
  if (O/'COMPLETE.json').exists():
   n=json.loads((O/'COMPLETE.json').read_text())['chunks']
   if len(done)==n:break
  time.sleep(15)
 (O/'ALL_HWE_COMPLETE.json').write_text(json.dumps({'status':'complete','chunks':len(done),'site_FILTER_pending':True}))
except Exception as e:
 (N/'status/X_HWE_FAILED.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc()}));raise
