from pathlib import Path
import json, hashlib, re, datetime, concurrent.futures, traceback
import numpy as np
import pandas as pd
N=Path(__file__).resolve().parents[1]
O=N/'submission_assets_20260911'
F=N/'freeze_review_20260912'
F.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+chr(10))
def check(v):
 p=Path(v['path']);assert sha(p)==v['sha256'],str(p)
 assert v['tested']==v['finite'] and v['unavailable']==0,str(p)
 log=p.with_name(p.name.replace('_psoriasis_case.regenie','.log'));s=log.read_text()
 assert 'End time:' in s,str(log)
 counts=re.findall(r'Number of failed tests\s*:\s*\((\d+)/(\d+)\)',s)
 assert counts and all(int(a)==0 for a,b in counts),str(log)
 ignored=re.findall(r'Number of ignored tests due to low MAC\s*:\s*(\d+)',s)
 assert ignored and all(int(a)==0 for a in ignored),str(log)
 warnings=[l for l in s.splitlines() if re.search(r'warning|error|fail|converg',l,re.I) and not re.search(r'Number of failed tests|--minMAC|fitting null',l,re.I)]
 return dict(v,bytes=p.stat().st_size,log=str(log),log_sha256=sha(log),warning_lines=warnings)
def main():
 sources={}
 for p in sorted(O.glob('*_source_receipts.private.json')):
  for v in json.loads(p.read_text()):
   if v['path'] in sources:assert sources[v['path']]==v
   sources[v['path']]=v
 print('SOURCE_FILES',len(sources),flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:results=list(ex.map(check,sources.values()))
 save(F/'source_and_convergence_audit.private.json',results)
 warnings=[v for v in results if v['warning_lines']]
 print('CHECKED_SOURCE_FILES',len(results),'LOGS_WITH_WARNINGS',len(warnings),flush=True)
 ranks=[]
 for p in sorted(O.glob('*_QQ_ranks.private.tsv.gz')):
  q=pd.read_csv(p,sep='	');r=q.full_output_rank.to_numpy();n=q.N_valid_tests.to_numpy();assert len(set(n))==1 and r[0]==1 and r[-1]==n[0] and (np.diff(r)>0).all()
  assert np.allclose(q.expected_log10P,-np.log10((r-.5)/n),rtol=0,atol=1e-12)
  assert (np.diff(q.observed_log10P)<=0).all() and np.isfinite(q.select_dtypes('number')).all().all()
  assert (q.null_95_low<=q.null_95_high).all()
  ranks.append({'path':str(p),'sha256':sha(p),'full_tests':int(n[0]),'display_points':len(q),'unit_weight_full_rank_expected_values_verified':True})
 save(F/'QQ_rank_audit.private.json',ranks)
 manifest=json.loads((O/'MANIFEST.json').read_text())
 for v in manifest['files']:
  p=O/v['name'];assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],str(p)
 files=[]
 for root in [N/'scripts',N/'diagnostics',N/'stratified_step1']:
  for p in sorted(root.rglob('*')):
   if p.is_file() and '__pycache__' not in p.parts:files.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
 freeze={'version':'corrected_autosomes_X_20260912_review_candidate','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'source integrity and log audit passed; author and release decisions pending','source_files':len(results),'qq_rank_files':len(ranks),'logs_with_warnings':len(warnings),'asset_manifest_sha256':sha(O/'MANIFEST.json'),'source_audit_sha256':sha(F/'source_and_convergence_audit.private.json'),'code_and_diagnostic_files':files,'scope':'autosomal historical panel re-QC, source-wide X and five overlapping diagnostic strata; no corrected Y analysis','release':'HELD; no export approval; no participant-level release','outstanding':'Review warning lines if any, complete cohort/input inventory, Y scope, result-table reconciliation, release/disclosure and author review'}
 save(F/'FREEZE_CANDIDATE.private.json',freeze)
 save(F/'COMPLETE.json',{'status':'validated','source_files':len(results),'rank_files':len(ranks),'logs_with_warnings':len(warnings),'freeze_manifest_sha256':sha(F/'FREEZE_CANDIDATE.private.json'),'release':'HELD'})
 print('FREEZE_AUDIT_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(F/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
