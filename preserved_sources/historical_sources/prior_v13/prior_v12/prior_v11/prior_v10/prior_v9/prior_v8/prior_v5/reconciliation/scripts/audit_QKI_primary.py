from pathlib import Path
import json,hashlib,traceback
import pandas as pd
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';A=N/'submission_assets_20260911';O=I/'QKI_primary_audit';O.mkdir(exist_ok=True);vid='chr6:163803171:CTTATT:C'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 rows=[];qc=[];region=[];sources=[]
 for group in ['pooled','female','male','EUR','AFR','AMR']:
  f=A/(group+'_autosomes_male02_source_receipts.private.json');rr=json.loads(f.read_text());rr=[r for r in rr if '/chr6_' in r['path']];assert rr
  for rec in rr:
   p=Path(rec['path']);assert sha(p)==rec['sha256'];t=pd.read_csv(p,sep=r'\s+');z=t.loc[t.ID.eq(vid)].copy();z['group']=group;z['source']=str(p);z['source_sha256']=rec['sha256'];rows.append(z)
   if group=='pooled':
    z=t.loc[t.GENPOS.between(163303171,164303171)].copy();z['source']=str(p);z['source_sha256']=rec['sha256'];region.append(z)
   q=p.parent/'variant_QC.private.tsv.gz';assert q.exists();v=pd.read_csv(q,sep=chr(9));z=v.loc[v.ID.eq(vid)].copy();z['group']=group;z['source']=str(q);z['source_sha256']=sha(q);qc.append(z);sources.append({'path':str(q),'sha256':sha(q)})
 t=pd.concat(rows,ignore_index=True);q=pd.concat(qc,ignore_index=True);r=pd.concat(region,ignore_index=True);assert not t.group.duplicated().any() and not q.group.duplicated().any() and not r.ID.duplicated().any();t.to_csv(O/'QKI_exact_associations.private.tsv',sep=chr(9),index=False);q.to_csv(O/'QKI_QC.private.tsv',sep=chr(9),index=False);r.to_csv(O/'full_QKI_regional_associations.private.tsv.gz',sep=chr(9),index=False);(O/'COMPLETE.private.json').write_text(json.dumps({'status':'QKI exact source lookup completed','variant':vid,'association_groups':t.group.tolist(),'QC_groups':q.group.tolist(),'full_region_pooled_tests':len(r),'region_build':'GRCh38','region':[163303171,164303171],'LD':'pending','source_hashes':sources,'script_sha256':sha(Path(__file__)),'release':'HELD'},indent=2));print('QKI_AUDIT_COMPLETE',flush=True)
try:main()
except Exception as e:(O/'FAILED.private.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc()}));raise
