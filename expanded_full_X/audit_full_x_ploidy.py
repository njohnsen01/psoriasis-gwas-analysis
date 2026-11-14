"""Read-only whole-scan check of retained male genotype coding in original PGENs."""
import os
from pathlib import Path
import pandas as pd,numpy as np,json,sys,datetime,time
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'full_x_ploidy_audit';O.mkdir(exist_ok=True)
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');Q=R.parent/'chen_scan_benchmark'
sys.path.insert(0,f'{os.environ["PGENLIB_SITE_DIR"]}');import pgenlib
c=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});assert len(c)==357530
records=[]
for pv in sorted((B/'X_projected_source').glob('chunk_*/input.pvar')):
 out=O/f'{pv.parent.name}.private.json'
 if out.exists():records.append(json.loads(out.read_text()));continue
 passed=set((B/'results/X_source'/pv.parent.name/'passed_variants.txt').read_text().split());v=pd.read_csv(pv,sep='\t',skiprows=sum(1 for line in pv.open() if line.startswith('##')));idx=v.index[v.ID.isin(passed)&(v.POS>2781479)&(v.POS<155701383)].to_numpy()
 if not len(idx):continue
 ps=pd.read_csv(pv.with_suffix('.psam'),sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});ps['row']=np.arange(len(ps));d=ps[['FID','IID','SEX','row']].merge(c[['FID','IID','S']],on=['FID','IID'],validate='one_to_one');assert len(d)==357530 and np.array_equal(d.S,(d.SEX==1).astype(int));mi=d.loc[d.S==1,'row'].to_numpy();arr=np.empty(len(ps),np.int8);bad=[];missing=0;called=0;t=time.time()
 reader=pgenlib.PgenReader(str(pv.with_suffix('.pgen')).encode(),raw_sample_ct=len(ps))
 for i in idx:
  reader.read(int(i),arr);g=arr[mi];invalid=np.isfinite(g)&(g>=0)&~np.isin(g,[0,2]);n=int(invalid.sum())
  if n:bad.append({'ID':v.ID.iloc[i],'invalid_male_cells':n})
  missing+=int((g<0).sum());called+=int((g>=0).sum())
 reader.close();rec={'chunk':pv.parent.name,'selected_sites':len(idx),'male_participants':len(mi),'bad_sites':bad,'missing_male_cells':missing,'called_male_cells':called,'seconds':time.time()-t};out.write_text(json.dumps(rec,indent=2)+'\n');records.append(rec);print('PLOIDY_AUDIT_CHUNK',pv.parent.name,len(idx),'valid' if not bad else 'FAILED',flush=True)
assert sum(r['selected_sites'] for r in records)==334926
receipt={'status':'passed' if not any(r['bad_sites'] for r in records) else 'FAILED_INVALID_MALE_CALLS','scope':'All 334926 QC-passed non-PAR X variants; every retained male in the 357530-person analysis cohort; original PGEN files read without modification','variants':334926,'chunks':len(records),'invalid_nonmissing_male_encoding_present':any(bool(r['bad_sites']) for r in records),'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'GLOBAL_PLOIDY_AUDIT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True);assert receipt['status']=='passed'
