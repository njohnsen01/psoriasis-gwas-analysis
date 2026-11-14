from pathlib import Path
import os,json,hashlib,subprocess,concurrent.futures,datetime
import numpy as np,pandas as pd
import rescue_x_stable as st
import rescue_x_firth as fx
R,O,B,Q=st.R,st.O,st.B,st.Q

def main():
 N=O/'native_validation';N.mkdir(exist_ok=True);diag=pd.read_csv(O/'diagnostics.private.tsv',sep='\t');rows=[]
 for lab,a in diag[diag.empty_outcome_cell & diag.full_genotype_rank].groupby('original_error'):
  ids=sorted(a.ID,key=lambda x:hashlib.sha256(('Firthrescuevalidation|'+x).encode()).hexdigest())[:4];rows += [dict(ID=x,ml_eligible=False) for x in ids]
 a=diag[~diag.empty_outcome_cell & diag.full_genotype_rank];ids=sorted(a.ID,key=lambda x:hashlib.sha256(('MLrescuevalidation|'+x).encode()).hexdigest())[:8];rows += [dict(ID=x,ml_eligible=True) for x in ids];sel=pd.DataFrame(rows);sel.assign(ml_eligible=sel.ml_eligible.map({True:'TRUE',False:'FALSE'})).to_csv(N/'selection.private.csv',index=False)
 cov=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});ph=pd.read_csv(Q/'phenotype.private.tsv',sep='\t',dtype={'FID':str,'IID':str});d=cov.merge(ph,on=['FID','IID'],validate='one_to_one');sources={}
 for pv in B.glob('chunk_*/input.pvar'):
  v=pd.read_csv(pv,sep='\t',skiprows=sum(l.startswith('##') for l in pv.open()));v['index']=np.arange(len(v))
  for t in v[v.ID.isin(sel.ID)].itertuples():sources[t.ID]=(pv.parent,t.index)
 assert len(sources)==len(sel)
 for j,vid in enumerate(sel.ID):
  folder,index=sources[vid];s=pd.read_csv(folder/'input.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});s['idx']=np.arange(len(s));a=d[['FID','IID','S']].merge(s[['FID','IID','SEX','idx']],on=['FID','IID'],sort=False,validate='one_to_one');assert len(a)==len(d) and np.array_equal(a.IID,d.IID) and np.array_equal(a.S,(a.SEX==1).astype(int));g=np.empty(len(s),np.int8)
  with st.pgenlib.PgenReader(str(folder/'input.pgen').encode(),raw_sample_ct=len(s)) as reader:reader.read(index,g)
  z=g[a.idx.to_numpy()].astype(float);z[z<0]=np.nan;d['g'+str(j)]=z
 d[['CASE','S']+['Q'+str(i) for i in range(1,25)]+['g'+str(j) for j in range(len(sel))]].to_csv(N/'selected.private.csv',index=False)
 lib=next(R.parent.glob('**/logistf/DESCRIPTION')).parent.parent
 with (N/'R.console.log').open('w') as f:r=subprocess.run(['/opt/conda/envs/jupyter/bin/Rscript',str(R/'validate_rescue_native.R'),str(N),str(lib)],stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0
 st.setup();jobs=fx.source_jobs(list(sel[sel.ml_eligible].ID),'native_ml_control')
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:ml=pd.concat(list(ex.map(st.one,jobs)),ignore_index=True)
 ml.to_csv(N/'python_ml_controls.private.tsv',sep='\t',index=False);native=pd.read_csv(N/'native_summary.csv');x=ml.merge(native[['ID','wald','status']],on='ID',suffixes=('_python','_R'));assert (x.status_python=='estimated').all() and (x.status_R=='estimated').all();diff=float(np.max(np.abs(x.WALD_CHISQ_3DF-x.wald)));assert diff<1e-4
 rec={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'native_controls':len(native),'status_counts':native.status.value_counts().to_dict(),'ML_controls':len(x),'max_ML_Wald_difference':diff};(O/'NATIVE_VALIDATION_COMPUTED.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
if __name__=='__main__':main()
