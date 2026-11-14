"""Prespecified Firth sensitivity for full-rank failed ML tests, protected outputs."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json,time,hashlib,sys,concurrent.futures,datetime
from scipy.special import expit
from scipy.linalg import cho_factor,cho_solve,eigvalsh
from scipy.linalg.blas import dsyrk
from scipy.stats import chi2
import rescue_x_stable as st
R,Q,O,B=st.R,st.Q,st.O,st.B

def evaluate(A,y,b,free,derivative=True):
 eta=A@b;p=expit(eta);w=p*(1-p);Aw=np.asfortranarray(A*np.sqrt(w)[:,None]);I=dsyrk(1.,Aw,trans=1,lower=0);cf=cho_factor(I,lower=False,check_finite=False);pll=float(np.sum(y*eta-np.logaddexp(0.,eta))+np.log(np.diag(cf[0])).sum())
 if not derivative:return pll
 invI=cho_solve(cf,np.eye(len(b)),check_finite=False);h=np.sum((A@invI)*A,axis=1)*w;U=A.T@(y-p+h*(.5-p));If=np.triu(I)+np.triu(I,1).T;step=np.zeros(len(b));step[free]=np.linalg.solve(If[np.ix_(free,free)],U[free]);return pll,U,step,If

def solve_firth(A,y,b=None,free=None,maxit=500):
 b=np.zeros(A.shape[1]) if b is None else b.copy();free=np.arange(len(b)) if free is None else np.array(free);last=0.
 for it in range(maxit):
  ll,U,step,I=evaluate(A,y,b,free);score=float(np.max(np.abs(U[free])));dec=float(U@step)
  if score<=1e-6 and np.max(np.abs(step))<=1e-7 and dec<=1e-10:
   ev=eigvalsh(I,check_finite=False);cond=float(ev[-1]/ev[0]);assert ev[0]>0 and cond<1e10,'firth_ill_conditioned';return b,dict(ll=ll,score=score,iterations=it,decrement=dec,condition=cond)
  assert np.all(np.isfinite(step)),'firth_nonfinite_step'
  alpha=min(1.,5./max(np.max(np.abs(step)),5.))
  for hs in range(50):
   bn=b+alpha*step;ln=evaluate(A,y,bn,free,False)
   if ln>=ll-1e-9:break
   alpha*=.5
  else:raise ValueError('firth_line_search_failed')
  b=bn
 raise ValueError('firth_not_converged')

def one(item):
 chunk,selected,stage=item;folder=B/chunk;outdir=O/stage;outdir.mkdir(exist_ok=True);out=outdir/(chunk+'.private.tsv')
 if out.exists():return pd.read_csv(out,sep='\t')
 pv=folder/'input.pvar';v=pd.read_csv(pv,sep='\t',skiprows=sum(l.startswith('##') for l in pv.open()));v['index']=np.arange(len(v));v=v[v.ID.isin(selected)];assert len(v)==len(selected)
 s=pd.read_csv(folder/'input.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});s['idx']=np.arange(len(s));a=st.d[['FID','IID','S']].merge(s[['FID','IID','SEX','idx']],on=['FID','IID'],sort=False,validate='one_to_one');assert len(a)==len(st.d) and np.array_equal(a.IID,st.d.IID) and np.array_equal(a.S,(a.SEX==1).astype(int));ix=a.idx.to_numpy();g=np.empty(len(s),np.int8);rows=[]
 with st.pgenlib.PgenReader(str(folder/'input.pgen').encode(),raw_sample_ct=len(s)) as rd:
  for row in v.itertuples():
   t0=time.time();rec=dict(ID=row.ID,status='pending')
   try:
    rd.read(row.index,g);z=g[ix];details=[];stats=[];N=0
    for sex in [0,1]:
     gs=z[st.sexinds[sex]];mask=gs>=0;gs=gs[mask];C=st.Cs[sex][mask];y=st.Ys[sex][mask];N+=len(y);levels=[0,1,2] if sex==0 else [0,2];ct=np.array([np.sum(gs==k) for k in levels]);assert np.all(ct>0),'rank_deficient_genotype';ref=levels[np.argmax(ct)];tested=[k for k in levels if k!=ref];G=np.column_stack([(gs==k).astype(float) for k in tested]);G/=G.std(axis=0);A=np.asfortranarray(np.column_stack([C,G]));init=np.r_[st.nulls[sex],np.zeros(len(tested))]
     b,full=solve_firth(A,y,init);start=b.copy();start[13:]=0;nb,null=solve_firth(A,y,start,free=np.arange(13));stat=2*(full['ll']-null['ll']);assert stat>=-1e-7 and np.isfinite(stat),'negative_PLRT';stats.append(max(0.,stat));details += [full,null]
    stat=sum(stats);rec.update(status='estimated',FIRTH_CHISQ_3DF=stat,P_FIRTH_3DF=float(chi2.sf(stat,3)),n_complete=N,max_score=max(x['score'] for x in details),max_condition=max(x['condition'] for x in details),max_decrement=max(x['decrement'] for x in details),iterations=max(x['iterations'] for x in details))
   except Exception as e:rec['status']=str(e)
   rec['seconds']=time.time()-t0;rows.append(rec)
 tmp=out.with_suffix('.tmp');pd.DataFrame(rows).to_csv(tmp,sep='\t',index=False);tmp.replace(out);return pd.DataFrame(rows)

def source_jobs(ids,stage):
 jobs=[];found=[]
 for pv in sorted(B.glob('chunk_*/input.pvar')):
  v=pd.read_csv(pv,sep='\t',skiprows=sum(l.startswith('##') for l in pv.open()));got=list(v[v.ID.isin(ids)].ID)
  if got:jobs.append((pv.parent.name,got,stage));found+=got
 assert len(found)==len(ids) and set(found)==set(ids);return jobs

def main():
 stage=sys.argv[1];assert stage in ['firth_benchmark','firth_recovery'];assert (O/'FIRTH_RESCUE_ADDENDUM.json').exists();st.setup()
 if stage=='firth_benchmark':
  lead=pd.read_csv(R/'full_x_followup/followup_named.csv');assert (lead.status=='estimated').all();diag=pd.read_csv(O/'diagnostics.private.tsv',sep='\t');extra=[]
  for lab,a in diag[diag.empty_outcome_cell & diag.full_genotype_rank].groupby('original_error'):
   ix=sorted(a.ID,key=lambda x:hashlib.sha256(('Firthrescuevalidation|'+x).encode()).hexdigest())[:4];extra+=ix
  ids=list(dict.fromkeys(list(lead.ID)+extra));(O/'firth_validation_selection.json').write_text(json.dumps({'ids':ids,'new_separated_controls':extra},indent=2)+'\n')
 else:
  assert json.loads((O/'FIRTH_BENCHMARK_VALIDATED.json').read_text())['status']=='passed';assert (O/'RECOVERY_COMPUTED.json').exists();a=pd.read_csv(O/'all_rescue_outcomes.private.tsv',sep='\t');audit=pd.read_csv(R/'full_x_estimability_audit/diagnostics.private.tsv',sep='\t');additional=audit[(audit.original_error=='.') & audit.empty_outcome_cell & audit.full_genotype_rank];ids=list(dict.fromkeys(list(a[a.status!='estimated'].ID)+list(additional.ID)));(O/'firth_recovery_scope.private.json').write_text(json.dumps({'ids':ids,'original_failed_subset':int((a.status!='estimated').sum()),'additional_native_valid_separated':len(additional)},indent=2)+'\n')
 jobs=source_jobs(ids,stage)
 with concurrent.futures.ThreadPoolExecutor(max_workers=8 if stage=='firth_benchmark' else 24) as ex:rows=list(ex.map(one,jobs))
 a=pd.concat(rows,ignore_index=True);assert len(a)==len(ids) and a.ID.is_unique;a.to_csv(O/(stage+'.private.tsv'),sep='\t',index=False)
 rec={'stage':stage,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'attempted':len(a),'status_counts':a.status.value_counts().to_dict(),'seconds_sum':float(a.seconds.sum()),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'status':'computed_pending_independent_review'}
 if stage=='firth_benchmark':
  x=a.merge(lead[['ID','firth_stat']],on='ID');rec['prior_R_controls']=len(x);rec['max_prior_R_stat_difference']=float(np.max(np.abs(x.FIRTH_CHISQ_3DF-x.firth_stat))) if 'FIRTH_CHISQ_3DF' in x else None
 (O/(stage.upper()+'_COMPUTED.json')).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
if __name__=='__main__':main()
