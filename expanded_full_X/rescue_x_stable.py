"""Same-model 3df Wald recovery via independent sex-specific categorical bases.
All individual records and cell counts remain in protected Verily.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,datetime,hashlib,time,concurrent.futures
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.linalg import qr,solve_triangular,cho_factor,cho_solve,eigvalsh
from scipy.linalg.blas import dsyrk
from scipy.stats import chi2
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');Q=R.parent/'chen_scan_benchmark';O=R/'full_x_rescue';B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/X_projected_source')
sys.path.insert(0,f'{os.environ["PGENLIB_SITE_DIR"]}');import pgenlib

def solve_ml(A,y,b=None,maxit=80):
 b=np.zeros(A.shape[1]) if b is None else b.copy()
 for it in range(maxit):
  eta=A@b;p=expit(eta);w=p*(1-p);score=A.T@(y-p);Xw=np.asfortranarray(A*np.sqrt(w)[:,None]);I=dsyrk(1.,Xw,trans=1,lower=0);cf=cho_factor(I,lower=False,check_finite=False);step=cho_solve(cf,score,check_finite=False)
  dec=float(score@step)
  if np.max(np.abs(score))<=1e-6 and dec<=1e-10 and np.max(np.abs(step))<=1e-7:
   fullI=np.triu(I)+np.triu(I,1).T;ev=eigvalsh(fullI,check_finite=False);cond=float(ev[-1]/ev[0]);assert ev[0]>0 and cond<1e10,'ill_conditioned_fisher'
   return b,cho_solve(cf,np.eye(len(b))),{'iterations':it,'score':float(np.max(np.abs(score))),'decrement':dec,'condition':cond,'ll':float(np.sum(y*eta-np.logaddexp(0.,eta)))}
  assert np.all(np.isfinite(step)) and np.max(np.abs(b))<100,'divergent_ML'
  ll=float(np.sum(y*eta-np.logaddexp(0.,eta)));a=1.
  for hs in range(40):
   bn=b+a*step;en=A@bn;ln=float(np.sum(y*en-np.logaddexp(0.,en)))
   if ln>=ll-1e-9:break
   a*=.5
  else:raise ValueError('line_search_failed')
  b=bn
 raise ValueError('not_converged')

def setup():
 global d,sexinds,Cs,Ys,nulls,baseline,alignment_cache
 cov=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});ph=pd.read_csv(Q/'phenotype.private.tsv',sep='\t',dtype={'FID':str,'IID':str});d=cov.merge(ph,on=['FID','IID'],validate='one_to_one');assert len(d)==357530
 baseline=pd.read_csv(R/'full_x/full_X_joint_Wald.tsv.gz',sep='\t').set_index('ID');sexinds={};Cs={};Ys={};nulls={};checks={};alignment_cache={}
 A=np.column_stack([np.ones(len(d)),d[['S']+[f'Q{i}' for i in range(1,25)]].to_numpy()])
 for sex in [0,1]:
  ix=np.where(d.S.to_numpy()==sex)[0];sexinds[sex]=ix;Z=A[ix];_,rr,piv=qr(Z,mode='economic',pivoting=True,check_finite=False);diag=np.abs(np.diag(rr));rank=int(np.sum(diag>diag[0]*1e-10));assert rank==13,(sex,rank)
  c=Z[:,piv[:13]];qq,tt=qr(c,mode='economic',check_finite=False);C=qq*np.sqrt(len(ix));res=float(np.max(np.abs(Z-C@np.linalg.lstsq(C,Z,rcond=None)[0])));assert res<1e-9,res
  Cs[sex]=np.asfortranarray(C);Ys[sex]=d.CASE.to_numpy(float)[ix]-1;nulls[sex]=solve_ml(Cs[sex],Ys[sex])[0];checks[str(sex)]={'nuisance_rank':rank,'span_max_error':res}
 (O/'NUISANCE_EQUIVALENCE.json').write_text(json.dumps(checks,indent=2)+'\n')

def one(item):
 chunk,selected,stage=item;folder=B/chunk;outdir=O/stage;outdir.mkdir(exist_ok=True);out=outdir/(chunk+'.private.tsv')
 if out.exists():return pd.read_csv(out,sep='\t')
 v=pd.read_csv(folder/'input.pvar',sep='\t',skiprows=sum(l.startswith('##') for l in (folder/'input.pvar').open()));v['index']=np.arange(len(v));v=v[v.ID.isin(selected)]
 assert len(v)==len(selected)
 s=pd.read_csv(folder/'input.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});s['idx']=np.arange(len(s));aligned=d[['FID','IID','S']].merge(s[['FID','IID','SEX','idx']],on=['FID','IID'],validate='one_to_one',sort=False);assert len(aligned)==len(d) and np.array_equal(aligned.IID,d.IID) and np.array_equal(aligned.S,(aligned.SEX==1).astype(int));ix=aligned.idx.to_numpy();g=np.empty(len(s),np.int8);rows=[]
 with pgenlib.PgenReader(str(folder/'input.pgen').encode(),raw_sample_ct=len(s)) as reader:
  for row in v.itertuples():
   t0=time.time();rec={'ID':row.ID,'original_error':baseline.loc[row.ID,'ERRCODE'],'status':'pending'}
   try:
    reader.read(row.index,g);z=g[ix];fits=[];stats=[];N=0
    for sex in [0,1]:
     gs=z[sexinds[sex]];mask=gs>=0;gs=gs[mask];y=Ys[sex][mask];C=Cs[sex][mask];N+=len(y);levels=[0,1,2] if sex==0 else [0,2];ct=np.array([[np.sum((gs==k)&(y==j)) for j in [0,1]] for k in levels]);assert np.all(ct.sum(axis=1)>0),'rank_deficient_genotype';assert np.all(ct>0),'separated_genotype_cell'
     ref=levels[np.argmax(ct.sum(axis=1))];tested=[k for k in levels if k!=ref];G=np.column_stack([(gs==k).astype(float) for k in tested]);G/=G.std(axis=0);A=np.asfortranarray(np.column_stack([C,G]));init=np.r_[nulls[sex],np.zeros(len(tested))];b,V,details=solve_ml(A,y,init);i=np.arange(13,len(b));stat=float(b[i]@np.linalg.solve(V[np.ix_(i,i)],b[i]));assert stat>=0 and np.isfinite(stat);stats.append(stat);fits.append(details)
    stat=sum(stats);rec.update(status='estimated',WALD_CHISQ_3DF=stat,P_CHISQ_3DF=float(chi2.sf(stat,3)),n_complete=N,max_score=max(t['score'] for t in fits),max_condition=max(t['condition'] for t in fits),max_decrement=max(t['decrement'] for t in fits),iterations=max(t['iterations'] for t in fits),female_wald=stats[0],male_wald=stats[1])
   except Exception as e:rec.update(status=str(e))
   rec['seconds']=time.time()-t0;rows.append(rec)
 tmp=out.with_suffix('.tmp');pd.DataFrame(rows).to_csv(tmp,sep='\t',index=False);tmp.replace(out);return pd.DataFrame(rows)

def main():
 O.mkdir(exist_ok=True);stage=sys.argv[1];assert stage in ['benchmark','recovery'];assert (O/'RESCUE_PROTOCOL.json').exists();setup()
 if stage=='benchmark':
  ids=[]
  for p in sorted((R/'full_x').glob('chunk_*/joint.tsv')):
   a=pd.read_csv(p,sep='\t');a=a[a.ERRCODE=='.'];a=a.assign(h=a.ID.map(lambda x:hashlib.sha256(('Xstablevalidation|'+x).encode()).hexdigest())).sort_values('h').head(1)
   ids += [(p.parent.name,list(a.ID),stage)]
 else:
  assert json.loads((O/'BENCHMARK_VALIDATED.json').read_text())['status']=='passed'
  diag=pd.read_csv(O/'diagnostics.private.tsv',sep='\t');assert len(diag)==21730;eligible=diag[diag.full_genotype_rank & ~diag.empty_outcome_cell];ids=[(chunk,list(g.ID),stage) for chunk,g in eligible.groupby('chunk')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:rr=list(ex.map(one,ids))
 a=pd.concat(rr,ignore_index=True);a.to_csv(O/(stage+'.private.tsv'),sep='\t',index=False);valid=a.status=='estimated'
 rec={'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':stage,'attempted':len(a),'status_counts':a.status.value_counts().to_dict(),'seconds_sum':float(a.seconds.sum()),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 if stage=='benchmark':
  assert valid.all(),a.status.value_counts().to_dict();err=np.abs(a.WALD_CHISQ_3DF.to_numpy()-baseline.loc[a.ID].WALD_CHISQ_3DF.to_numpy());rec['max_wald_difference']=float(err.max());assert err.max()<.005;rec['status']='passed';(O/'BENCHMARK_VALIDATED.json').write_text(json.dumps(rec,indent=2)+'\n')
 else:
  remaining=diag[~diag.ID.isin(a.ID)][['ID','original_error']].copy();remaining['status']=np.where(diag.set_index('ID').loc[remaining.ID,'full_genotype_rank'].to_numpy(),'separated_genotype_cell','rank_deficient_genotype');a=pd.concat([a,remaining],ignore_index=True);assert len(a)==21730 and a.ID.is_unique;a.to_csv(O/'all_rescue_outcomes.private.tsv',sep='\t',index=False);rec['all_status_counts']=a.status.value_counts().to_dict();rec['status']='computed_pending_independent_review';(O/'RECOVERY_COMPUTED.json').write_text(json.dumps(rec,indent=2)+'\n')
 print(json.dumps(rec),flush=True)
if __name__=='__main__':main()
