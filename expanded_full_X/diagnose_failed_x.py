"""Protected-only diagnostic audit; no association-result rescue in this stage."""
from pathlib import Path
import os,json,datetime,sys,concurrent.futures
import pandas as pd,numpy as np
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');Q=R.parent/'chen_scan_benchmark';O=R/'full_x_rescue';B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/X_projected_source')
sys.path.insert(0,f'{os.environ["PGENLIB_SITE_DIR"]}');import pgenlib

def one(p):
 name=p.parent.name;out=O/(name+'.diagnostics.private.tsv')
 if out.exists():return pd.read_csv(out,sep='\t')
 v=pd.read_csv(p,sep='\t',skiprows=sum(l.startswith('##') for l in p.open()));v['index']=np.arange(len(v));v=v[v.ID.isin(failed.index)]
 if len(v)==0:return pd.DataFrame()
 s=pd.read_csv(p.parent/'input.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});s['idx']=np.arange(len(s))
 a=d.merge(s[['FID','IID','SEX','idx']],on=['FID','IID'],validate='one_to_one');assert len(a)==len(d) and np.array_equal(a.S,(a.SEX==1).astype(int));ix=a.idx.to_numpy();S=a.S.to_numpy(int);y=a.CASE.to_numpy(int)-1
 g=np.empty(len(s),np.int8);rows=[]
 with pgenlib.PgenReader(str(p.parent/'input.pgen').encode(),raw_sample_ct=len(s)) as rd:
  for t in v.itertuples():
   rd.read(t.index,g);z=g[ix];assert not np.any((S==1)&(z==1))
   cell=np.bincount((S[z>=0]*6+y[z>=0]*3+z[z>=0]),minlength=12).reshape(2,2,3)
   f=cell[0].sum(axis=0);m=cell[1].sum(axis=0);rank=bool(np.all(f>0) and m[0]>0 and m[2]>0)
   rows.append(dict(ID=t.ID,chunk=name,index=t.index,original_error=failed.loc[t.ID,'ERRCODE'],full_genotype_rank=rank,empty_outcome_cell=bool(np.any(cell[0]==0) or np.any(cell[1][:,[0,2]]==0)),female_hom_alt_zero=bool(f[2]==0),female_hom_ref_zero=bool(f[0]==0),min_genotype_cell=int(min(f.min(),m[[0,2]].min())),min_outcome_cell=int(min(cell[0].min(),cell[1][:,[0,2]].min())),**{f'cell_{i}_{j}_{k}':int(cell[i,j,k]) for i in range(2) for j in range(2) for k in range(3)}))
 outtmp=out.with_suffix('.tmp');pd.DataFrame(rows).to_csv(outtmp,sep='\t',index=False);outtmp.replace(out);return pd.DataFrame(rows)

if __name__=='__main__':
 O.mkdir(exist_ok=True)
 x=pd.read_csv(R/'full_x/full_X_joint_Wald.tsv.gz',sep='\t');failed=x[x.ERRCODE!='.'].set_index('ID');assert len(failed)==21730
 cov=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});ph=pd.read_csv(Q/'phenotype.private.tsv',sep='\t',dtype={'FID':str,'IID':str});d=cov[['FID','IID','S']].merge(ph,on=['FID','IID'],validate='one_to_one');assert len(d)==357530
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:rows=list(ex.map(one,sorted(B.glob('chunk_*/input.pvar'))))
 a=pd.concat(rows,ignore_index=True);assert len(a)==21730 and a.ID.is_unique;a.to_csv(O/'diagnostics.private.tsv',sep='\t',index=False)
 rec={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_tests':len(a),'privacy':'Counts below are numbers of tests, not participant counts. Exact genotype/outcome counts remain protected.','by_original_error':{}}
 for label,g in a.groupby('original_error'):
  rec['by_original_error'][label]={'tests':len(g),'full_genotype_rank':int(g.full_genotype_rank.sum()),'rank_deficient':int((~g.full_genotype_rank).sum()),'empty_outcome_cell':int(g.empty_outcome_cell.sum()),'female_hom_alt_zero':int(g.female_hom_alt_zero.sum()),'female_hom_ref_zero':int(g.female_hom_ref_zero.sum()),'minimum_genotype_cell_gt20':int((g.min_genotype_cell>20).sum()),'minimum_case_control_cell_gt20':int((g.min_outcome_cell>20).sum())}
 (O/'DIAGNOSTICS_COMPLETE.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
