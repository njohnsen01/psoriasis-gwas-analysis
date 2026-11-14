"""Select prespecified post-recovery leads, independently validate, preserve private inputs."""
from pathlib import Path
import json,datetime,subprocess,concurrent.futures
import pandas as pd,numpy as np
import rescue_x_stable as st
import rescue_x_firth as fx
R,O,B,Q=st.R,st.O,st.B,st.Q

def spaced(df,pcol):
 selected=[]
 for t in df.sort_values(pcol).itertuples():
  if all(abs(t.POS-p)>1000000 for p in selected):selected.append(t.POS)
  if len(selected)==10:break
 return list(df[df.POS.isin(selected)].sort_values(pcol).drop_duplicates('POS').ID)

def main():
 N=O/'final_lead_checks';N.mkdir(exist_ok=True);base=pd.read_csv(R/'full_x/full_X_joint_Wald.tsv.gz',sep='\t');audit=pd.read_csv(R/'full_x_estimability_audit/diagnostics.private.tsv',sep='\t');rec=pd.read_csv(O/'all_rescue_outcomes.private.tsv',sep='\t');firth=pd.read_csv(O/'firth_recovery.private.tsv',sep='\t');a=base.merge(audit[['ID','full_genotype_rank','empty_outcome_cell']],on='ID',validate='one_to_one');orig=a[(a.ERRCODE=='.')&a.full_genotype_rank&~a.empty_outcome_cell][['ID','POS','WALD_CHISQ_3DF','P_CHISQ_3DF']].copy();orig['method']='original_PLINK_Wald'
 recovered=rec[rec.status=='estimated'][['ID','WALD_CHISQ_3DF','P_CHISQ_3DF']].merge(base[['ID','POS']],on='ID',validate='one_to_one');recovered['method']='equivalent_stable_ML_Wald';combined=pd.concat([orig,recovered],ignore_index=True);assert combined.ID.is_unique and np.isfinite(combined.P_CHISQ_3DF).all();combined.to_csv(O/'combined_eligible_Wald.private.tsv',sep='\t',index=False)
 ff=firth[firth.status=='estimated'].merge(base[['ID','POS']],on='ID',validate='one_to_one');prior=['chrX:5995921:G:A','chrX:5996494:GCT:G','chrX:5998724:A:G'];waldids=list(dict.fromkeys(spaced(combined,'P_CHISQ_3DF')+list(combined[combined.P_CHISQ_3DF<5e-8].ID)+prior));fids=list(dict.fromkeys(spaced(ff,'P_FIRTH_3DF')+list(ff[ff.P_FIRTH_3DF<5e-8].ID)));ids=list(dict.fromkeys(waldids+fids));assert len(ids)<=100,'Unexpectedmanythresholdcrossers:scientificcalibrationreviewbeforeextensiveselectedfits'
 sel=pd.DataFrame({'ID':ids,'ml_eligible':[v in set(combined.ID) for v in ids],'selected_Wald':[v in waldids for v in ids],'selected_Firth':[v in fids for v in ids]});sel.assign(ml_eligible=sel.ml_eligible.map({True:'TRUE',False:'FALSE'})).to_csv(N/'selection.private.csv',index=False)
 (N/'SELECTION.json').write_text(json.dumps({'selected':len(sel),'Wald_selected':len(waldids),'Firth_selected':len(fids),'rule':'top10>1Mbspacing+allP<5e-8foreachappropriateanalysis,plus3priorNLGN4X;notreplication'},indent=2)+'\n')
 cov=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});ph=pd.read_csv(Q/'phenotype.private.tsv',sep='\t',dtype={'FID':str,'IID':str});d=cov.merge(ph,on=['FID','IID'],validate='one_to_one');sources={}
 for pv in B.glob('chunk_*/input.pvar'):
  v=pd.read_csv(pv,sep='\t',skiprows=sum(l.startswith('##') for l in pv.open()));v['index']=np.arange(len(v))
  for t in v[v.ID.isin(sel.ID)].itertuples():sources[t.ID]=(pv.parent,t.index)
 assert len(sources)==len(sel);arrays={};old=pd.read_csv(R/'full_x_followup/followup_named.csv').set_index('ID');reused=0
 for j,vid in enumerate(sel.ID):
  folder,index=sources[vid];s=pd.read_csv(folder/'input.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});s['idx']=np.arange(len(s));z=d[['FID','IID','S']].merge(s[['FID','IID','SEX','idx']],on=['FID','IID'],sort=False,validate='one_to_one');assert len(z)==len(d) and np.array_equal(z.IID,d.IID) and np.array_equal(z.S,(z.SEX==1).astype(int));g=np.empty(len(s),np.int8)
  with st.pgenlib.PgenReader(str(folder/'input.pgen').encode(),raw_sample_ct=len(s)) as reader:reader.read(index,g)
  arr=g[z.idx.to_numpy()].astype(float);arr[arr<0]=np.nan;arrays['g'+str(j)]=arr
  if False and vid in old.index and old.loc[vid,'status']=='estimated' and not (N/('g'+str(j)+'.csv')).exists():
   t=old.loc[vid];pd.DataFrame([dict(key='g'+str(j),ID=vid,status='estimated',warnings='',wald=t.wald,firth_stat=t.firth_stat,firth_p=t.firth_p)]).to_csv(N/('g'+str(j)+'.csv'),index=False);reused+=1
 pd.concat([d[['CASE','S']+['Q'+str(i) for i in range(1,25)]],pd.DataFrame(arrays)],axis=1).to_csv(N/'selected.private.csv',index=False);lib=next(R.parent.glob('**/logistf/DESCRIPTION')).parent.parent
 with (N/'R.console.log').open('w') as log:p=subprocess.run(['/opt/conda/envs/jupyter/bin/Rscript',str(R/'validate_rescue_native_exact_information.R'),str(N),str(lib)],stdout=log,stderr=subprocess.STDOUT)
 assert p.returncode==0;native=pd.read_csv(N/'native_summary.csv');assert (native.status=='estimated').all(),native.status.value_counts().to_dict()
 st.setup()
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:ml=pd.concat(list(ex.map(st.one,fx.source_jobs(list(sel[sel.ml_eligible].ID),'final_ml_check'))),ignore_index=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:fp=pd.concat(list(ex.map(fx.one,fx.source_jobs(ids,'final_firth_check'))),ignore_index=True)
 assert (ml.status=='estimated').all() and (fp.status=='estimated').all();mx=ml.merge(native[['ID','wald']],on='ID');px=fp.merge(native[['ID','firth_stat']],on='ID');dm=float(np.max(np.abs(mx.WALD_CHISQ_3DF-mx.wald)));df=float(np.max(np.abs(px.FIRTH_CHISQ_3DF-px.firth_stat)));assert dm<1e-4 and df<1e-4,(dm,df)
 final=sel.merge(native,on='ID',validate='one_to_one');final.to_csv(N/'verified_leads.private.tsv',sep='\t',index=False)
 receipt={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'lead_checks_passed_pending_final_review','selected':len(sel),'reused_old_checks_this_run':reused,'ML_verified':len(mx),'Firth_verified':len(px),'max_ML_stat_difference':dm,'max_Firth_stat_difference':df,'ML_threshold_crossers_in_combined_scan':int((combined.P_CHISQ_3DF<5e-8).sum()),'Firth_sensitivity_threshold_crossers':int((ff.P_FIRTH_3DF<5e-8).sum()),'limitation':'Selectedverificationisnotindependentreplication;FirthsensitivityisnotprimaryWald.'};(O/'FINAL_LEAD_CHECKS_COMPUTED.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
