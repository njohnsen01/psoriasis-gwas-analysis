from pathlib import Path
import json, hashlib, datetime, traceback, importlib.metadata
import numpy as np, pandas as pd
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
A=N/'submission_assets_20260911'; I=N/'reconciliation_20260912'; O=I/'estimate_audit'
R=Path('/home/jupyter/psoriasis_reproducibility_20260910')
O.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def tsv(p):return pd.read_csv(p,sep=chr(9),dtype={'FID':str,'IID':str,'ID':str})
def main():
 assert not (O/'COMPLETE.private.json').exists()
 sources={}; stages=[]; receipts=[]; qcs=[]; rows=[]; qa=[]
 def pin(p):
  p=Path(p);h=sha(p);sources[str(p)]={'sha256':h,'bytes':p.stat().st_size};return h
 ph=tsv(N/'inputs/phenotype.tsv')
 totals={'pooled':{'n':len(ph),'cases':int(ph.psoriasis_case.sum()),'controls':int(ph.psoriasis_case.eq(0).sum())}}
 totals.update(json.loads((N/'stratified_step1/strata.json').read_text()))
 for p in [N/'inputs/phenotype.tsv',N/'inputs/covariates.tsv',N/'stratified_step1/strata.json']:pin(p)
 mp=R/'reconstructed_inputs_v2/variant_map.tsv';M=tsv(mp);pin(mp);historic=set(M.variant_id)
 xp=A/'X_lead_estimates.private.tsv';X=tsv(xp);pin(xp)
 xids={f'chrX:{int(r.GRCh38_position)}:{r.REF}:{r.ALT_effect_allele}' for r in X.itertuples()}
 numeric=['GENPOS','N','N_CASES','N_CONTROLS','BETA','SE','LOG10P','OR','CI95_Wald_low','CI95_Wald_high']
 for group in ['pooled','female','male','EUR','AFR','AMR']:
  for domain in ['autosomes','X']:
   models=['male02'] if domain=='autosomes' or group=='female' else ['male02','male01']
   for model in models:
    rp=A/f'{group}_{domain}_{model}_source_receipts.private.json';pin(rp);rr=json.loads(rp.read_text())
    if domain=='autosomes':
     fp=A/f'{group}_autosomal_significant.private.tsv.gz';tab=tsv(fp);pin(fp)
     fp=A/f'{group}_autosomal_reporting_groups.private.tsv';pin(fp)
     try:report=tsv(fp)
     except pd.errors.EmptyDataError:report=pd.DataFrame(columns=tab.columns)
     want=set(tab.ID)|historic
    else:
     tab=X.loc[X.group.eq(group)&X.model.eq(model)].copy()
     tab['ID']=[f'chrX:{int(r.GRCh38_position)}:{r.REF}:{r.ALT_effect_allele}' for r in tab.itertuples()];want=xids
    got=[];total_read=0;finite=0
    for j,rec in enumerate(rr,1):
     p=Path(rec['path']);h=pin(p);assert h==rec['sha256'],str(p)
     d=pd.read_csv(p,sep=r'\s+',dtype={'ID':str,'ALLELE0':str,'ALLELE1':str})
     assert len(d)==rec['tested'] and d.ID.is_unique and d.TEST.eq('ADD').all()
     v=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce')
     ok=np.isfinite(v).all(axis=1)&v.SE.gt(0)&v.LOG10P.ge(0)
     assert int(ok.sum())==rec['finite'] and int((~ok).sum())==rec['unavailable']
     total_read+=len(d);finite+=int(ok.sum())
     z=d.loc[d.ID.isin(want)].copy()
     if len(z):
      z['source_file']=str(p);z['source_sha256']=h;z['group']=group;z['domain']=domain;z['model']=model;got.append(z)
     if model=='male02':
      qpath=p.parent/('variant_QC.private.tsv.gz' if domain=='autosomes' else 'final_variant_QC.private.tsv.gz')
      if qpath.exists():
       q=tsv(qpath);q=q.loc[q.ID.isin(want)].copy()
       if len(q):q['group']=group;q['domain']=domain;q['QC_source_file']=str(qpath);q['QC_source_sha256']=pin(qpath);qcs.append(q)
     receipts.append({'group':group,'domain':domain,'model':model,**rec})
     if j%200==0:print('AUDIT_READING',group,domain,model,j,flush=True)
    g=pd.concat(got,ignore_index=True) if got else pd.DataFrame(columns=['ID','LOG10P'])
    assert g.ID.is_unique
    if len(g):
     assert g.N.eq(g.N_CASES+g.N_CONTROLS).all() and g.N.le(totals[group]['n']).all()
     assert g.N_CASES.le(totals[group]['cases']).all() and g.N_CONTROLS.le(totals[group]['controls']).all()
     parts=g.ID.str.split(':',expand=True)
     assert parts[1].astype(int).eq(g.GENPOS).all() and parts[2].eq(g.ALLELE0).all() and parts[3].eq(g.ALLELE1).all()
     g['P_from_REGENIE_LOG10P']=np.power(10.,-g.LOG10P)
     g['OR']=np.exp(g.BETA);g['CI95_Wald_low']=np.exp(g.BETA-1.95996398454*g.SE);g['CI95_Wald_high']=np.exp(g.BETA+1.95996398454*g.SE)
     g['N_total_primary_cohort']=totals['pooled']['n'];g['N_total_stratum']=totals[group]['n'];g['cases_total_stratum']=totals[group]['cases'];g['controls_total_stratum']=totals[group]['controls']
     g['effect_allele']=g.ALLELE1
     g['dosage_units']='diploid ALT dosage 0/1/2' if domain=='autosomes' else ('female 0/1/2' if group=='female' else 'female 0/1/2; non-PAR male 0/'+('2' if model=='male02' else '1')+'; PAR diploid')
     g['N_definition']='REGENIE reported N; observed nonmissing-genotype counts retained separately in QC'
     g['P_definition']='exact REGENIE LOG10P; not reconstructed from beta/SE'
     g['CI_definition']='95% Wald from the same beta/SE row, not profile-likelihood intervals'
     rows.append(g)
    gi=g.set_index('ID')
    if domain=='autosomes':
     assert set(tab.ID)==set(g.loc[g.LOG10P.gt(-np.log10(5e-8))].ID)
     for label,t in [('significant',tab),('reporting_groups',report)]:
      if len(t):
       z=gi.loc[t.ID]
       assert z.ALLELE0.tolist()==t.ALLELE0.tolist() and z.ALLELE1.tolist()==t.ALLELE1.tolist()
       for col in numeric:assert np.allclose(z[col].to_numpy(float),t[col].to_numpy(float),rtol=2e-12,atol=1e-14,equal_nan=True),(group,label,col)
      qa.append({'group':group,'domain':domain,'model':model,'table':label,'rows':len(t),'source_values_match':True})
    else:
     t=tab.loc[tab.BETA.notna()];assert set(t.ID)==set(g.ID)
     if len(t):
      z=gi.loc[t.ID]
      for col in numeric[1:]:assert np.allclose(z[col].to_numpy(float),t[col].to_numpy(float),rtol=2e-12,atol=1e-14,equal_nan=True),(group,model,col)
     qa.append({'group':group,'domain':domain,'model':model,'table':'X_targets','rows':len(tab),'estimated':len(t),'source_values_match':True})
    stages.append({'group':group,'domain':domain,'model':model,'files':len(rr),'tested':total_read,'finite':finite,'target_rows':len(g)})
    save(O/'PROGRESS.private.json',stages);print('AUDIT_STAGE_VALIDATED',group,domain,model,flush=True)
 pd.concat(rows,ignore_index=True).to_csv(O/'exact_association_rows.private.tsv.gz',sep=chr(9),index=False)
 q=pd.concat(qcs,ignore_index=True);assert not q.duplicated(['group','domain','ID']).any();q.to_csv(O/'target_QC.private.tsv.gz',sep=chr(9),index=False)
 pd.DataFrame(qa).to_csv(O/'existing_table_values_QA.private.tsv',sep=chr(9),index=False)
 save(O/'source_receipts.private.json',receipts);save(O/'input_files.private.json',sources)
 save(O/'COMPLETE.private.json',{'status':'validated','created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':stages,'source_values_all_match':True,'script_sha256':sha(__file__),'versions':{p:importlib.metadata.version(p) for p in ['numpy','pandas']},'release':'HELD: private outputs; coordinated disclosure and approved export required','limits':['Follow-up model/input audit is separate','REGENIE N is distinct from observed genotype calls','Physical reporting groups do not establish independent loci']})
 print('EXACT_ESTIMATE_AUDIT_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise
