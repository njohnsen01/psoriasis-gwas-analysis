from pathlib import Path
import argparse,hashlib,json,os,sys,time,traceback
import numpy as np,pandas as pd
N=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(N.parent/'X_source_reanalysis_20260910/python_dependencies'));import pgenlib
from reconcile_X_site_filters import sha,save
LABELS=['EUR','AFR','AMR','female','male']
def initialize(folder):
 ps=pd.read_csv(folder/'input.psam',sep=chr(9),dtype={'IID':str})
 ph=pd.read_csv(N/'inputs/phenotype.tsv',sep=chr(9),dtype={'IID':str}).set_index('IID')
 assert len(ps)==len(ph)==386069 and ps.IID.is_unique and set(ps.IID)==set(ph.index)
 y=ph.psoriasis_case.reindex(ps.IID).to_numpy(int);male=ps.SEX.eq(1).to_numpy();assert np.isin(ps.SEX,[1,2]).all()
 an=np.full(len(ps),3,dtype=np.int8);groups=json.loads((N/'stratified_step1/strata.json').read_text())
 for i,label in enumerate(LABELS[:3]):
  ids=pd.read_csv(N/'stratified_step1'/label/'phenotype.tsv',sep=chr(9),dtype={'IID':str}).IID
  idx=ps.IID.isin(ids).to_numpy();assert np.all(an[idx]==3);an[idx]=i
 for label in LABELS:
  mask=(an==LABELS.index(label)) if label in LABELS[:3] else (male if label=='male' else ~male)
  ids=pd.read_csv(N/'stratified_step1'/label/'phenotype.tsv',sep=chr(9),dtype={'IID':str}).IID
  assert set(ps.loc[mask,'IID'])==set(ids) and int(mask.sum())==groups[label]['n'] and int(y[mask].sum())==groups[label]['cases']
 # Histogram cells represent ancestry x sex x case/control x genotype, including missing.
 base=((an.astype(np.int32)*2+male.astype(np.int32))*2+y)*4
 return ps,groups,base

def summarize(h,nonpar,ncase,ncontrol):
 vals={}
 for name,yi,total in [('case',1,ncase),('control',0,ncontrol)]:
  z=h[:,yi,:];fcalled=int(z[0,:3].sum());mcalled=int(z[1,:3].sum())
  fac=int(z[0,1]+2*z[0,2]);mac=int(z[1,2] if nonpar else z[1,1]+2*z[1,2]);ac=fac+mac;an=2*fcalled+(1 if nonpar else 2)*mcalled
  vals.update({name+'_AC':ac,name+'_AN':an,name+'_MAC':min(ac,an-ac),name+'_called':fcalled+mcalled,name+'_callrate':(fcalled+mcalled)/total})
 ac=vals['case_AC']+vals['control_AC'];an=vals['case_AN']+vals['control_AN'];mac=min(ac,an-ac);maf=mac/an if an else 0.
 vals.update({'total_MAC':mac,'MAF':maf,'count_call_maf_PASS':bool(ncase>=200 and ncontrol>=500 and vals['case_MAC']>=20 and vals['control_MAC']>=20 and mac>=100 and maf>=.01 and vals['case_callrate']>=.98 and vals['control_callrate']>=.98)})
 return vals

def process(folder,outroot):
 assert (folder/'COMPLETE.json').exists();out=outroot/folder.name;out.mkdir(parents=True,exist_ok=True)
 if (out/'COMPLETE.json').exists():
  rec=json.loads((out/'COMPLETE.json').read_text())
  for lab,digest in rec['output_sha256'].items():assert sha(out/(lab+'.count_QC.private.tsv.gz'))==digest
  return rec
 ps,groups,base=initialize(folder);pv=pd.read_csv(folder/'input.pvar',sep=chr(9));q=pd.read_csv(folder/'variant_QC.private.tsv.gz',sep=chr(9)).set_index('ID');assert pv.ID.tolist()==q.index.tolist()
 reader=pgenlib.PgenReader(str(folder/'input.pgen').encode());assert reader.get_raw_sample_ct()==len(ps) and reader.get_variant_ct()==len(pv)
 buf=np.empty(len(ps),np.int8);rows={lab:[] for lab in LABELS};pooled_match=0
 try:
  for j,v in enumerate(pv.itertuples(index=False)):
   reader.read(j,buf);assert np.isin(buf,[-9,0,1,2]).all()
   hist=np.bincount(base+np.where(buf<0,3,buf).astype(np.int32),minlength=64).reshape(4,2,2,4)
   nonpar=2781479<int(v.POS)<155701383
   if nonpar:assert int(hist[:,1,:,1].sum())==0,'Unexpected male non-PAR heterozygote'
   pooled=summarize(hist.sum(axis=0),nonpar,5705,380364)
   for key,value in pooled.items():
    if key=='count_call_maf_PASS':assert bool(q.loc[v.ID,key])==value
    else:assert np.isclose(float(q.loc[v.ID,key]),value,rtol=1e-12,atol=0,equal_nan=True),(v.ID,key)
   pooled_match+=1
   for label in LABELS:
    if label in LABELS[:3]:h=hist[LABELS.index(label)]
    else:
     h=hist.sum(axis=0);h[0 if label=='male' else 1,:,:]=0
    z=summarize(h,nonpar,groups[label]['cases'],groups[label]['controls'])
    rows[label].append({'ID':v.ID,'nonPAR':nonpar,'N_total_stratum':groups[label]['n'],'N_cases_total_stratum':groups[label]['cases'],**z})
 finally:reader.close()
 hashes={};passed={}
 for lab in LABELS:
  result=pd.DataFrame(rows[lab]);path=out/(lab+'.count_QC.private.tsv.gz');result.to_csv(path,sep=chr(9),index=False);hashes[lab]=sha(path);passed[lab]=int(result.count_call_maf_PASS.sum())
 rec={'status':'counts_validated','variants':len(pv),'independent_pooled_count_matches':pooled_match,'passed_subgroup_counts_only':passed,'output_sha256':hashes,'input_sha256':{ext:sha(folder/('input.'+ext)) for ext in ['pgen','pvar','psam']},'script_sha256':sha(Path(__file__)),'scope':'Subgroup count/call/MAF diagnostics only. Association additionally requires the pooled site FILTER and female-control ancestry HWE gates plus each subgroup Step 1 refit. Eligibility does not establish statistical power.'}
 save(out/'COMPLETE.json',rec);return rec

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');a=ap.parse_args()
 if a.pilot:
  r=process(N/'X_projection_pilot/chunk_0000',N/'X_stratified_QC_pilot');print('STRATIFIED_X_PILOT_COUNTS_VALIDATED',r['variants']);return
 root=N/'X_projected_source';out=N/'X_stratified_QC';out.mkdir(exist_ok=True);done={}
 while True:
  if (root/'FAILED.json').exists():raise RuntimeError('Projection failure requires review')
  for folder in sorted(root.glob('chunk_*')):
   if folder.name in done or not (folder/'COMPLETE.json').exists():continue
   save(out/'status.json',{'status':'running','completed_chunks':len(done),'current':folder.name,'pid':os.getpid()})
   done[folder.name]=process(folder,out);print('STRATIFIED_X_COUNTS_COMPLETE',folder.name,flush=True)
  if (root/'COMPLETE.json').exists():
   expected=json.loads((root/'COMPLETE.json').read_text())['chunks']
   if len(done)==expected:break
  time.sleep(30)
 save(out/'COMPLETE.json',{'status':'complete','chunks':len(done),'variants':sum(x['variants'] for x in done.values()),'groups':LABELS,'subgroup_QC_only':True})
if __name__=='__main__':
 try:main()
 except Exception as e:save(N/'status/X_stratified_QC_failed.json',{'error':str(e),'traceback':traceback.format_exc()});raise
