from pathlib import Path
import json,hashlib,sys,traceback,datetime
import numpy as np,pandas as pd
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');R=Path('/home/jupyter/psoriasis_reproducibility_20260910');I=N/'reconciliation_20260912';A=N/'submission_assets_20260911';O=I/'regional_X';O.mkdir(exist_ok=True);X=N.parent/'X_source_reanalysis_20260910';sys.path.insert(0,str(X/'python_dependencies'));import pgenlib
regions={'NLGN4X':(5498724,6498724,['chrX:5998724:A:G']),'PTCHD1-AS':(22417619,23417619,['chrX:22917619:C:T']),'Xq12':(67236589,68342269,['chrX:67736589:T:C','chrX:67842269:G:T'])}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def main():
 protocol={'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'regions_GRCh38':regions,'association':'Every finite QC-passing frozen pooled X test in the declared windows, both male coding models; no enriched extract','LD':'Pairwise complete Pearson genotype correlation squared in provider-pruned controls, separately EUR/AFR/AMR and female/male. Female diploid dosage r2 is unphased; male nonPAR r2 is haploid. No cross-sex pooled LD and no D-prime claim. Target ALT vs all other alleles; whole-genotype male heterozygotes masked.','LD_availability':'>=500 pair-complete controls, biological MAC>=20 for each member, each member call rate>=98%; unavailable LD remains visible. These are reporting thresholds, not demonstrated statistical power.','QC_failed_leads':'Reference positions and diagnostic LD tags only; never restored as primary associations','genotypes':'Retained target-ALT PGEN projections from unaffected provider source; whole-genotype male heterozygotes were masked before projection; projection source indices and allele indices retained','fine_mapping':'No credible-set or causal-variant inference','release':'HELD'};save('PROTOCOL.private.json',protocol)
 sources={};allqc=[];assocs=[]
 for model in ['male02','male01']:
  recs=json.loads((A/f'pooled_X_{model}_source_receipts.private.json').read_text())
  if isinstance(recs,dict):recs=recs.get('files',recs.get('sources',recs))
  assert isinstance(recs,list)
  for rec in recs:
   p=Path(rec['path']);assert sha(p)==rec['sha256'];sources[str(p)]={'sha256':rec['sha256'],'role':'frozen association'};d=pd.read_csv(p,sep=chr(92)+'s+');pos=d.GENPOS.astype(int)
   for region,(lo,hi,refs) in regions.items():
    z=d.loc[pos.between(lo,hi)].copy()
    if len(z):z['region']=region;z['model']=model;z['source_file']=str(p);z['source_sha256']=rec['sha256'];assocs.append(z)
   if model=='male02':
    qp=p.parent/'final_variant_QC.private.tsv.gz';q=pd.read_csv(qp,sep=chr(9));qp_sha=sha(qp);sources[str(qp)]={'sha256':qp_sha,'role':'final source QC'};qpos=q.ID.str.split(':').str[1].astype(int)
    for region,(lo,hi,refs) in regions.items():
     z=q.loc[qpos.between(lo,hi)].copy()
     if len(z):z['region']=region;z['source_file']=str(qp);z['source_sha256']=qp_sha;allqc.append(z)
 ass=pd.concat(assocs,ignore_index=True);qc=pd.concat(allqc,ignore_index=True);assert not qc.duplicated(['region','ID']).any();assert not ass.duplicated(['region','model','ID']).any()
 for model in ['male02','male01']:
  tested=set(ass.loc[ass.model.eq(model),'ID']);expected=set(qc.loc[qc.PASS,'ID']);assert tested==expected,(model,'full regional coverage')
 ass.to_csv(O/'full_regional_associations.private.tsv.gz',sep=chr(9),index=False);qc.to_csv(O/'all_regional_QC.private.tsv.gz',sep=chr(9),index=False)
 D=pd.read_pickle(R/'reconstructed_x_v2/X_inputs.private.pkl');D.index=D.index.astype(str);ps=pd.read_csv(N/'X_projected_source/chunk_0009/input.psam',sep=chr(9),dtype={'IID':str});sub=D.loc[D.unrelated&D.psoriasis_case.eq(0)&D.ancestry_pred.str.upper().isin(['EUR','AFR','AMR'])];ix=np.flatnonzero(ps.IID.isin(sub.index)).astype(np.uint32);d=sub.reindex(ps.IID.iloc[ix]);assert d.index.is_unique and d.SEX.eq(ps.SEX.iloc[ix].to_numpy()).all();male=d.SEX.eq(1).to_numpy();groups={(anc,sex):d.ancestry_pred.str.upper().eq(anc).to_numpy()&d.SEX.eq(1 if sex=='male' else 2).to_numpy() for anc in ['EUR','AFR','AMR'] for sex in ['female','male']}
 readers={};lookups={};buf=np.empty(len(ix),np.int8)
 def readg(r):
  chunk=Path(r.source_file).parent.name;folder=N/'X_projected_source'/chunk
  if chunk not in readers:
   pv=pd.read_csv(folder/'input.pvar',sep=chr(9));pp=pd.read_csv(folder/'input.psam',sep=chr(9),dtype={'IID':str});assert pp.IID.equals(ps.IID) and pp.SEX.equals(ps.SEX);assert pv.ID.is_unique
   lookups[chunk]=pd.Series(np.arange(len(pv)),index=pv.ID);readers[chunk]=pgenlib.PgenReader(str(folder/'input.pgen').encode(),sample_subset=ix)
   assert readers[chunk].get_variant_ct()==len(pv)
   for name in ['input.pgen','input.pvar','input.psam']:
    p=folder/name;sources[str(p)]={'sha256':sha(p),'role':'retained verified source projection for regional LD'}
  readers[chunk].read(int(lookups[chunk].loc[r.name]),buf);z=buf.astype(float);z[buf<0]=np.nan;assert np.isin(z[male&np.isfinite(z)],[0,2]).all();return z
 ld=[];refs_cache={}
 try:
  for region,(lo,hi,refs) in regions.items():
   z=qc.loc[qc.region.eq(region)].copy().set_index('ID');selected=z.loc[z.PASS|z.index.isin(refs)]
   for ref in refs:
    assert ref in z.index;refs_cache[ref]=readg(z.loc[ref]);rs={'chrX:5998724:A:G':'rs183018930','chrX:22917619:C:T':'rs5970952','chrX:67736589:T:C':'rs5964607','chrX:67842269:G:T':'rs5919457'}[ref];assert np.allclose(refs_cache[ref],d[rs],rtol=0,atol=0,equal_nan=True),'reference projection equality'
   for k,(vid,v) in enumerate(selected.iterrows()):
    g=readg(v)
    for ref in refs:
     base=refs_cache[ref]
     for (anc,sex),use in groups.items():
      gx=g[use];gy=base[use];ok=np.isfinite(gx)&np.isfinite(gy);n=int(ok.sum());ploidy=1 if sex=='male' else 2;mult=.5 if sex=='male' else 1.;acx=float(gx[ok].sum()*mult);acy=float(gy[ok].sum()*mult);an=ploidy*n;macx=min(acx,an-acx);macy=min(acy,an-acy);crx=float(np.isfinite(gx).mean());cry=float(np.isfinite(gy).mean());why=[]
      if n<500:why.append('pair_complete_n')
      if min(macx,macy)<20:why.append('MAC')
      if min(crx,cry)<.98:why.append('call_rate')
      r2=np.nan
      if not why:
       xx=gx[ok]-gx[ok].mean();yy=gy[ok]-gy[ok].mean();den=float(np.dot(xx,xx)*np.dot(yy,yy))
       if den<=0:why.append('zero_variance')
       else:r2=float(np.dot(xx,yy)**2/den);assert -.000001<=r2<=1.000001;r2=min(1.,max(0.,r2))
      ld.append({'region':region,'variant':vid,'reference':ref,'ancestry':anc,'sex':sex,'r2':r2,'status':'available' if not why else ','.join(why),'pair_complete_n':n,'variant_MAC':macx,'reference_MAC':macy,'variant_callrate':crx,'reference_callrate':cry,'variant_primary_PASS':bool(v.PASS),'reference_primary_PASS':bool(z.loc[ref,'PASS'])})
    if k%200==0:print('REGIONAL_LD',region,k,len(selected),flush=True)
 finally:
  for reader in readers.values():reader.close()
 ld=pd.DataFrame(ld);ld.to_csv(O/'ancestry_sex_control_LD.private.tsv.gz',sep=chr(9),index=False);save('input_files.private.json',sources)
 save('COMPLETE.private.json',{'status':'full_regional_association_and_LD_extraction_verified','association_records':len(ass),'QC_records':len(qc),'LD_records':len(ld),'windows':regions,'all_QC_passing_regional_tests_present':True,'source_coding':'female0/1/2 male0/2 genotypes; within-sex correlation invariant to male linear rescaling','plot_generation':'pending','release':'HELD','script_sha256':sha(Path(__file__))});print('REGIONAL_DATA_COMPLETE',flush=True)
try:main()
except Exception as e:save('FAILED.private.json',{'error':repr(e),'traceback':traceback.format_exc(),'release':'HELD'});raise
