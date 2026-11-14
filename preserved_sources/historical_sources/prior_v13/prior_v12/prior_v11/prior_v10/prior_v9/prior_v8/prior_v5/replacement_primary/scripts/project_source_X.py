from pathlib import Path
import argparse,hashlib,json,os,shutil,sys,time,traceback
import numpy as np,pandas as pd
N=Path(__file__).resolve().parents[1];X=N.parent/'X_source_reanalysis_20260910';sys.path.insert(0,str(X/'python_dependencies'));import pgenlib
p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');a=p.parse_args()
O=N/('X_projection_pilot' if a.pilot else 'X_projected_source');O.mkdir(exist_ok=True)
c=pd.read_csv(N/'source_X_candidates/count_candidates.tsv.gz',sep='\t')
if a.pilot:c=pd.concat([c.iloc[:40],c.iloc[np.linspace(0,len(c)-1,40,dtype=int)],c.loc[c.POS.isin([5998724,22917619,67736589,67842269])&c.allele_index.eq(1)]]).drop_duplicates(['source_index','allele_index']).sort_values(['source_index','allele_index'])
c['ID']=['chrX:'+str(int(r.POS))+':'+r.REF+':'+r.ALT for r in c.itertuples(index=False)]
# Preserve provider representations; do not silently pick one duplicate identity.
c['duplicate_identity']=c.ID.duplicated(keep=False);c.loc[c.duplicate_identity].to_csv(O/'duplicate_identities.private.tsv.gz',sep='\t',index=False);c=c.loc[~c.duplicate_identity].reset_index(drop=True)
ps=pd.read_csv(N/'inputs/provider_X_aligned.psam',sep='\t',dtype={'#FID':str,'IID':str});ph=pd.read_csv(N/'inputs/phenotype.tsv',sep='\t',dtype={'FID':str,'IID':str}).set_index('IID')
idx=np.flatnonzero(ps.IID.isin(ph.index)).astype(np.uint32);ps=ps.iloc[idx].reset_index(drop=True);assert len(ps)==len(ph)==386069 and ps.SEX.isin([1,2]).all()
y=ph.psoriasis_case.reindex(ps.IID).to_numpy(int);male=ps.SEX.eq(1).to_numpy();assert int(y.sum())==5705
reader=pgenlib.PgenReader(str(X/'private_provider/acaf_threshold.chrX.pgen').encode(),sample_subset=idx,allele_idx_offsets=np.load(X/'source_allele_offsets.npy'))
buf=np.empty(len(idx),np.int8);gt=np.empty(2*len(idx),np.int32);case=y==1;ctrl=y==0;parts=[];start=time.time()
try:
 for k in range(0,len(c),5000):
  folder=O/('chunk_'+str(k//5000).zfill(4));folder.mkdir(exist_ok=True);done=folder/'COMPLETE.json'
  if done.exists():parts.append(json.loads(done.read_text()));continue
  assert shutil.disk_usage(N).free>128*1024**3,'scratch floor reached'
  chunk=c.iloc[k:k+5000];ps.to_csv(folder/'input.psam',sep='\t',index=False)
  pv=chunk[['CHROM','POS','ID','REF','ALT']].copy();pv.columns=['#CHROM','POS','ID','REF','ALT'];pv.to_csv(folder/'input.pvar',sep='\t',index=False)
  writer=pgenlib.PgenWriter(str(folder/'input.pgen').encode(),sample_ct=len(idx),variant_ct=len(chunk),nonref_flags=False);rows=[];sample_checks=[]
  try:
   for j,r in enumerate(chunk.itertuples(index=False)):
    reader.read(int(r.source_index),buf,allele_idx=int(r.allele_index));nonpar=2781479<int(r.POS)<155701383
    masked=0
    if nonpar:
     reader.read_alleles(int(r.source_index),gt);alleles=gt.reshape(-1,2);bad=male & (alleles[:,0]!=alleles[:,1]);masked=int((bad & (buf>=0)).sum());buf[bad]=-9
     assert np.isin(buf[male & (buf>=0)],[0,2]).all()
    ok=buf>=0;ploidy=np.where(male,1,2) if nonpar else np.full(len(idx),2);dos=buf/np.where(male,2.,1.) if nonpar else buf.astype(float)
    vals={}
    for lab,grp in [('case',case),('control',ctrl)]:
     valid=grp & ok;ac=int(dos[valid].sum());an=int(ploidy[valid].sum());vals.update({lab+'_AC':ac,lab+'_AN':an,lab+'_MAC':min(ac,an-ac),lab+'_called':int(valid.sum()),lab+'_callrate':float(valid.sum()/grp.sum())})
    ac=vals['case_AC']+vals['control_AC'];an=vals['case_AN']+vals['control_AN'];mac=min(ac,an-ac);maf=mac/an if an else 0.
    passed=vals['case_MAC']>=20 and vals['control_MAC']>=20 and mac>=100 and maf>=.01 and vals['case_callrate']>=.98 and vals['control_callrate']>=.98
    rows.append({'ID':r.ID,'source_index':int(r.source_index),'allele_index':int(r.allele_index),'nonPAR':nonpar,'male_heterozygous_calls_masked':masked,**vals,'total_MAC':mac,'MAF':maf,'count_call_maf_PASS':bool(passed)})
    if a.pilot:sample_checks.append(buf.copy())
    writer.append_biallelic(buf)
    if j%1000==0:print('X_SOURCE_PROJECT',k+j,len(c),flush=True)
  finally:writer.close()
  qc=pd.DataFrame(rows);qc.to_csv(folder/'variant_QC.private.tsv.gz',sep='\t',index=False)
  verify=pgenlib.PgenReader(str(folder/'input.pgen').encode());assert verify.get_variant_ct()==len(chunk) and verify.get_raw_sample_ct()==len(idx)
  if a.pilot:
   check=np.empty(len(idx),np.int8)
   for j,original in enumerate(sample_checks):verify.read(j,check);assert np.array_equal(check,original)
  verify.close()
  rec={'status':'projected_and_count_QC_complete','variants':len(chunk),'passed_count_call_maf':int(qc.count_call_maf_PASS.sum()),'n':len(idx),'pilot_roundtrip_verified':a.pilot,'scope':'site FILTER and female-control ancestry HWE still required; no association tests here','elapsed_seconds':time.time()-start}
  done.write_text(json.dumps(rec,indent=2));parts.append(rec)
 (O/'COMPLETE.json').write_text(json.dumps({'status':'projection_complete','candidate_variants':len(c),'chunks':len(parts),'passed_count_call_maf':sum(z['passed_count_call_maf'] for z in parts),'pilot':a.pilot,'site_and_HWE_QC_still_required':True},indent=2))
except Exception as e:
 (O/'FAILED.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc()}));raise
finally:reader.close()
