from pathlib import Path
import argparse,csv,json,sys,re,traceback,time,hashlib,subprocess
import numpy as np,pandas as pd
Y=Path(__file__).resolve().parents[1];N=Y.parent
sys.path.insert(0,str(N.parent/'X_source_reanalysis_20260910/python_dependencies'))
import pgenlib
pa=argparse.ArgumentParser();pa.add_argument('--pilot',action='store_true');args=pa.parse_args()
O=Y/('projection_pilot' if args.pilot else 'projected');O.mkdir(exist_ok=True)
def save(p,v):p.write_text(json.dumps(v,indent=2)+chr(10))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 assert (Y/'counts/COMPLETE.json').exists() and (Y/'site_filters/COMPLETE.json').exists()
 p=Y/'source/acaf_threshold.chrY.pvar'
 with p.open() as f:meta=sum(1 for l in f if l.startswith('##'))
 pv=pd.read_csv(p,sep=chr(9),skiprows=meta,dtype={'#CHROM':str,'ID':str,'REF':str,'ALT':str})
 sites=pd.read_csv(Y/'site_filters/chrY_site_filters.tsv.bgz',sep=chr(9),compression='gzip')
 assert sites.is_y_nonpar.dtype.kind=='b'
 site_map={}
 for r in sites.itertuples(index=False):
  alleles=json.loads(r.alleles);key=(str(r.locus),alleles[0],tuple(sorted(alleles[1:])))
  assert key not in site_map,('duplicate_Hail_site_key',key)
  site_map[key]=(len(json.loads(r.filters))==0,bool(r.is_y_nonpar),alleles)
 tables=[pd.read_csv(Y/'counts'/f,sep=chr(9),dtype={'ID':str}) for f in ['cases.acount','controls.acount','cases.vmiss','controls.vmiss']]
 assert all(len(t)==len(pv) for t in tables)
 candidates=[];mapping=[]
 for i,(pr,ca,co,cm,om) in enumerate(zip(pv.to_dict('records'),*[t.to_dict('records') for t in tables],strict=True)):
  assert pr['ID']==ca['ID']==co['ID']==cm['ID']==om['ID']
  assert pr['REF']==ca['REF']==co['REF'] and pr['ALT']==ca['ALT']==co['ALT']
  key=(pr['#CHROM']+':'+str(pr['POS']),pr['REF'],tuple(sorted(pr['ALT'].split(','))))
  assert key in site_map,('unmatched_source_site',i)
  passed,nonpar,ha=site_map[key];mapping.append({'source_index':i,'site_FILTER_PASS':passed,'is_y_nonpar':nonpar,'ALT_order_differs':ha[1:]!=pr['ALT'].split(',')})
  if not passed or not nonpar:continue
  ac=list(map(float,str(ca['ALT_CTS']).split(',')));bc=list(map(float,str(co['ALT_CTS']).split(',')));alts=pr['ALT'].split(',')
  assert len(ac)==len(bc)==len(alts)
  an=float(ca['OBS_CT']);bn=float(co['OBS_CT'])
  cr=1-float(cm['MISSING_CT'])/float(cm['OBS_CT']);br=1-float(om['MISSING_CT'])/float(om['OBS_CT'])
  if cr<.98 or br<.98:continue
  for ai,(alt,a,b) in enumerate(zip(alts,ac,bc),1):
   if not re.fullmatch('[ACGT]+',pr['REF']) or not re.fullmatch('[ACGT]+',alt):continue
   if min(a,an-a)<20 or min(b,bn-b)<20 or min(a+b,an+bn-a-b)<100:continue
   candidates.append({'source_index':i,'allele_index':ai,'CHROM':pr['#CHROM'],'POS':int(pr['POS']),'REF':pr['REF'],'ALT':alt,'ID':'chrY:'+str(pr['POS'])+':'+pr['REF']+':'+alt})
 pd.DataFrame(mapping).to_csv(O/'source_site_mapping.private.tsv.gz',sep=chr(9),index=False)
 c=pd.DataFrame(candidates);assert len(c)>0,'No count/site-qualified candidates; review before model stage'
 duplicates=c.ID.duplicated(keep=False);c.loc[duplicates].to_csv(O/'duplicate_identities.private.tsv.gz',sep=chr(9),index=False);c=c.loc[~duplicates].reset_index(drop=True)
 total_candidates=len(c)
 if args.pilot:c=c.iloc[np.unique(np.concatenate([np.arange(min(30,len(c))),np.linspace(0,len(c)-1,min(30,len(c)),dtype=int)]))].reset_index(drop=True)
 c.to_csv(O/'candidate_alleles.private.tsv.gz',sep=chr(9),index=False)
 ps=pd.read_csv(Y/'aligned_provider.psam',sep=chr(9),dtype={'#FID':str,'IID':str})
 ph=pd.read_csv(N/'stratified_step1/male/phenotype.tsv',sep=chr(9),dtype={'FID':str,'IID':str}).set_index('IID')
 idx=np.flatnonzero(ps.IID.isin(ph.index)).astype(np.uint32);ps=ps.iloc[idx].reset_index(drop=True)
 assert len(idx)==len(ph)==151282 and ps.SEX.eq(1).all()
 case=ph.psoriasis_case.reindex(ps.IID).eq(1).to_numpy();ctrl=~case;assert int(case.sum())==2259
 offsets=np.concatenate([[0],np.cumsum(1+pv.ALT.str.count(',').to_numpy()+1)]).astype(np.uintp)
 reader=pgenlib.PgenReader(str(Y/'source/acaf_threshold.chrY.pgen').encode(),sample_subset=idx,allele_idx_offsets=offsets)
 assert reader.get_variant_ct()==len(pv)
 ps.to_csv(O/'input.psam',sep=chr(9),index=False)
 c[['CHROM','POS','ID','REF','ALT']].rename(columns={'CHROM':'#CHROM'}).to_csv(O/'input.pvar',sep=chr(9),index=False)
 writer=pgenlib.PgenWriter(str(O/'input.pgen').encode(),sample_ct=len(idx),variant_ct=len(c),nonref_flags=False)
 buf=np.empty(len(idx),np.int8);gt=np.empty(2*len(idx),np.int32);rows=[];checks=[]
 try:
  for k,r in enumerate(c.itertuples(index=False)):
   reader.read(int(r.source_index),buf,allele_idx=int(r.allele_index));reader.read_alleles(int(r.source_index),gt)
   alleles=gt.reshape(-1,2);bad=(alleles[:,0]!=alleles[:,1])|(alleles[:,0]<0)|(alleles[:,1]<0)
   masked=int((bad & (buf>=0)).sum());buf[bad]=-9
   assert np.isin(buf[buf>=0],[0,2]).all()
   ok=buf>=0;dos=buf/2.;vals={}
   for lab,grp in [('case',case),('control',ctrl)]:
    valid=grp&ok;ac=int(dos[valid].sum());an=int(valid.sum());vals.update({lab+'_AC':ac,lab+'_AN':an,lab+'_MAC':min(ac,an-ac),lab+'_called':an,lab+'_callrate':float(an/grp.sum())})
   mac=min(vals['case_AC']+vals['control_AC'],vals['case_AN']+vals['control_AN']-vals['case_AC']-vals['control_AC']);an=vals['case_AN']+vals['control_AN'];maf=mac/an if an else 0
   passed=vals['case_MAC']>=20 and vals['control_MAC']>=20 and mac>=100 and maf>=.01 and vals['case_callrate']>=.98 and vals['control_callrate']>=.98
   rows.append(dict(ID=r.ID,source_index=int(r.source_index),allele_index=int(r.allele_index),heterozygous_or_invalid_calls_masked=masked,**vals,total_MAC=mac,MAF=maf,PASS=bool(passed)))
   writer.append_biallelic(buf)
   if args.pilot:checks.append(buf.copy())
   if k%1000==0:print('Y_PROJECTED',k,len(c),flush=True)
 finally:writer.close();reader.close()
 q=pd.DataFrame(rows);q.to_csv(O/'final_variant_QC.private.tsv.gz',sep=chr(9),index=False);q.loc[q.PASS,'ID'].to_csv(O/'passed_variants.txt',header=False,index=False)
 verify=pgenlib.PgenReader(str(O/'input.pgen').encode());assert verify.get_variant_ct()==len(c) and verify.get_raw_sample_ct()==len(idx)
 if args.pilot:
  for k,a in enumerate(checks):verify.read(k,buf);assert np.array_equal(a,buf)
 verify.close()
 for lab,keep in [('case','cases'),('control','controls')]:
  prefix=O/('independent_'+lab)
  cmd=['plink2','--pfile',str(O/'input'),'--keep',str(N/'stratified_step1/male'/(keep+'.autosomal.keep.tsv')),'--freq','counts','--missing','variant-only','--threads','4','--memory','4096','--out',str(prefix)]
  save(O/(lab+'_independent.command.json'),cmd)
  with (O/(lab+'_independent.console.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
  assert r.returncode==0
  z=pd.read_csv(str(prefix)+'.acount',sep=chr(9));m=pd.read_csv(str(prefix)+'.vmiss',sep=chr(9));assert z.ID.tolist()==q.ID.tolist()==m.ID.tolist()
  assert np.array_equal(pd.to_numeric(z.ALT_CTS).to_numpy(),q[lab+'_AC'].to_numpy())
  assert np.array_equal(z.OBS_CT.to_numpy(),q[lab+'_AN'].to_numpy())
  assert np.array_equal((m.OBS_CT-m.MISSING_CT).to_numpy(),q[lab+'_called'].to_numpy())
 save(O/'COMPLETE.json',{'status':'validated','pilot':args.pilot,'source_sites':len(pv),'candidate_alleles_before_pilot':total_candidates,'projected_alleles':len(c),'passed':int(q.PASS.sum()),'n_males':len(idx),'independent_PLINK_count_and_missingness_checks':'all projected variants agree exactly','genotype_roundtrip':args.pilot,'inputs':{p.name:sha(p) for p in [O/'input.pgen',O/'input.pvar',O/'input.psam']},'script_sha256':sha(Path(__file__)),'site_metadata_sha256':sha(Y/'site_filters/chrY_site_filters.tsv.bgz'),'release':'HELD'})
 print('Y_PROJECTION_QC_COMPLETE',len(c),int(q.PASS.sum()),flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
