"""Match provider X sites to Hail site FILTER without assuming ALT order."""
from pathlib import Path
import csv,gzip,hashlib,json,os,time,traceback
import numpy as np,pandas as pd,psutil
N=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as handle:
  for b in iter(lambda:handle.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,d):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.partial');tmp.write_text(json.dumps(d,indent=2));tmp.replace(p)
def sitekey(locus,alleles):
 assert isinstance(locus,str) and locus.startswith('chrX:')
 assert isinstance(alleles,list) and len(alleles)>=2 and all(isinstance(a,str) and a for a in alleles)
 assert len(set(alleles))==len(alleles),'duplicate allele within site'
 # FILTER is site-level; preserve REF but allow non-reference allele reordering.
 return (locus,alleles[0],tuple(sorted(alleles[1:])))
def index_rows(rows):
 index={};duplicates=set();n=0
 for r in rows:
  n+=1;a=r['alleles'];fs=r['filters']
  a=json.loads(a) if isinstance(a,str) else a
  fs=json.loads(fs) if isinstance(fs,str) and fs not in ['NA',''] else fs
  assert fs is None or (isinstance(fs,list) and all(isinstance(x,str) for x in fs)),'unsupported FILTER encoding'
  k=sitekey(r['locus'],a)
  if k in index:duplicates.add(k)
  else:index[k]=(a,fs)
 return index,duplicates,n
def lookup(locus,alleles,index,duplicates):
 k=sitekey(locus,alleles)
 if k in duplicates:return {'site_mapping':'ambiguous_duplicate','site_FILTER_PASS':False,'site_filters':None,'ALT_order_differs':False}
 if k not in index:return {'site_mapping':'unmatched','site_FILTER_PASS':False,'site_filters':None,'ALT_order_differs':False}
 aa,fs=index[k]
 return {'site_mapping':'matched' if isinstance(fs,list) else 'FILTER_unavailable','site_FILTER_PASS':isinstance(fs,list) and len(fs)==0,'site_filters':fs,'ALT_order_differs':aa!=alleles}
def main():
 O=N/'source_X_sites';receipt=O/'COMPLETE.json';pid=json.loads((N/'status/X_sites_pid.json').read_text())['pid']
 while not receipt.exists():
  if not psutil.pid_exists(pid) or psutil.Process(pid).status()==psutil.STATUS_ZOMBIE:raise RuntimeError('Site FILTER extraction ended without completion')
  time.sleep(30)
 source=json.loads(receipt.read_text());assert source['status']=='site_filters_exported'
 path=O/'chrX_site_filters.tsv.bgz'
 with gzip.open(path,'rt') as handle:
  rd=csv.DictReader(handle,delimiter='\t');assert rd.fieldnames==['locus','alleles','filters'],rd.fieldnames
  index,duplicates,hn=index_rows(rd)
 P=N.parent/'X_source_reanalysis_20260910/private_provider/acaf_threshold.chrX.pvar'
 rows=[]
 with P.open() as handle:
  rd=csv.DictReader((line for line in handle if not line.startswith('##')),delimiter='\t')
  for i,r in enumerate(rd):
   alleles=[r['REF']]+r['ALT'].split(',');locus=r['#CHROM']+':'+r['POS']
   match=lookup(locus,alleles,index,duplicates)
   rows.append({'source_index':i,'CHROM':r['#CHROM'],'POS':int(r['POS']),'REF':r['REF'],'ALT':r['ALT'],**match})
 df=pd.DataFrame(rows);assert len(df)==2153825
 df['site_filters']=df.site_filters.map(lambda x:json.dumps(x,separators=(',',':')))
 df.to_csv(O/'provider_X_site_QC.tsv.gz',sep='\t',index=False)
 cc=pd.read_csv(N/'source_X_candidates/count_candidates.tsv.gz',sep='\t',usecols=['source_index']);eligible=df.loc[df.source_index.isin(cc.source_index)]
 summary={'status':'reconciled','hail_rows':hn,'provider_sites':len(df),'candidate_sites':len(eligible),'all_site_mapping_counts':df.site_mapping.value_counts().to_dict(),'candidate_site_mapping_counts':eligible.site_mapping.value_counts().to_dict(),'ALT_order_reconciled_sites':int(df.ALT_order_differs.sum()),'site_filter_source':source,'input_sha256':{'hail_export':sha(path),'provider_pvar':sha(P)},'output_sha256':sha(O/'provider_X_site_QC.tsv.gz'),'key_rule':'exact locus and REF and unordered complete ALT set; genotype allele indices stay in PGEN order; duplicate or unmatched sites never PASS'}
 save(O/'RECONCILIATION_COUNTS.json',summary)
 assert eligible.site_mapping.eq('matched').all(),'Candidate site mappings require manual review; inspect mapping counts and QC table'
 historical=json.loads((N.parent/'sample_lead_qc_20260910/hail_rows.json').read_text())
 for r in historical:
  v=df.loc[df.CHROM.eq(r['locus'].split(':')[0])&df.POS.eq(int(r['locus'].split(':')[1]))]
  matches=[t for t in v.itertuples(index=False) if sitekey(r['locus'],[t.REF]+t.ALT.split(','))==sitekey(r['locus'],r['alleles'])]
  assert len(matches)==1 and bool(matches[0].site_FILTER_PASS)==(r['filters']==[])
 save(O/'RECONCILIATION_COMPLETE.json',summary);print('PSORIASIS_X_SITE_RECONCILIATION_COMPLETE',summary['candidate_sites'],flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:
  save(N/'status/X_site_reconciliation_failed.json',{'error':str(e),'traceback':traceback.format_exc()});raise
