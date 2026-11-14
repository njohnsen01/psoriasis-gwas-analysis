from pathlib import Path
import json,pandas as pd,numpy as np,hashlib
O=Path(__file__).resolve().parent;P=O.parent;B=P.parents[2]
e=pd.read_csv(P/'three_sites_entries.private.tsv.bgz',sep='\t',compression='gzip',dtype={'s':str,'GT':str})
k=pd.read_csv(B/'inputs/keep.tsv',sep=r'\s+',header=None,dtype=str);e=e.loc[e.s.isin(k[1])]
a=pd.read_csv(B.parent/'sample_lead_qc_20260910/sample_audit.private.tsv',sep='\t',dtype={'IID':str})
e=e.merge(a[['IID','SEX','psoriasis_case','ancestry_pred']],left_on='s',right_on='IID',validate='many_to_one')
rows=[]
for r in e.itertuples():
 try:
  gt=[int(x) for x in str(r.GT).replace('|','/').split('/')]
  al=json.loads(r.alleles);idx=al.index('A' if r.locus=='chrX:5995921' else 'G')
  if idx not in gt:continue
  ad=json.loads(r.AD) if isinstance(r.AD,str) else r.AD
  if not isinstance(ad,list) or len(ad)!=len(al) or any(v is None or v<0 for v in ad) or sum(ad)<=0:continue
  kind='heterozygous carrier' if len(gt)==2 and len(set(gt))==2 else 'single-allele carrier'
  rows.append({'variant':r.locus,'sex':r.SEX,'ancestry':str(r.ancestry_pred).upper(),'case_status':r.psoriasis_case,'carrier_type':kind,'target_AD_fraction':ad[idx]/sum(ad),'summed_AD':sum(ad)})
 except (ValueError,TypeError,IndexError):continue
f=pd.DataFrame(rows);f.to_pickle(O/'ad_values.private.pkl')
private=[];public=[]
for keys,z in f.groupby(['variant','sex','ancestry','case_status','carrier_type'],dropna=False):
 row=dict(zip(['variant','sex','ancestry','case_status','carrier_type'],keys));n=len(z)
 pr={**row,'N_AD_available':n};pu={**row,'status':'reported' if n>=21 else 'withheld'}
 for metric in ['target_AD_fraction','summed_AD']:
  for q,label in [(0.25,'q25'),(0.5,'median'),(0.75,'q75')]:
   val=float(z[metric].quantile(q));pr[metric+'_'+label]=val;pu[metric+'_'+label]=val if n>=21 else None
 private.append(pr);public.append(pu)
pd.DataFrame(private).to_csv(O/'allele_depth_summary.private.tsv',sep='\t',index=False)
pd.DataFrame(public).to_csv(O/'allele_depth_summary_reviewed.tsv',sep='\t',index=False)
audit={'complete':True,'source':'provider acaf_threshold multiMT','groups_total':len(public),'group_counts_withheld':True,'min_release_group':21,'no_participant_records_in_reviewed_file':True,'DP_available':False,'MQ_available':False,'note':'Summed AD is not complete read depth. Descriptive assessment only; no participants excluded and no association refit. Sex codes follow existing sample audit; verify labels before manuscript use.'}
(O/'AD_AUDIT.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit))

