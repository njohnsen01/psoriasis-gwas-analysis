from pathlib import Path
import time,json,pandas as pd,numpy as np
O=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912/joint_G_sex_20260914/provider_quality');B=O.parents[2]
while not (O/'hail_status.json').exists():time.sleep(30)
s=json.loads((O/'hail_status.json').read_text());assert s['status']=='complete','Provider extraction failed; inspect restricted logs'
e=pd.read_csv(O/'three_sites_entries.private.tsv.bgz',sep='\t',compression='gzip',dtype={'s':str,'GT':str});k=pd.read_csv(B/'inputs/keep.tsv',sep=r'\s+',header=None,dtype=str);e=e[e.s.isin(k[1])]
a=pd.read_csv(B.parent/'sample_lead_qc_20260910/sample_audit.private.tsv',sep='\t',dtype={'IID':str});e=e.merge(a[['IID','SEX','psoriasis_case','ancestry_pred','site_id']],left_on='s',right_on='IID',validate='many_to_one');rows=[]
def category(row):
 gt=str(row.GT).replace('|','/').split('/')
 if not all(x.isdigit() for x in gt):return 'missing'
 allele='A' if row.locus=='chrX:5995921' else 'G'
 ix=json.loads(row.alleles).index(allele)
 count=sum(int(x)==ix for x in gt)
 return 'target_carrier' if count else 'no_target_allele'
e['target_category']=e.apply(category,axis=1)
for keys,z in e.groupby(['locus','ancestry_pred','SEX','psoriasis_case','target_category'],dropna=False):
 r=dict(zip(['locus','ancestry','sex','case_status','target_category'],keys));r['N']=len(z)
 for c in ['GQ','RGQ']:
  v=pd.to_numeric(z[c],errors='coerce').dropna();r[c+'_available']=len(v);r[c+'_median']=float(v.median()) if len(v) else None;r[c+'_below20']=int(v.lt(20).sum())
 r['AD_note']='Summed allelic depths are not complete DP';r['mapping_quality']='not in extracted fields';rows.append(r)
pd.DataFrame(rows).to_csv(O/'quality_by_carrier.private.tsv',sep='\t',index=False);(O/'CARRIER_SUMMARY_COMPLETE.json').write_text(json.dumps({'complete':True,'no_mapping_quality_claim':True,'counts_remain_private':True}));print('QUALITY_SUMMARY_COMPLETE')
