from pathlib import Path
import json,re
import pandas as pd,numpy as np
R=Path(__file__).resolve().parent
D=pd.read_csv(R/'sample_audit.private.tsv',sep='\t');G=np.load(R/'lead_G.private.npy');H=np.load(R/'lead_H.private.npy');D['gindex']=np.arange(len(D))
E=pd.read_csv(R/'four_sites_entries.private.tsv.bgz',sep='\t',compression='gzip',dtype={'s':int,'GT':str})
E=E.merge(D[['IID','SEX','psoriasis_case','ancestry_pred','strict_eligible','gindex']],left_on='s',right_on='IID',validate='many_to_one');E=E[E.strict_eligible].copy()
leads=[('rs183018930',5998724,'G'),('rs5970952',22917619,'T'),('rs5964607',67736589,'C'),('rs5919457',67842269,'T')]
rows=[];validation=[]
for j,(rs,pos,ea) in enumerate(leads):
 e=E[E.locus.eq('chrX:'+str(pos))].copy();alleles=json.loads(e.alleles.iloc[0]);idx=alleles.index(ea)
 calls=e.GT.str.replace('|','/',regex=False).str.split('/').map(lambda x:list(map(int,x)))
 e['expected']=calls.map(lambda x:sum(v==idx for v in x)*(2 if len(x)==1 else 1));e['expected_het']=calls.map(lambda x:len(x)==2 and x[0]!=x[1])
 e['pg']=G[j,e.gindex.to_numpy(int)];e['ph']=H[j,e.gindex.to_numpy(int)]
 assert (e.expected==e.pg).all(),rs+' dosage mismatch';assert (e.expected_het==e.ph).all(),rs+' heterozygote mismatch'
 present=set(e.IID);absent=D[D.strict_eligible & ~D.IID.isin(present)];assert (G[j,absent.gindex.to_numpy(int)]<0).all(),rs+' missingness mismatch'
 e['ADsum']=e.AD.map(lambda x:sum(json.loads(x)) if isinstance(x,str) else np.nan)
 e['valid']=(e.pg>=0)&~((e.SEX==1)&e.ph)
 validation.append({'lead':rs,'source_ref':alleles[0],'effect_allele':ea,'genotypes_match_PGEN':True,'missingness_match_PGEN':True})
 for ancestry in ['EUR','AFR','AMR']:
  for sex,label in [(2,'Female'),(1,'Male')]:
   for y,status in [(1,'Case'),(0,'Control')]:
    x=e[e.ancestry_pred.str.upper().eq(ancestry)&e.SEX.eq(sex)&e.psoriasis_case.eq(y)&e.valid];eligible=int((D.strict_eligible & D.ancestry_pred.str.upper().eq(ancestry) & D.SEX.eq(sex) & D.psoriasis_case.eq(y)).sum())
    out={'lead':rs,'ancestry':ancestry,'sex':label,'status':status,'eligible_N':eligible,'called_N':len(x) if eligible-len(x)==0 or eligible-len(x)>=20 else 'S'}
    for metric in ['GQ','RGQ','ADsum']:
     a=pd.to_numeric(x[metric],errors='coerce').dropna();n=len(a);missing=len(x)-n
     out[metric+'_n']=n if (n==0 or n>=20) and (missing==0 or missing>=20) else 'S'
     qs=a.quantile([.25,.5,.75]).tolist() if n>=20 else [None]*3
     for name,value in zip(['q25','median','q75'],qs):out[metric+'_'+name]=None if value is None else round(float(value),2)
    if out['called_N']=='S':
     for metric in ['GQ','RGQ','ADsum']:out[metric+'_n']='S'
    out['DP']='not present';out['MQ']='not present';rows.append(out)
(R/'quality_by_stratum.public.json').write_text(json.dumps(rows,indent=2));(R/'quality_validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps({'status':'validated','stratum_rows':len(rows),'genotypes_match':True,'source_missingness_match':True,'GQ_median_range':[min(r['GQ_median'] for r in rows if r['GQ_median'] is not None),max(r['GQ_median'] for r in rows if r['GQ_median'] is not None)],'ADsum_median_range':[min(r['ADsum_median'] for r in rows if r['ADsum_median'] is not None),max(r['ADsum_median'] for r in rows if r['ADsum_median'] is not None)]},indent=2))
