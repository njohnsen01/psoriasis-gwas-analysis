from pathlib import Path
import json,hashlib,subprocess
import numpy as np,pandas as pd
from scipy.stats import fisher_exact
R=Path(__file__).resolve().parent
D=pd.read_csv(R/'sample_audit.private.tsv',sep='	');G=np.load(R/'lead_G.private.npy');H=np.load(R/'lead_H.private.npy')
leads=[('rs183018930','chrX:5998724:A:G','G'),('rs5970952','chrX:22917619:C:T','T'),('rs5964607','chrX:67736589:T:C','C'),('rs5919457','chrX:67842269:G:T','T')]
assert G.shape==H.shape==(4,len(D))
rows=[];public=[];den=[]
for ancestry in ['EUR','AFR','AMR']:
 for sex,label in [(2,'Female'),(1,'Male')]:
  ix=np.flatnonzero(D.strict_eligible & D.ancestry_pred.str.upper().eq(ancestry) & D.SEX.eq(sex))
  y=D.iloc[ix].psoriasis_case.to_numpy();assert set(y)=={0,1}
  den.append({'ancestry':ancestry,'sex':label,'cases':int(y.sum()),'controls':int((y==0).sum())})
  for j,(rs,key,ea) in enumerate(leads):
   g=G[j,ix];het=H[j,ix];ok=(g>=0)&(~het if sex==1 else True);ploidy=1 if sex==1 else 2
   ac=g/(2 if sex==1 else 1);missing=[];parts=[];pr={'lead':rs,'ancestry':ancestry,'sex':label,'EA':ea}
   for status,v in [('case',1),('control',0)]:
    m=y==v;n=int(m.sum());nc=int((m&ok).sum());nm=n-nc;a=int(ac[m&ok].sum());an=nc*ploidy
    carrier=int((m&ok&(ac>0)).sum());noncarrier=int((m&ok&(ac==0)).sum());nonEAcarrier=int((m&ok&(ac<ploidy)).sum())
    parts.append({'n':n,'called':nc,'missing':nm,'AC':a,'AN':an,'carriers':carrier,'nonEAcarriers':nonEAcarrier,'call_rate':nc/n})
    small=lambda z:0<z<20
    hide=small(nm) or any(small(z) for z in [a,an-a,carrier,nonEAcarrier])
    ac_hidden=any(small(z) for z in [a,an-a,carrier,nonEAcarrier])
    pr[status+'_AC_AN']=('S' if ac_hidden else str(a))+'/'+('S' if small(nm) else str(an))
    pr[status+'_call_pct']=('≥98' if nc/n>=.98 else '<98') if small(nm) else format(100*nc/n,'.2f')
    missing.append(nm)
   p=float(fisher_exact([[parts[0]['missing'],parts[0]['called']],[parts[1]['missing'],parts[1]['called']]]).pvalue)
   pr['missingness_P']='S' if any(0<n<20 for n in missing) else format(p,'.3g')
   rows.append({'lead':rs,'key':key,'EA':ea,'ancestry':ancestry,'sex':label,'case':parts[0],'control':parts[1],'missingness_p':p})
   public.append(pr)
(R/'lead_qc.private.json').write_text(json.dumps(rows,indent=2));(R/'lead_qc.public.json').write_text(json.dumps(public,indent=2));(R/'lead_denominators.public.json').write_text(json.dumps(den,indent=2))
D.loc[D.strict_eligible,['#FID','IID']].to_csv(R/'strict.keep.private.tsv',sep='	',index=False,header=False)
for ancestry in ['EUR','AFR','AMR']:
 d=D[D.strict_eligible & D.SEX.eq(2) & D.psoriasis_case.eq(0) & D.ancestry_pred.str.upper().eq(ancestry)]
 keep=R/(ancestry+'.female_controls.keep.private.tsv');d[['#FID','IID']].to_csv(keep,sep='	',index=False,header=False)
 cmd=['plink2','--pfile',str(R.parent/'X_source_reanalysis_20260910/private_provider_leads'),'--keep',str(keep),'--hardy','--out',str(R/(ancestry+'.female_controls'))]
 (R/(ancestry+'.hwe.command.json')).write_text(json.dumps(cmd));res=subprocess.run(cmd,capture_output=True,text=True);assert res.returncode==0,res.stdout+res.stderr
print(json.dumps({'status':'complete','rows':len(rows),'min_call_rate':min(z[s]['call_rate'] for z in rows for s in ['case','control']),'min_missingness_p':min(z['missingness_p'] for z in rows),'diagnostic_tests':len(rows)},indent=2))
