"""Bounded PLINK feasibility benchmark; not a chromosome-wide discovery scan."""
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json,time,subprocess,datetime
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
R=Path('/home/jupyter/x_followup_20261013');O=R/'chen_scan_benchmark';O.mkdir(exist_ok=True)
P='/opt/workbench-tools/binaries/bin/plink2'
def read(p):
    d=pd.read_csv(p,sep=r'\s+',dtype={'FID':str,'#FID':str,'IID':str});d.rename(columns={'#FID':'FID'},inplace=True);return d
c=read(B/'inputs/covariates.tsv');ph=read(B/'inputs/phenotype.tsv');s=read(B/'inputs/provider_X_aligned.psam')
k=pd.read_csv(B/'reconciliation_20260912/joint_G_sex_20260914/provider_pruned.keep.tsv',sep=r'\s+',header=None,names=['FID','IID'],dtype=str)
k=k[k.FID.str.lstrip('#')!='FID'];assert not k.duplicated(['FID','IID']).any()
d=c.merge(ph,on=['FID','IID'],validate='one_to_one').merge(s[['FID','IID','SEX']],on=['FID','IID'],validate='one_to_one').merge(k,on=['FID','IID'],validate='one_to_one')
base=['age','age2','age_sex']+[f'PC{i}' for i in range(1,11)]
assert set(d['sex'].unique())=={0,1} and set(d['psoriasis_case'].unique())=={0,1}
assert np.array_equal(d['sex'].to_numpy(),(d['SEX']==1).astype(int).to_numpy())
assert d[base].notna().all().all()
S=d['sex'].to_numpy(float);V=d[base].to_numpy(float)
V=np.column_stack([V,d[[f'PC{i}' for i in range(1,11)]+['age2']].to_numpy(float)*S[:,None]])
# Stable change of nuisance basis; retain sex as first covariate for one interaction.
A=np.column_stack([np.ones(len(S)),S]);coef=np.linalg.lstsq(A,V,rcond=None)[0];res=V-A@coef
Q,T=np.linalg.qr(res,mode='reduced');assert Q.shape[1]==24
assert np.min(np.abs(np.diag(T)))>1e-7
assert np.max(np.abs(res-Q@T))<1e-6
C=d[['FID','IID']].copy();C['S']=S
for i in range(24):C[f'Q{i+1}']=Q[:,i]*np.sqrt(len(S))
C.to_csv(O/'covariates.private.tsv',sep='\t',index=False)
F=d[['FID','IID']].copy();F['CASE']=d['psoriasis_case']+1;F.to_csv(O/'phenotype.private.tsv',sep='\t',index=False)
d[['FID','IID']].to_csv(O/'keep.private.tsv',sep='\t',index=False,header=False)
targets={'chrX:5995921:G:A','chrX:5996494:GCT:G','chrX:5998724:A:G'};chosen=None
for p in sorted((B/'X_projected_source').glob('chunk_*/input.pvar')):
    ids=set(line.split()[2] for line in p.open() if not line.startswith('#'))
    if targets<=ids:chosen=p.parent;break
assert chosen is not None
passed=set((B/'results/X_source'/chosen.name/'passed_variants.txt').read_text().split());assert targets<=passed
extra=sorted(passed-targets,key=lambda x:hashlib.sha256(('X-omnibus-benchmark-v1|'+x).encode()).digest())[:97]
(O/'extract100.txt').write_text('\n'.join(sorted(targets|set(extra)))+'\n')
(O/'targets.txt').write_text('\n'.join(sorted(targets))+'\n')
start={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'100-variant feasibility benchmark; not null calibration or a full scan','cohort':'Provider kinship>0.1 pruned analysis keep; complete original covariates','nuisance':'Same column span as target expanded nuisance model, using stable QR coordinates, with sex retained first','method':'PLINK double-precision ordinary logistic, joint Wald test of ADD + DOMDEV + ADDxS; not the earlier likelihood-ratio statistic','parameters':'1-28 selects ADD, DOMDEV, 25 covariates, ADDxS; tests 1,2,28; no automatic extra sex term','source_chunk':chosen.name,'n':len(d),'cases':int(d.psoriasis_case.sum()),'correction':'No discovery claims; benchmark and implementation validation only','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(O/'BENCHMARK_PROTOCOL.json').write_text(json.dumps(start,indent=2)+'\n')
for coding,extract in [(2,'extract100.txt'),(1,'targets.txt')]:
    prefix=O/f'coding{coding}'
    cmd=[P,'--pfile',str(chosen/'input'),'--keep',str(O/'keep.private.tsv'),'--covar',str(O/'covariates.private.tsv'),'--pheno',str(O/'phenotype.private.tsv'),'--extract',str(O/extract),'--glm','genotypic','interaction','no-x-sex','no-firth','--parameters','1-28','--tests','1,2,28','--xchr-model',str(coding),'--threads','8','--memory','16000','--out',str(prefix)]
    (O/f'coding{coding}.command.json').write_text(json.dumps(cmd));t=time.time()
    with (O/f'coding{coding}.console.log').open('w') as f:status=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=900).returncode
    receipt={'coding':coding,'returncode':status,'elapsed_seconds':time.time()-t};(O/f'coding{coding}.run.json').write_text(json.dumps(receipt));print(json.dumps(receipt),flush=True)
    if status:raise RuntimeError('PLINK feasibility run failed; inspect retained logs')
print('BENCHMARK_COMPLETED; output still requires term and numerical validation',flush=True)
