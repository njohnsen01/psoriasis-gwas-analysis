from pathlib import Path
import json,subprocess,hashlib,shutil
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;B=O.parent;N=B.parent;S=N/'stratified_step1/female';G=N/'X_projected_source/chunk_0009/input';T=B/'joint_G_sex_20260914/targets.txt';exe=shutil.which('regenie')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert (O/'MODEL_INPUT_VERIFIED.json').exists();assert not (O/'REGENIE_PLAN_v2.json').exists();meta=json.loads((S/'VALIDATED.json').read_text());assert sha(exe)==meta['regenie_binary_sha256']
targets=T.read_text().splitlines();prior=pd.read_csv(N/'results/X_stratified/female/chunk_0009/male02_psoriasis_case.regenie',sep=r'\s+');prior=prior[prior.ID.isin(targets)];assert set(prior.ID)==set(targets)
commands={}
for name,keep,approx in [('full_approx',S/'keep.tsv',True),('pruned_approx',O/'pruned_female.keep.private.tsv',True),('pruned_exact',O/'pruned_female.keep.private.tsv',False)]:
 cmd=[exe,'--step','2','--bt','--pgen',str(G),'--keep',str(keep),'--extract',str(T),'--phenoFile',str(S/'phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(S/'covariates.tsv'),'--bsize','200','--firth','--pThresh','.01','--pred',str(S/'model_pred.list'),'--par-region','hg38','--threads','2','--out',str(O/('regenie_v2_'+name))]
 if approx:cmd.append('--approx')
 commands[name]=cmd
(O/'REGENIE_PLAN_v2.json').write_text(json.dumps(dict(commands=commands,scope='Only three fixed variants. Fixed original female predictor/covariates/phenotype. Full versus pruned sample; approximate versus exact Firth on pruned women.',hashes={str(p):sha(p) for p in [Path(exe),S/'model_1.loco',S/'covariates.tsv',S/'phenotype.tsv',S/'keep.tsv',O/'pruned_female.keep.private.tsv',T,Path(__file__)]}),indent=2))
results=[]
for name,cmd in commands.items():
 print('START',name,flush=True)
 with (O/('regenie_v2_'+name+'.console.private.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=1800)
 f=O/('regenie_v2_'+name+'_psoriasis_case.regenie');r=pd.read_csv(f,sep=r'\s+');assert set(r.ID)==set(targets) and len(r)==3 and r.TEST.eq('ADD').all();assert np.isfinite(r[['BETA','SE','LOG10P']].to_numpy()).all()
 if name=='full_approx':
  x=r.set_index('ID').loc[targets];y=prior.set_index('ID').loc[targets];assert np.allclose(x[['BETA','SE','LOG10P']],y[['BETA','SE','LOG10P']],rtol=1e-5,atol=1e-5)
 for x in r.itertuples():results.append(dict(target='g'+str(targets.index(x.ID)),model=name,beta=x.BETA,se=x.SE,OR=float(np.exp(x.BETA)),lower=float(np.exp(x.BETA-1.95996398454*x.SE)),upper=float(np.exp(x.BETA+1.95996398454*x.SE)),p=float(10**(-x.LOG10P)),verified=True))
 pd.DataFrame(results).to_csv(O/'MATCHED_REGENIE_RESULTS.csv',index=False);print('VERIFIED',name,flush=True)
(O/'MATCHED_REGENIE_VERIFIED').write_text('All nine target estimates verified; original full-female scan reproduced.')
