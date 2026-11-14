"""Synthetic-only compatibility check. No research participants or results."""
from pathlib import Path
import numpy as np,pandas as pd,subprocess,json,os,datetime
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'synthetic_regenie_smoke';O.mkdir(exist_ok=True);rng=np.random.default_rng(9031)
n=2000;s=np.r_[np.full(1000,2),np.full(1000,1)];g=np.empty((8,n),dtype=np.uint8)
for j in range(8):g[j,:1000]=rng.binomial(2,.008,1000);g[j,1000:]=2*rng.binomial(1,.008,1000)
y=rng.binomial(1,np.where(g[:4].max(axis=0)>0,.6,.07));ids=[f'SYN{i:05d}' for i in range(n)]
pd.DataFrame({'FID':ids,'IID':ids,'psoriasis_case':y}).to_csv(O/'phenotype.tsv',sep='\t',index=False)
age=rng.normal(size=n);pd.DataFrame({'FID':ids,'IID':ids,'age':age,'age2':age**2,'PC1':rng.normal(size=n)}).to_csv(O/'covariates.tsv',sep='\t',index=False)
(O/'input.fam').write_text(''.join(f'{v} {v} 0 0 {s[i]} -9\n' for i,v in enumerate(ids)))
vid=[f'SYN_X_{i}' for i in range(8)];(O/'input.bim').write_text(''.join(f'23 {v} 0 {3000000+i} C A\n' for i,v in enumerate(vid)))
with (O/'input.bed').open('wb') as f:
 f.write(b'\x6c\x1b\x01')
 for j in range(8):
  code=np.choose(g[j],[3,2,0]).astype(np.uint8).reshape(-1,4);f.write(np.sum(code<<np.array([0,2,4,6],dtype=np.uint8),axis=1).astype(np.uint8).tobytes())
(O/'annotations.txt').write_text(''.join(f'{v} SYNGENE '+('PTV' if i<4 else 'DAMAGING_MISSENSE')+'\n' for i,v in enumerate(vid)))
(O/'setlist.txt').write_text('SYNGENE 23 3000000 '+','.join(vid)+'\n');(O/'masks.txt').write_text('PTV PTV\nPTV_DAMAGING PTV,DAMAGING_MISSENSE\n');(O/'aaf.txt').write_text(''.join(f'{v} 0.008\n' for v in vid))
with (O/'null.loco').open('w') as f:
 f.write('FID_IID '+' '.join(v+'_'+v for v in ids)+'\n')
 for ch in range(1,24):f.write(str(ch)+' '+' '.join('0' for _ in ids)+'\n')
(O/'pred.list').write_text('psoriasis_case '+str(O/'null.loco')+'\n')
for sex,code in [('female',2),('male',1)]:
 (O/f'{sex}.keep').write_text(''.join(f'{v} {v}\n' for i,v in enumerate(ids) if s[i]==code))
 cmd=['/opt/workbench-tools/binaries/bin/regenie','--step','2','--bed',str(O/'input'),'--keep',str(O/f'{sex}.keep'),'--phenoFile',str(O/'phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(O/'covariates.tsv'),'--pred',str(O/'pred.list'),'--bt','--bsize','200','--firth','--approx','--pThresh','0.05','--minMAC','1','--anno-file',str(O/'annotations.txt'),'--set-list',str(O/'setlist.txt'),'--mask-def',str(O/'masks.txt'),'--aaf-file',str(O/'aaf.txt'),'--aaf-bins','0.01','--vc-tests','skato','--vc-maxAAF','0.01','--vc-MACthr','10','--singleton-carrier','--check-burden-files','--strict-check-burden','--write-mask','--write-mask-snplist','--par-region','hg38','--threads','1','--out',str(O/sex)]
 with (O/f'{sex}.console.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
 print(sex,'returncode',p.returncode,flush=True);assert p.returncode==0
 print((O/f'{sex}_psoriasis_case.regenie').read_text(),flush=True)
(O/'SMOKE_COMPLETE.json').write_text(json.dumps({'status':'synthetic_native_compatibility_check_passed','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'warning':'These are simulated data with artificially enriched cases; not study results or calibration validation'},indent=2)+'\n')
