from pathlib import Path
import json,subprocess
Y=Path(__file__).resolve().parents[1];N=Y.parent;O=Y/'counts';O.mkdir(exist_ok=True)
for g in ['cases','controls']:
 cmd=['plink2','--pgen',str(Y/'source/acaf_threshold.chrY.pgen'),'--pvar',str(Y/'source/acaf_threshold.chrY.pvar'),'--psam',str(Y/'aligned_provider.psam'),'--keep',str(N/'stratified_step1/male'/(g+'.autosomal.keep.tsv')),'--keep-males','--freq','counts','--missing','variant-only','--threads','4','--memory','8192','--out',str(O/g)]
 (O/(g+'.command.json')).write_text(json.dumps(cmd,indent=2))
 with (O/(g+'.console.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,(g,r.returncode)
 print('Y_SOURCE_COUNTS_COMPLETE',g,flush=True)
(O/'COMPLETE.json').write_text(json.dumps({'status':'complete','scope':'preliminary native source counts; exact masked count/call/MAF and site gates remain'}))
