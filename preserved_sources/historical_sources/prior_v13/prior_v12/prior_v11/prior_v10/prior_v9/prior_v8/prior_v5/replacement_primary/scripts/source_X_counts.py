from pathlib import Path
import os,json,subprocess
N=Path(__file__).resolve().parents[1];X=N.parent/'X_source_reanalysis_20260910';O=N/'source_X_counts';O.mkdir(exist_ok=True)
for g in ['cases','controls']:
 cmd=['plink2','--pgen',str(X/'private_provider/acaf_threshold.chrX.pgen'),'--pvar',str(X/'private_provider/acaf_threshold.chrX.pvar'),'--psam',str(N/'inputs/provider_X_aligned.psam'),'--split-par','hg38','--keep',str(N/'inputs'/(g+'.keep.tsv')),'--freq','counts','--missing','variant-only','--threads','4','--memory','8192','--out',str(O/g)]
 (O/(g+'.command.json')).write_text(json.dumps(cmd,indent=2))
 with (O/(g+'.console.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError((g,r.returncode))
 print('SOURCE_X_COUNTS_COMPLETE',g,flush=True)
(O/'COMPLETE.json').write_text(json.dumps({'stage':'source_wide_X_counts','complete':True}))
