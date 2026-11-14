"""Validate a disposable diploid-coded X model input. Original X files unchanged."""
from pathlib import Path
import json,subprocess,time,hashlib,datetime,os
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
R=Path('/home/jupyter/x_followup_20261013');Q=R/'chen_scan_benchmark';O=R/'continuation_20261013';O.mkdir(exist_ok=True)
D=O/'proxy_benchmark';D.mkdir(exist_ok=True)
source=B/'X_projected_source/chunk_0009/input'
for suffix in ['.pgen','.psam']:
 dst=D/('proxy'+suffix)
 if not dst.exists():dst.symlink_to(str(source)+suffix)
with Path(str(source)+'.pvar').open() as f,(D/'proxy.pvar').open('w') as w:
 for l in f:
  if l.startswith('##'):continue
  if l.startswith('#'):w.write(l);continue
  a=l.rstrip().split('\t'); assert a[0] in ['X','chrX','23'];a[0]='1';w.write('\t'.join(a)+'\n')
protocol={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'100-variant implementation benchmark only','model':'Additive + female dominance + additive-by-sex; expanded nuisance model','representation':'Disposable PVAR chromosome 1 wrapper; original IDs, positions, alleles, pgen and psam unchanged. Explicit sex covariate. No diploid biological claim.','test':'Joint 3df Wald; compare independently with R glm Wald before scale-up','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(D/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2)+'\n')
cmd=['/opt/workbench-tools/binaries/bin/plink2','--pfile',str(D/'proxy'),'--keep',str(Q/'keep.private.tsv'),'--covar',str(Q/'covariates.private.tsv'),'--pheno',str(Q/'phenotype.private.tsv'),'--extract',str(Q/'extract100.txt'),'--glm','genotypic','interaction','no-x-sex','no-firth','--parameters','1-28','--tests','1,2,28','--threads','8','--memory','16000','--out',str(D/'fit')]
(D/'COMMAND.json').write_text(json.dumps(cmd,indent=2)); t=time.time()
with (D/'console.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
receipt={'returncode':p.returncode,'elapsed_seconds':time.time()-t,'status':'needs_validation' if p.returncode==0 else 'failed'}
(D/'RUN.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
