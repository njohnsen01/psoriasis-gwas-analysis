import os
import subprocess,json,hashlib
from pathlib import Path
workspace_bucket = os.environ["WORKSPACE_BUCKET"].removeprefix("gs://").strip("/")
u=[f'gs://{workspace_bucket}/psoriasis_regenie_loco_indexed_range_v1/production_{os.environ["PSO_PRIMARY_RUN_ID"]}/chrX/secure_step2_chunks/chrX_pso_chunk_0022.pgen', f'gs://{workspace_bucket}/psoriasis_regenie_loco_indexed_range_v1/production_{os.environ["PSO_PRIMARY_RUN_ID"]}/chrX/secure_step2_chunks/chrX_pso_chunk_0022.pvar.zst', f'gs://{workspace_bucket}/psoriasis_regenie_loco_indexed_range_v1/production_{os.environ["PSO_PRIMARY_RUN_ID"]}/chrX/secure_step2_chunks/chrX_pso_step2.psam']
d=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/Xq12_conditional/private_genotypes')
for x in u:
 p=d/x.rsplit("/",1)[-1]
 subprocess.run(["/usr/bin/gsutil","-q","cp",x,str(p)],check=True)
 print(p.name,p.stat().st_size,flush=True)
(d/"restored.json").write_text(json.dumps({p.name:{"size":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in d.iterdir() if p.suffix in [".pgen",".zst",".psam"]},indent=2))
