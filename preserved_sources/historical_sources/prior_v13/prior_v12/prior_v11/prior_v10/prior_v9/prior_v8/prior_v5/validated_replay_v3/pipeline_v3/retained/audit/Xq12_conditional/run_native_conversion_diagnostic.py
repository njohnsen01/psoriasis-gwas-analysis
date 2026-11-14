import os
import subprocess, pathlib, json
p = pathlib.Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/Xq12_conditional/private_genotypes/source_qc')
root = p.parent.parent
cfg = json.loads(pathlib.Path(f'{os.environ["PSO_WORKSPACE_MOUNT"].rstrip("/")}/psoriasis_regenie_loco_indexed_range_v1/production_{os.environ["PSO_PRIMARY_RUN_ID"]}/chrX/indexed_range_config.json').read_text())
commands = [
 ['/opt/workbench-tools/samtools/bin/bcftools','norm','--multiallelics','-any','--multi-overlaps','0','-Ov','-o',str(p/'two_native_split.vcf'),str(p/'two_native.vcf')],
 ['/opt/workbench-tools/samtools/bin/bcftools','annotate','--set-id','%CHROM:%POS:%REF:%FIRST_ALT','-Ov','-o',str(p/'two_native_identified.vcf'),str(p/'two_native_split.vcf')],
 ['/opt/workbench-tools/binaries/bin/plink2','--vcf',str(p/'two_native_identified.vcf'),'--double-id','--extract',str(root/'two_ids.txt'),'--update-sex',cfg['local_sex_update'],'--keep',cfg['local_keep'],'--split-par','hg38','--make-pgen','--threads','2','--memory','4000','--out',str(p/'native_model')],
 ['/opt/workbench-tools/binaries/bin/plink2','--pfile',str(p/'native_model'),'--export','A','--threads','2','--memory','4000','--out',str(p/'native_model')]
]
for i, command in enumerate(commands):
 r = subprocess.run(command,capture_output=True,text=True)
 (p/f'native_conversion_{i}.log').write_text(r.stdout+'\n'+r.stderr)
 if r.returncode:
  raise RuntimeError(f'Native source diagnostic step {i} failed; review private log')
(root/'native_conversion_commands.json').write_text(json.dumps(commands,indent=2))
print('Native VCF conversion completed; private genotype files retained in Workbench.')
