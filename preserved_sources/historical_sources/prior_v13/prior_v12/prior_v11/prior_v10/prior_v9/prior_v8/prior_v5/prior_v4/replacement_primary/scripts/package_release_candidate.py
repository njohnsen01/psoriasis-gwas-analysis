from pathlib import Path
import json,zipfile,hashlib,datetime,importlib.metadata as im,re
N=Path(__file__).resolve().parents[1];R=Path('/home/jupyter/psoriasis_reproducibility_20260910');O=N/'code_release_candidate_20260911'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 O.mkdir(exist_ok=False)
 old=R/'psoriasis_reproducibility_code_v3.zip';assert sha(old)=='c83a97f3dfab572abd1640edea7b42d9160eaef1651f7003b058df3d78f14b0f'
 validation=N/'strict_input_bootstrap_validation_attempt2/VALIDATION.private.json';v=json.loads(validation.read_text());assert v['status']=='validated' and len(v['comparison'])==8
 files={'historical_replay/psoriasis_reproducibility_code_v3.zip':old}
 for p in sorted((N/'scripts').glob('*.py')):
  if '.' in p.stem:continue
  text=p.read_text();compile(text,str(p),'exec')
  assert not re.search(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|ya29\.[A-Za-z0-9_-]{30,}|AIza[A-Za-z0-9_-]{35}',text),'Possible credential material requires review'
  files['replacement_primary/scripts/'+p.name]=p
 for name in ['analysis_amendment.json','autosomal_tasks.json','pilot_sources.json','autosomal_stratified_analysis_plan.json','X_stratified_analysis_plan.json']:
  files['replacement_primary/'+name]=N/name
 files['replacement_primary/recorded_primary_step1_command.json']=N/'step1/command.json'
 files['replacement_primary/input_reconstruction_validation.private.json']=validation
 versions={}
 for p in ['numpy','pandas','scipy','matplotlib','google-cloud-storage','google-cloud-bigquery','zstandard','pgenlib','statsmodels','hail']:
  try:versions[p]=im.version(p)
  except im.PackageNotFoundError:versions[p]='not installed in packaging environment; see execution-specific environment records'
 (O/'versions.json').write_text(json.dumps(versions,indent=2));files['versions.json']=O/'versions.json'
 readme='''Corrected psoriasis analysis code release candidate

This archive remains inside the All of Us Workbench and is HELD for release review. It is not evidence of an approved export route or public code availability. It adds the current corrected primary-analysis scripts to the unchanged, separately identified v3 cohort, raw-laboratory and follow-up replay archive. Historical results are not relabeled as corrected primary results.

The replacement scripts depend on authorized Controlled Tier v8 access, the retained study phenotype/covariates, provider QC files, the complete echo_v4_r2 ancestry file, the documented Step 1 genotype panel and the frozen autosomal checkpoints. These participant/genotype input files are not included. Source and workspace paths are intentionally retained for review; authorized users must supply equivalent approved inputs and adapt their workspace root paths. A clean external-environment replay has not been performed.

Input reconstruction was tested against all eight retained primary phenotype, covariate and keep files. reconstruct_strict_primary_inputs.py takes --output NEW_DIRECTORY and optionally --reference RETAINED_INPUT_DIRECTORY. It restores only the exact pinned coverage-exclusion source if absent, verifies its checksum, applies sample QC and compares every participant key and covariate value. Participant input files it creates must remain in Workbench.

The recorded_primary_step1_command.json records the executed pooled Step 1 command. refit_stratified_step1.py builds and validates group-specific prediction files. The corrected cohort audit supplies group assignment and observation summaries. The X sequence is candidate selection, count extraction, genotype projection, female-control HWE, exact source-site FILTER reconciliation, pooled association, subgroup counts and association. Existing driver scripts enforce their prerequisites. The autosomal runner re-QCs the frozen archived common-variant panel; it is not a maximal source-wide autosomal scan. Pilot and per-chunk receipts are required. Do not launch duplicate runners against an active output directory.

finalize_pooled_calibration.py calculates full-output diagnostics. run_autosomes_stratified.py completes five group analyses and their diagnostics. build_submission_assets.py produces pooled and subgroup Manhattan/QQ plots, complete regional X plots, exact private result tables and disclosure-screened review views. Every plotted association is read from checksum-validated results. QQ ranking precedes deterministic display thinning; lambda uses all finite tests at unit weight. Regional plots do not encode LD or establish causality. The cohort-minus-analysis denominator check suppresses all three analyzed-N fields for a row if any such difference is 1-20; cross-document disclosure review remains necessary.

At packaging, the pooled and X analyses were complete and the autosomal subgroup run was still active. The final results, visual review, table concordance, author approvals and approved export route remain separate requirements. Freeze a final results manifest after the run completes. This archive contains source and provenance, not participant/genotype data files, notebook outputs or a genome-wide summary-statistic release.

Codex assisted with code and figure preparation. The authors must verify the analyses, interpretation and final AI assistance disclosure. No author attestation is made here.
'''
 (O/'README.md').write_text(readme);files['README.md']=O/'README.md'
 records=[{'path':name,'bytes':p.stat().st_size,'sha256':sha(p)} for name,p in sorted(files.items())]
 manifest={'status':'HELD release candidate','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':records,'source_script_count':sum(k.endswith('.py') for k in files),'input_reconstruction':'8 retained input files matched','legacy_archive_sha256':sha(old),'release_review_required':['participant-data exclusion audit of nested historical archive','credential and controlled-path review','approved organization export process','final results and methods concordance']}
 (O/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));files['MANIFEST.json']=O/'MANIFEST.json'
 target=O/'psoriasis_code_release_candidate_v4.zip'
 with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
  for name,p in sorted(files.items()):z.write(p,name)
 with zipfile.ZipFile(target) as z:
  assert z.testzip() is None
  for r in records:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
 (O/'ARCHIVE_RECEIPT.json').write_text(json.dumps({'path':str(target),'bytes':target.stat().st_size,'sha256':sha(target),'status':'created and checksum-verified inside Workbench; not exported'},indent=2))
 print('CODE_RELEASE_CANDIDATE_VALIDATED',len(files),target.stat().st_size,sha(target),flush=True)
if __name__=='__main__':main()
