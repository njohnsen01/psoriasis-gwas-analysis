from pathlib import Path
import json,zipfile,hashlib,datetime,io,importlib.metadata as im,re
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';O=N/'code_release_candidate_20260912_v5'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 O.mkdir(exist_ok=False);files={};orig=N/'code_release_candidate_20260911/psoriasis_code_release_candidate_v4.zip'
 with zipfile.ZipFile(orig) as z:
  assert z.testzip() is None
  for n in z.namelist():
   if n.endswith('/'):continue
   assert not n.startswith('/') and '..' not in Path(n).parts
   files['prior_v4/'+n]=z.read(n)
  with zipfile.ZipFile(io.BytesIO(z.read('historical_replay/psoriasis_reproducibility_code_v3.zip'))) as v3:
   assert v3.testzip() is None
   for n in v3.namelist():
    if n.endswith('/'):continue
    assert not n.startswith('/') and '..' not in Path(n).parts
    files['validated_replay_v3/'+n]=v3.read(n)
 for base,label in [(N/'scripts','replacement_primary/scripts'),(N/'Y_source_20260912/scripts','Y_source/scripts'),(I/'scripts','reconciliation/scripts')]:
  for p in sorted(base.iterdir()):
   if p.is_file() and (p.suffix=='.py' or p.name.endswith('.py.txt')):
    source=p.read_text();assert not re.search(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|ya29\.[A-Za-z0-9_-]{30,}|AIza[A-Za-z0-9_-]{35}',source)
    if p.suffix=='.py':compile(source,str(p),'exec')
    files[label+'/'+p.name]=p.read_bytes()
 for folder,paths in [('replacement_primary',['analysis_amendment.json','autosomal_tasks.json','pilot_sources.json','autosomal_stratified_analysis_plan.json','X_stratified_analysis_plan.json','step1/command.json']),('Y_source',['ANALYSIS_PLAN.json','provider_manifest.json','README_FOR_REVIEW.md'])]:
  base=N if folder=='replacement_primary' else N/'Y_source_20260912'
  for name in paths:files[folder+'/'+name]=(base/name).read_bytes()
 checkpoints=[N/'checklist_item1_complete_20260912/COMPLETE.private.json',N/'freeze_review_20260912/SNAPSHOT_INDEX.private.json',N/'Y_source_20260912/freeze_review/SNAPSHOT_INDEX.private.json',I/'COHORT_CHECKPOINT.private.json',I/'estimate_audit/COMPLETE.private.json',I/'followup_input_audit/COMPLETE.private.json',I/'followup_model_audit/COMPLETE.private.json',I/'followup_model_audit/with_QKI_source_QC/QKI_source_validation.private.json',I/'EHR_final_cohort_sensitivity/COMPLETE.private.json',I/'regional_X/COMPLETE.private.json',I/'regional_X/figures_v2/VISUAL_QA.private.json',I/'QKI_primary_audit/COMPLETE.private.json',I/'followup_input_audit/AUTOSOMAL_ALLELE_CHECK.private.json']
 receipts=[]
 for p in checkpoints:
  assert p.exists(),str(p);receipts.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
 versions={'python':__import__('sys').version}
 for name in ['numpy','pandas','scipy','statsmodels','matplotlib','pgenlib','hail','google-cloud-storage','google-cloud-bigquery','threadpoolctl']:
  try:versions[name]=im.version(name)
  except im.PackageNotFoundError:versions[name]='not installed here; see stage-specific execution records'
 files['versions_at_packaging.json']=json.dumps(versions,indent=2).encode()
 files['verification_checkpoint_references.private.json']=json.dumps(receipts,indent=2).encode()
 active=['audit_exact_estimates.py','audit_followup_inputs.py','audit_followup_models_v2.py','EHR_final_cohort.py','regional_X_data_v3.py','plot_full_X_regions_v2.py','audit_QKI_primary.py','audit_autosomal_followup_alleles_v3.py']
 assert all('reconciliation/scripts/'+x in files for x in active)
 files['ACTIVE_RECONCILIATION_SCRIPTS.json']=json.dumps(active,indent=2).encode()
 files['README.md']=b'''Psoriasis executable source collection v5 for secure review

STATUS HELD. This bundle remains inside the authorized All of Us Workbench. It is not approved for export, a public repository, a journal supplement, or an assertion of end-to-end portability. No participant/genotype matrix or result table is included. Exact workspace paths and source identifiers remain in code for provenance and need release review.

Contents and order
1. validated_replay_v3 contains the unchanged validated replay package, including the original cohort SQL/notebook cells, raw-laboratory extraction/processing, X/ancestry/biomarker/SII models, input reconstruction, versions and replay receipts. Read pipeline_v3/README.md. run.py provides check, cohort, labs, lab-processing, x, biomarkers, ancestry and sii stages. Its required --out must be a NEW empty directory. Supply the approved inputs and the documented --audit-root, --input-root, --cohort-dir, --classified-labs, --ancestry-file and --model-input-dir arguments as appropriate. Historical target selection and source versions remain identified.
2. replacement_primary contains the corrected strict-input reconstruction, pooled and group Step 1, source X projection/count/HWE/site FILTER stages, pooled and subgroup X, fixed historical-panel autosomal QC/association, full-output calibration, cohort audit and figure-generation scripts. The pooled Step 1 executed command is retained. Use the stage plans and prerequisite checks. The autosomal analysis is not a maximal source-wide scan. Original provider X PGEN was removed after retained exact projections were made; these verified projections are needed by later regional scripts.
3. Y_source supplies the later corrected non-PAR Y provider-panel analysis. Read ANALYSIS_PLAN and README_FOR_REVIEW. run_Y_pipeline.py is a verification driver by default; --run-missing and --download-missing are explicit options. It is not a source-wide Y scan.
4. reconciliation adds the September 12 exact-row, input, model, EHR, QKI and full X regional-LD/figure scripts. ACTIVE_RECONCILIATION_SCRIPTS.json identifies successful current versions. Older attempts and .executed_console.py.txt fragments are provenance, not independent successful pipelines. The latter require their recorded console imports/paths and are not standalone command-line tools. Use with_QKI_source_QC as the current model table annotation. Failed-QC targets remain diagnostic only. Original and current heterogeneity families must not be interchanged.
5. verification_checkpoint_references.private.json pins the completed primary, Y, cohort and reconciliation checkpoints by exact secure path and hash. Scientific files remain at those paths; this code bundle does not substitute for them. Prior_v4 preserves the earlier package unchanged, whose completion wording is historical.

Replay and safety
Authorized Controlled Tier v8 access, pinned provider data, retained original phenotype/covariates, complete provider ancestry, the specified Step 1 panel, frozen autosomal panel and retained X projections are prerequisites. No credentials are included. Do not run concurrent drivers against a shared output directory or overwrite frozen results. Check existing completion receipts first. Follow stage-specific environment records; versions_at_packaging supplements, rather than replaces, recorded REGENIE 4.1.2 and PLINK2 versions.

Validation and remaining work
Earlier v3 replay and current same-source reconciliation checkpoints are present. Ordinary newly collected Python scripts compile, and archive hashes are checked. A fresh clean-environment replay of this combined v5 collection has NOT been executed. Paths remain study-specific; there is no claimed universal one-command pipeline. Final S1/S2, all promised numerical/figure integration, additional read/mapping QC, author scientific review, disclosure checks, permitted transfer and actual public access instructions remain separate gates. Full regional plots are descriptive and do not constitute fine-mapping. The X regional figures use DejaVu Sans; final JID typography needs certification.

AI assistance
Codex assisted with analysis code, verification and figure preparation. Authors must review the source, results, limitations and final AI-use statement. No author approval or release attestation is made.
'''
 manifest={'status':'HELD source collection; clean external replay not performed','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':[{'path':n,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()} for n,b in sorted(files.items())]}
 files['MANIFEST.json']=json.dumps(manifest,indent=2).encode();archive=O/'psoriasis_code_source_collection_v5_HELD.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
  for n,b in files.items():z.writestr(n,b)
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  for row in manifest['files']:assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256']
 for n in ['README.md','MANIFEST.json','versions_at_packaging.json','verification_checkpoint_references.private.json','ACTIVE_RECONCILIATION_SCRIPTS.json']:(O/n).write_bytes(files[n])
 (O/'ARCHIVE_RECEIPT.json').write_text(json.dumps({'archive':str(archive),'sha256':sha(archive),'bytes':archive.stat().st_size,'files':len(files),'all_internal_hashes_verified':True,'release':'HELD'},indent=2));print('V5_SECURE_SOURCE_COLLECTION_COMPLETE',sha(archive))
if __name__=='__main__':main()
