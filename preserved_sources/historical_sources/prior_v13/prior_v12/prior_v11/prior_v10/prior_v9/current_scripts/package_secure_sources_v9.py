from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile,importlib.metadata as im,re
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');N=I.parent;O=N/'code_release_candidate_20260912_v9';assert not O.exists()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
old=N/'code_release_candidate_20260912_v8/psoriasis_code_source_collection_v8_HELD.zip'
assert sha(old)=='a45e421cdcad8c344f128bd8471179b7259982192f648fc8eda1ab819d93afd3'
files={}
with zipfile.ZipFile(old) as z:
    assert z.testzip() is None
    manifest=json.loads(z.read('MANIFEST.json'))
    for r in manifest['files']:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
    for name in z.namelist():
        assert not name.startswith('/') and '..' not in Path(name).parts
        if not name.endswith('/'):files['prior_v8/'+name]=z.read(name)
active=['audit_supplement_v4.py','build_supplement_v4.py','finalize_supplement_v4.py','render_supplement_v4_proofs.py','build_S1_secure_staging_v3.py','build_S2_secure_integration_v2.py','render_S1_secure_previews_v3.py','package_secure_sources_v9.py']
for name in active:
    p=I/'scripts'/name;code=p.read_text();compile(code,str(p),'exec')
    assert not re.search(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|AIza[0-9A-Za-z_-]{35}',code)
    files['current_scripts/'+name]=p.read_bytes()
paths=[I/'center_refit_v2/COMPLETE.private.json',I/'center_dependent_staging/COMPLETE.private.json',I/'supplement_review_v4/INDEPENDENT_AUDIT.private.json',I/'supplement_staging_v4/BUILD.private.json',I/'supplement_staging_v4/S2_source_manifest.private.json',I/'supplement_staging_v4/review_proofs/RENDER_MANIFEST.private.json',I/'supplement_staging_v4/review_proofs/VISUAL_QA.private.json']
for p in paths:json.loads(p.read_text())
assert json.loads(paths[-1].read_text())['status']=='REVIEWED_WITH_STATED_LIMITS'
refs=[{'path':str(p),'sha256':sha(p)} for p in paths]
files['verification_checkpoint_references.private.json']=json.dumps(refs,indent=2).encode()
files['ACTIVE_SCRIPTS.json']=json.dumps({'current':active,'center_refit_status':'completed and independently audited; do not rerun','supplement_status':'v4 secure integration and review complete; final scientific scope, full S1/S2 content and native format approval remain open','preserved_history':'prior_v8 contains its original time-stamped instructions and prior_v5 source collection; earlier running-status statements are historical','execution':'Run Python files directly; runpy.run_path(...,run_name=__main__) also supplies __file__. Bare exec in python -c does not.','scope':'Source code and checkpoint references, no new scientific tables, figures, reads, genotypes or identifiers. Entire archive remains HELD.'},indent=2).encode()
files['versions_at_packaging.json']=json.dumps({k:im.version(k) for k in ['numpy','pandas','scipy','openpyxl','Pillow','matplotlib','google-cloud-storage']},indent=2).encode()
files['README.md']='''# Psoriasis secure source collection v9 - HELD

This collection preserves v8 byte-for-byte by archive entry under prior_v8 and adds the executed center-dependent v3 builders, independent v4 audit, v4 supplement integration, metadata finalizer and saved-cell review renderer. The center refit and dependent staging have completed. Do not launch duplicate completed analyses. Earlier running-status prose under prior_v8 is historical.

Use ACTIVE_SCRIPTS.json and verification_checkpoint_references.private.json for the current stage. Earlier cohort construction, laboratory processing, model replay, corrected primary/Y code and study-specific protocols remain in the preserved source collection. A source collection is not a fresh-machine executable replay certification.

Within the authorized Workbench, use the existing study roots recorded in these scripts. Verify the frozen source files and hashes before any rerun. New builders create new output directories and intentionally refuse to overwrite an existing draft. To replay, arrange a separate approved working directory and explicitly adapt paths; do not rerun against the frozen primary directories. Execute scripts as Python files, not bare exec from python -c. The v4 metadata finalizer records and verifies the original invocation's receipt-writing interruption; it is not an additional statistical fit.

Current v4 checks reconcile the center estimates with exact raw REGENIE rows and source hashes, dosage flags, covariate/prediction paths, variant-specific N, primary QC flags, zero failed targeted tests, OR/Wald CI and REGENIE P. The later sex family is 16 post-hoc exploratory contrasts with nominal P and Bonferroni/BH correction, source-specific confidence intervals and dosage units. It is not the original sex-test family or independent replication. Estimability does not prove power; the covariance-zero and provider-pruning limitations remain explicit. The preserved 57-sheet map distinguishes available numerical sources, original evidence and unresolved conditional scope.

The v4 S1 is a genuine new XLSX and preserves all v3 cells/styles; it is not a conversion of the original XLSM. S2 remains a secure Markdown integration draft. Review proofs display saved values and are not native Excel or Word certification. Full final-content coverage, clinical/lab methods, any retained conditional/functional claims, full-cohort mapping/fine-batch scope and manuscript/figure integration remain open as recorded. Earlier fresh-process preflight does not establish fresh-machine full replay.

All new scientific data, code and figures remain under the existing organization export restriction. This archive is HELD. Coordinated disclosure review, permitted transfer/access instructions, final formatting and PI/author approval are required before release. No accession, repository deposit or release permission is invented. Eczema files and jobs are outside this workflow.
'''.encode()
manifest={'status':'HELD source collection','created_UTC':datetime.now(timezone.utc).isoformat(),'files':[{'path':n,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)} for n,b in sorted(files.items())]}
files['MANIFEST.json']=json.dumps(manifest,indent=2).encode();O.mkdir();archive=O/'psoriasis_code_source_collection_v9_HELD.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for n,b in files.items():z.writestr(n,b)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for r in manifest['files']:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
for name in ['README.md','MANIFEST.json','ACTIVE_SCRIPTS.json','verification_checkpoint_references.private.json']:(O/name).write_bytes(files[name])
receipt={'archive':str(archive),'sha256':sha(archive),'files':len(files),'bytes':archive.stat().st_size,'all_internal_hashes_verified':True,'release':'HELD'}
(O/'ARCHIVE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+chr(10));print('V9_SOURCE_ARCHIVE_VERIFIED; remains HELD')
