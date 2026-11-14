from pathlib import Path
import csv,hashlib,json,datetime,zipfile,sys,importlib.metadata as M
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');N=I.parent;V=I/'supplement_staging_v5';F=V/'review_proofs_final'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();now=datetime.datetime.now(datetime.timezone.utc).isoformat()
audit=json.loads((V/'INDEPENDENT_METHODS_AUDIT.private.json').read_text());assert audit['typed_unit_table_verified'] and audit['all_original_v4_cells_preserved']
proofs=sorted(F.glob('*.png'));assert len(proofs)==18
qa={'verified_UTC':now,'reviewed_in_authorized_browser':True,'scope':'All 18 saved-cell/text proof pages visually checked for legibility, overlap, clipping and source-grounded wording. Six unchanged laboratory/phenotype pages reused only after byte-identical hashes were verified. These are selected-column review proofs, not native Excel/Word certification or every-cell visual review.','images':[{'name':p.name,'sha256':sha(p),'review_basis':'earlier visual inspection plus byte-identical final image' if p.name in audit['reviewed_initial_images_identical_to_final'] else 'direct final-image visual inspection'} for p in proofs],'corrections':'Strict-unit standard-name precedence was corrected to match actual source; factor header made explicit and saved factors numeric. Original association results unchanged.','initial_proof_directory':'review_proofs is superseded by review_proofs_final.'}
(F/'VISUAL_QA.private.json').write_text(json.dumps(qa,indent=2))
rows=list(csv.DictReader((V/'methods_sources/Methods_provenance.private.tsv').open(),delimiter='\t'))
for r in rows:
 paths=[Path(v) for v in r.values() if isinstance(v,str) and v.startswith('/home/jupyter/')];hashes=[v for k,v in r.items() if 'sha256' in k.lower()];assert len(paths)==len(hashes)==1;assert paths[0].is_file() and sha(paths[0])==hashes[0]
r=json.loads((V/'BUILD.private.json').read_text());r['final_verification_UTC']=now;r['final_status']='VERIFIED SECURE V5 METHODS INTEGRATION; HELD, NOT FINAL SUBMISSION';r['independent_audit']=str(V/'INDEPENDENT_METHODS_AUDIT.private.json');r['visual_QA']=str(F/'VISUAL_QA.private.json');r['all_provenance_source_hashes_rechecked']=len(rows);r['files']={str(p.relative_to(V)):sha(p) for p in V.rglob('*') if p.is_file() and p.name!='BUILD.private.json' and 'review_proofs/' not in str(p)};(V/'BUILD.private.json').write_text(json.dumps(r,indent=2))
A=N/'code_release_candidate_20260912_v10';assert not A.exists();A.mkdir()
prior=N/'code_release_candidate_20260912_v9/psoriasis_code_source_collection_v9_HELD.zip';assert sha(prior)=='23d04891c5052de67e3d9fce149a7d16fde8af3fed6d0466ada25a6e8f4e9048'
files={}
with zipfile.ZipFile(prior) as z:
 assert z.testzip() is None
 for item in z.infolist():
  if not item.is_dir():files['prior_v9/'+item.filename]=z.read(item.filename)
script_names=['build_methods_integration_v5.py','query_lab_vocabulary_v5.py','audit_lab_concept_usage_v5.py','finalize_methods_integration_v5.py','render_methods_v5_proofs.py','correct_methods_unit_wording_v5.py','audit_methods_v5.py','finalize_secure_methods_package_v5.py']
for name in script_names:
 p=I/'scripts'/name;assert p.is_file();compile(p.read_text(),str(p),'exec');files['methods_v5/scripts/'+name]=p.read_bytes()
for p in [V/'BUILD.private.json',V/'INDEPENDENT_METHODS_AUDIT.private.json',F/'VISUAL_QA.private.json',F/'RENDER_MANIFEST.private.json',V/'methods_sources/lab_vocabulary_query.sql',V/'methods_sources/lab_vocabulary_query_receipt.private.json']:
 files['methods_v5/checkpoints/'+p.name]=p.read_bytes()
refs={str(p):sha(p) for p in [V/'S1_reconciled_tables_HELD.xlsx',V/'S2_reconciliation_integration.private.md',V/'Manuscript_cohort_laboratory_insertions.private.md',V/'BUILD.private.json',V/'INDEPENDENT_METHODS_AUDIT.private.json',F/'VISUAL_QA.private.json',I/'CONTINUATION_CHECKPOINT_v3.private.json']}
files['verification_checkpoint_references.private.json']=json.dumps(refs,indent=2).encode()
versions={'verified_UTC':now,'python':sys.version,'packages':{n:M.version(n) for n in ['numpy','pandas','openpyxl','Pillow']},'scope':'Integration runtime only; prior_v9 retains model/tool versions and prior execution records.'};files['runtime_versions.json']=json.dumps(versions,indent=2).encode()
readme='''# Psoriasis secure source collection v10 — HELD

This source collection preserves every v9 archive entry under prior_v9 and adds executed cohort/laboratory Methods integration, vocabulary/usage checks, rendering and independent verification scripts. It is not an approved public release or a certified portable end-to-end pipeline. No new association fit was performed for v5 Methods integration.

The original cohort/laboratory code and prior frozen analysis/replay records remain in prior_v9. Read its README and active-script list before running any analysis. methods_v5 scripts are execution/provenance records with absolute secure paths and prerequisite files; do not blindly replay them. The initial builder correctly stopped at an unmatched classifier-ID discrepancy, then at empty-string-to-XLSX-blank serialization. The discrepancy was investigated and both issues resolved with saved-source checks. The finalizer and unit-wording correction are one-time migrations; the final builder is guarded against overwriting existing S1. The independent unit audit selects the actual nested AST assignments and verifies the saved unit pairs. Early audit helper errors occurred before a passing receipt and changed no scientific results.

Current authoritative secure integration is supplement_staging_v5: S1 XLSX, S2 Markdown and main-text insertion candidate, with BUILD, independent audit and final visual-QA receipts. Scientific files stay in Workbench. verification_checkpoint_references.private.json gives exact paths/hashes; large workbooks, participant-level data and scientific proof images are not embedded in this source collection.

Release/disclosure review, final scientific coverage and native manuscript/figure integration, full portability/replay and PI/author decisions remain open. Source-code presence and compilation are not evidence that every stage was freshly re-executed. Exact-timestamp SII is a concurrency proxy, not specimen-confirmed same-draw validation. Original source results and the corrected negative primary-X conclusion remain unchanged. Eczema work was untouched.
'''
files['README.md']=readme.encode();active={'current_verification':'methods_v5/scripts/audit_methods_v5.py','one_time_integration_records':script_names,'prior_analysis_instructions':'prior_v9/README.md','release':'HELD','native_format_certification':False,'fresh_machine_full_replay':False};files['ACTIVE_SCRIPTS.json']=json.dumps(active,indent=2).encode()
manifest={'created_UTC':now,'prior_archive_sha256':sha(prior),'files':{k:{'sha256':hashlib.sha256(v).hexdigest(),'bytes':len(v)} for k,v in files.items()},'release':'HELD'};files['MANIFEST.json']=json.dumps(manifest,indent=2).encode()
zpath=A/'psoriasis_code_source_collection_v10_HELD.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for name,data in files.items():z.writestr(name,data)
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None and set(z.namelist())==set(files)
 assert all(hashlib.sha256(z.read(k)).hexdigest()==v['sha256'] for k,v in manifest['files'].items())
for name in ['README.md','MANIFEST.json','ACTIVE_SCRIPTS.json','runtime_versions.json','verification_checkpoint_references.private.json']:(A/name).write_bytes(files[name])
receipt={'verified_UTC':now,'archive':str(zpath),'sha256':sha(zpath),'bytes':zpath.stat().st_size,'entries':len(files),'CRC_verified':True,'all_manifest_hashes_verified':True,'prior_archive_unchanged':True,'release':'HELD'};(A/'ARCHIVE_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
checkpoint={'created_UTC':now,'study':'PSORIASIS ONLY','latest_secure_supplement':str(V),'code_collection':receipt,'artifacts':refs,'completed':'Source-grounded cohort/laboratory Methods integrated in secure v5; all original cells preserved; new table values and source hashes checked; unit rule evaluated from actual source; 18 proof pages reviewed.','remaining':['Approved export and coordinated disclosure review','Final original-57-sheet scientific coverage and retained-locus scope','Required full-cohort mapping/fine-batch and clinical validation or justified author-approved limitations','Native manuscript/figure/supplement integration and final-format QA','Fresh-machine/full replay and approved access/repository instructions','PI/author scope and declarations'],'active_job_for_this_step':None,'Drive_scientific_files_updated':False,'eczema_untouched':True};(I/'CONTINUATION_CHECKPOINT_v4.private.json').write_text(json.dumps(checkpoint,indent=2))
print('V5_FINAL_PACKAGE_VERIFIED',len(proofs),'proofs;',len(rows),'source hashes;',len(files),'archive entries');print('CHECKPOINT',sha(I/'CONTINUATION_CHECKPOINT_v4.private.json'));print('ARCHIVE',receipt['sha256'],receipt['bytes'])
