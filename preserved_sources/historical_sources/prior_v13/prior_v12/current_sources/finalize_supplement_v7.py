from pathlib import Path
import datetime,hashlib,json
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');V=I/'supplement_staging_v7';O=V/'review_proofs'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
b=json.loads((V/'BUILD.private.json').read_text());m=json.loads((O/'RENDER_MANIFEST.private.json').read_text());prior=json.loads((O/'RENDER_MANIFEST.before_heading_keep.private.json').read_text());old={Path(x['file']).name:x['sha256'] for x in prior['proofs']};assert len(m['proofs'])==8;unchanged=[];new=[]
for name,h in b['outputs'].items():assert sha(V/name)==h
for x in m['proofs']:
 p=Path(x['file']);assert sha(p)==x['sha256']
 if p.name in old:assert old[p.name]==x['sha256'];unchanged.append(p.name)
 else:new.append(p.name)
assert len(unchanged)==6 and set(new)=={'prose_final_01.png','prose_final_02.png'}
qa={'status':'selected_proofs_visually_reviewed','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'render_manifest_sha256':sha(O/'RENDER_MANIFEST.private.json'),'all_current_proofs':[Path(x['file']).name for x in m['proofs']],'reviewed_unchanged_proof_hashes':unchanged,'newly_reviewed_after_heading_correction':new,'observations':'Both conditional-results pages, the model-diagnostic page and three LD pages reviewed at full-page scale. Availability reasons and all displayed values/labels were legible without overlap or clipping. Methods/results headings were kept with their next paragraph, and both revised prose pages reviewed.','scope':'Selected columns of all newly added table rows and all manuscript insertion text. Not native Excel/Word or all-column layout certification.','workbook_sha256':sha(V/'S1_reconciled_tables_HELD.xlsx'),'source_tables_verified_by_build':True,'finalizer_sha256':sha(__file__)}
(O/'VISUAL_QA.private.json').write_text(json.dumps(qa,indent=2));print('V7_VISUAL_QA_VERIFIED',len(m['proofs']))
