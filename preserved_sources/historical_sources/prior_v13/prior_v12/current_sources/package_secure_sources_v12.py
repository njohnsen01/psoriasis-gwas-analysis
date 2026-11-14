from pathlib import Path
import datetime,hashlib,json,sys,zipfile
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';V=I/'supplement_staging_v7';C=I/'conditional_current_v1';O=N/'code_release_candidate_20260912_v12'
def digest(data):return hashlib.sha256(data).hexdigest()
def sha(p):return digest(Path(p).read_bytes())
prior=N/'code_release_candidate_20260912_v11/psoriasis_code_source_collection_v11_HELD.zip';r=json.loads((prior.parent/'ARCHIVE_RECEIPT.json').read_text());assert sha(prior)==r['sha256'];qa=json.loads((V/'review_proofs/VISUAL_QA.private.json').read_text());assert qa['status']=='selected_proofs_visually_reviewed';assert qa['render_manifest_sha256']==sha(V/'review_proofs/RENDER_MANIFEST.private.json');assert json.loads((C/'AUDIT.private.json').read_text())['status']=='source_model_and_summary_verified'
with zipfile.ZipFile(prior) as z:
 assert z.testzip() is None;entries={'prior_v11/'+name:z.read(name) for name in z.namelist()}
source_names=['conditional_current_v1.py','audit_conditional_current_v1.py','audit_conditional_current_v1_initial_warningcheck.py','build_supplement_v7.py','review_supplement_v7.py','finalize_supplement_v7.py','package_secure_sources_v12.py']
for name in source_names:entries['current_sources/'+name]=(I/'scripts'/name).read_bytes()
for name in ['PROTOCOL.private.json','runtime_versions.json','BASELINE_VERIFIED.private.json','GENOTYPE_CONCORDANCE.private.json','AUDIT.private.json']:
 entries['conditional_v1/receipts/'+name]=(C/name).read_bytes()
for name in ['BUILD.private.json','review_proofs/RENDER_MANIFEST.private.json','review_proofs/VISUAL_QA.private.json']:
 entries['supplement_v7/receipts/'+name.replace('/','__')]=(V/name).read_bytes()
entries['README.md']=b'''# Psoriasis source collection v12 - HELD

Preserves every v11 archive entry under prior_v11. Adds the executed current-primary conditional analysis, independent source/model audit, supplement v7 builder, proof reviewer and finalizer, with receipts. The original warning-check source is preserved for provenance; the current auditor accepts only the exact expected extract/exclude intersection notice after verifying the tested IDs and every conditioning n_used.

Run only inside the authorized Workbench against the frozen primary source files and existing verified Step 1 predictions. The conditional runner refuses to overwrite its output directory and gates conditional testing on reproduction of unconditioned primary estimates. A completed supplement build also refuses overwrite; its interrupted text-match attempt was corrected against the actual source wording.

New payloads are source code, runtime/protocol descriptions and integrity receipts. No new participant genotype/phenotype files, scientific tables or images are included. Exact working outputs remain separately in conditional_current_v1 and supplement_staging_v7 inside Workbench. Their disclosure/transport and final native-document integration remain open. This collection does not certify fresh-machine or complete end-to-end replay, authorized public release, or author approval.
'''
manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'release':'HELD','prior_archive_sha256':sha(prior),'Python':sys.version,'files':[{'path':name,'bytes':len(data),'sha256':digest(data)} for name,data in sorted(entries.items())]};entries['MANIFEST.json']=json.dumps(manifest,indent=2).encode();assert not O.exists();O.mkdir();out=O/'psoriasis_code_source_collection_v12_HELD.zip'
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for name,data in sorted(entries.items()):z.writestr(name,data)
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None and len(z.namelist())==len(entries)
 for x in manifest['files']:
  data=z.read(x['path']);assert len(data)==x['bytes'] and digest(data)==x['sha256']
receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'entries':len(entries),'all_CRC_and_manifest_hashes_verified':True,'status':'HELD; not full replay or public release','packaging_script_sha256':sha(__file__)};(O/'ARCHIVE_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
checkpoint={'status':'secure_v7_conditional_source_verified_proofs_reviewed_v12_packaged_HELD','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'supplement':str(V),'conditional_source':str(C),'conditional_audit_sha256':sha(C/'AUDIT.private.json'),'build_sha256':sha(V/'BUILD.private.json'),'visual_QA_sha256':sha(V/'review_proofs/VISUAL_QA.private.json'),'workbook_sha256':sha(V/'S1_reconciled_tables_HELD.xlsx'),'source_archive':receipt,'completed':'13 unconditioned targets reproduced exactly; 22 current relatedness-adjusted conditional estimates verified; 60 ancestry-LD comparisons, 55 available and five unavailable due to reference call-rate gate. Prior 41 worksheets preserved; current 44. All eight current selected-column/text proofs reviewed.','interpretation':'Reporting windows do not establish independent signals; current conditional results remain exploratory. Corrected primary X remains without a genome-wide significant QC-passing association.','release':'HELD; no new scientific export; Drive Word/PDF/ZIP still earlier author-review snapshot','remaining':['permitted transfer resolution and coordinated disclosure','final native manuscript/S1/S2/figure integration and JID formatting','retained-claim functional source checks','full-cohort mapping/fine-batch and clinical validation or explicit justified author-approved limitations','full-pipeline/fresh-machine replay and real approved access identifiers','PI scope/journal fit and author declarations'],'support':'Existing ticket 67164 freshly checked Open with only original request and no reply; no further support request/reply sent.','eczema':'untouched'};(I/'CONTINUATION_CHECKPOINT_20260913_v7.private.json').write_text(json.dumps(checkpoint,indent=2));print('V7_CONTINUATION_CHECKPOINT_SAVED')
