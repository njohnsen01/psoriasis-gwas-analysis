from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile,sys
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
I=N/'reconciliation_20260912'; V=I/'supplement_staging_v6'
O=N/'code_release_candidate_20260912_v11'
def digest(data): return hashlib.sha256(data).hexdigest()
def sha(p): return digest(Path(p).read_bytes())
build=json.loads((V/'BUILD.private.json').read_text())
qa=json.loads((V/'review_proofs/VISUAL_QA.private.json').read_text())
assert qa['status']=='selected_proofs_visually_reviewed'
assert qa['render_manifest_sha256']==sha(V/'review_proofs/RENDER_MANIFEST.private.json')
for name,h in build['output_hashes'].items(): assert sha(V/name)==h
prior=N/'code_release_candidate_20260912_v10/psoriasis_code_source_collection_v10_HELD.zip'
assert sha(prior)=='6cf33af23ff554aad6790ad65e65ee17da29b32fea38fb7c9aecd9270cbb6240'
entries={}
with zipfile.ZipFile(prior) as z:
 assert z.testzip() is None
 for name in z.namelist():
  assert not name.startswith('/') and '..' not in Path(name).parts
  if not name.endswith('/'): entries['prior_v10/'+name]=z.read(name)
for name in ['build_supplement_v6.py','build_supplement_v6_revision.py','review_supplement_v6.py','package_secure_sources_v11.py']:
 entries['supplement_v6/scripts/'+name]=(I/'scripts'/name).read_bytes()
for name in ['BUILD.private.json','review_proofs/RENDER_MANIFEST.private.json','review_proofs/VISUAL_QA.private.json']:
 entries['supplement_v6/receipts/'+name.replace('/','__')]=(V/name).read_bytes()
entries['README.md']=b'''# Psoriasis source collection v11 - HELD
Preserves every v10 source entry under prior_v10 and adds the executed v6 source-reconciliation, ancestry leave-one-group-out, figure-inventory and proof-rendering scripts with audit receipts. The active builder is build_supplement_v6_revision.py, which includes QC-excluded historical X targets and caches equivalent style-pair checks while comparing every original cell. The earlier builder is retained for provenance; the failed initial style-proxy check remains separately preserved in Workbench. Scientific tables, images and participant-level data are not included as new v6 payloads. Existing historical entries remain unchanged.

Run stage scripts only inside the authorized Workbench against their pinned dependencies. They deliberately refuse to overwrite existing staging directories. This is a verified source collection, not a fresh-machine or end-to-end replay certification. Native manuscript/S1/S2 integration, coordinated disclosure, remaining science and author approvals remain open. No public release is authorized by this archive.
'''
manifest={'created_UTC':datetime.now(timezone.utc).isoformat(),'release':'HELD','prior_archive_sha256':sha(prior),'Python':sys.version,'files':[{'path':name,'bytes':len(data),'sha256':digest(data)} for name,data in sorted(entries.items())]}
entries['MANIFEST.json']=json.dumps(manifest,indent=2).encode()
assert not O.exists(); O.mkdir()
out=O/'psoriasis_code_source_collection_v11_HELD.zip'
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for name,data in sorted(entries.items()): z.writestr(name,data)
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None and len(z.namelist())==len(entries)
 for record in manifest['files']:
  data=z.read(record['path']); assert len(data)==record['bytes'] and digest(data)==record['sha256']
receipt={'created_UTC':datetime.now(timezone.utc).isoformat(),'archive':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'entries':len(entries),'all_CRC_and_manifest_hashes_verified':True,'status':'HELD; not full replay or public release','packaging_script_sha256':sha(__file__)}
(O/'ARCHIVE_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))

