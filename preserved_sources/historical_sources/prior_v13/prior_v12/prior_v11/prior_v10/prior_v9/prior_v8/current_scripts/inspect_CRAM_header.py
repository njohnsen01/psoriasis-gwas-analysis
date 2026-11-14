from pathlib import Path
import os,json,subprocess,traceback
import google.auth
from google.auth.transport.requests import Request
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');O=N/'reconciliation_20260912/read_mapping_QC'
try:
 rec=json.loads((O/'first_CRAM_source.private.json').read_text());creds,project=google.auth.default(scopes=['https://www.googleapis.com/auth/devstorage.read_only']);creds.refresh(Request())
 env=dict(os.environ,GCS_OAUTH_TOKEN=creds.token,GCS_REQUESTER_PAYS_PROJECT=os.environ['GOOGLE_PROJECT'])
 with (O/'pilot_CRAM_header.private.sam').open('w') as out,(O/'pilot_CRAM_header.private.log').open('w') as err:
  p=subprocess.run(['samtools','view','-H',rec['uri']],stdout=out,stderr=err,env=env,timeout=60)
 result={'status':'header_read' if p.returncode==0 else 'failed','returncode':p.returncode,'read_records_retrieved':False}
 if p.returncode==0:
  hdr=(O/'pilot_CRAM_header.private.sam').read_text().splitlines();result['read_group_field_names']=sorted({f.split(':',1)[0] for line in hdr if line.startswith('@RG') for f in line.split(chr(9))[1:]});result['reference_fields']=[line for line in hdr if line.startswith('@SQ') and chr(9)+'SN:chrX'+chr(9) in line]
 (O/'pilot_CRAM_header_status.private.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='reference_fields'}),flush=True)
except Exception as e:
 (O/'pilot_CRAM_header_status.private.json').write_text(json.dumps({'status':'failed','error':repr(e),'traceback':traceback.format_exc()},indent=2));raise
