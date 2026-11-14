
import json, os, datetime, traceback
from pathlib import Path
from google.cloud import storage
w=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/X_source_reanalysis_20260910')
out=w/'private_provider';out.mkdir(exist_ok=True)
c=storage.Client();bk=c.bucket('vwb-aou-datasets-controlled',user_project=os.environ.get('GOOGLE_PROJECT'))
state={'status':'running','completed':[]}
def save():
    state['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(w/'download_status.json').write_text(json.dumps(state,indent=2))
try:
    for item in json.loads((w/'provider_manifest.json').read_text()):
        name=item['name'];dst=out/Path(name).name;state['current']=name;save()
        blob=bk.blob(name,generation=int(item['generation']))
        if not dst.exists():
            tmp=Path(str(dst)+'.partial');blob.download_to_filename(str(tmp),if_generation_match=int(item['generation']),checksum='crc32c',timeout=600);assert tmp.stat().st_size==item['bytes'];tmp.rename(dst)
        assert dst.stat().st_size==item['bytes']
        state['completed'].append({'name':name,'bytes':dst.stat().st_size,'generation':item['generation'],'download_checksum':'google-cloud-storage verified'});save()
    state['status']='complete';save()
except Exception as e:
    state.update(status='failed',error=str(e),traceback=traceback.format_exc());save();raise
