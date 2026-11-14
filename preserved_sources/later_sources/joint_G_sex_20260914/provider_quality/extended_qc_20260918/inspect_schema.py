from pathlib import Path
import os,json,traceback
R=Path(__file__).resolve().parent
try:
 import hail as hl
 hl.init(local='local[2]',default_reference='GRCh38',log=str(R/'hail_lead_qc_standalone.log'),gcs_requester_pays_configuration=os.environ['GOOGLE_PROJECT'],spark_conf={'spark.jars':'/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/sample_lead_qc_20260910/gcs-connector-hadoop3-latest.jar','spark.hadoop.fs.gs.impl':'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem','spark.hadoop.fs.gs.requester.pays.mode':'AUTO','spark.hadoop.fs.gs.requester.pays.project.id':os.environ['GOOGLE_PROJECT']})
 uri=os.environ['CDR_STORAGE_PATH'].rstrip('/')+'/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt'
 m=hl.read_matrix_table(uri)
 fields={'entry':list(m.entry.dtype.fields),'row':list(m.row.dtype.fields),'hail':hl.__version__}
 if 'info' in fields['row']: fields['info']=list(m.row.info.dtype.fields)
 (R/'schema.json').write_text(json.dumps(fields,indent=2));print(json.dumps(fields))
except Exception:
 import traceback; (R/'schema_error.txt').write_text(traceback.format_exc());raise
