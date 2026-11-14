from pathlib import Path
import os,json,hashlib
import hail as hl
N=Path(__file__).resolve().parents[1];Q=N.parent/'sample_lead_qc_20260910';O=N/'source_X_sites';O.mkdir(exist_ok=True)
jar=Q/'gcs-connector-hadoop3-latest.jar'
hl.init(local='local[2]',default_reference='GRCh38',log=str(O/'hail.log'),gcs_requester_pays_configuration=os.environ['GOOGLE_PROJECT'],spark_conf={'spark.jars':str(jar),'spark.hadoop.fs.gs.impl':'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem','spark.hadoop.fs.gs.requester.pays.mode':'AUTO','spark.hadoop.fs.gs.requester.pays.project.id':os.environ['GOOGLE_PROJECT']})
uri=os.environ['CDR_STORAGE_PATH'].rstrip('/')+'/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt'
mt=hl.filter_intervals(hl.read_matrix_table(uri),[hl.parse_locus_interval('chrX',reference_genome='GRCh38')])
t=mt.rows().select('filters');t.export(str(O/'chrX_site_filters.tsv.bgz'))
(O/'COMPLETE.json').write_text(json.dumps({'source':uri,'hail_version':hl.__version__,'status':'site_filters_exported','participant_fields_exported':False}))
hl.stop();print('SOURCE_X_SITE_FILTERS_COMPLETE',flush=True)
