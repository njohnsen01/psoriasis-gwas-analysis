from pathlib import Path
import os,json,hashlib
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[-2:]))
import hail as hl
Y=Path(__file__).resolve().parents[1];N=Y.parent;Q=N.parent/'sample_lead_qc_20260910';O=Y/'site_filters';O.mkdir(exist_ok=True)
jar=Q/'gcs-connector-hadoop3-latest.jar'
hl.init(local='local[2]',default_reference='GRCh38',log=str(O/'hail.log'),gcs_requester_pays_configuration=os.environ['GOOGLE_PROJECT'],spark_conf={'spark.master':'local[2]','spark.jars':str(jar),'spark.hadoop.fs.gs.impl':'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem','spark.hadoop.fs.gs.requester.pays.mode':'AUTO','spark.hadoop.fs.gs.requester.pays.project.id':os.environ['GOOGLE_PROJECT']})
uri=os.environ['CDR_STORAGE_PATH'].rstrip('/')+'/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt'
mt=hl.filter_intervals(hl.read_matrix_table(uri),[hl.parse_locus_interval('chrY',reference_genome='GRCh38')])
t=mt.rows();t=t.select('filters',is_y_nonpar=t.locus.in_y_nonpar());t.export(str(O/'chrY_site_filters.tsv.bgz'))
(O/'COMPLETE.json').write_text(json.dumps({'source':uri,'hail_version':hl.__version__,'status':'site_filters_exported','participant_fields_exported':False}))
hl.stop();print('SOURCE_Y_SITE_FILTERS_COMPLETE',flush=True)
