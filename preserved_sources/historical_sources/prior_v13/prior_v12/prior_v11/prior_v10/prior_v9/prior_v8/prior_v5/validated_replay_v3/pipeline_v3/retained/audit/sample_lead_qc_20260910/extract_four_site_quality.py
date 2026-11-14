from pathlib import Path
import os,json,traceback
R=Path(__file__).resolve().parent
try:
 import hail as hl
 hl.init(local='local[4]',default_reference='GRCh38',log=str(R/'hail_lead_qc_standalone.log'),gcs_requester_pays_configuration=os.environ['GOOGLE_PROJECT'],spark_conf={'spark.jars':str(R/'gcs-connector-hadoop3-latest.jar'),'spark.hadoop.fs.gs.impl':'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem','spark.hadoop.fs.gs.requester.pays.mode':'AUTO','spark.hadoop.fs.gs.requester.pays.project.id':os.environ['GOOGLE_PROJECT']})
 uri=os.environ['CDR_STORAGE_PATH'].rstrip('/')+'/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt'
 ints=[hl.Interval(hl.Locus('chrX',p,'GRCh38'),hl.Locus('chrX',p+1,'GRCh38'),includes_start=True,includes_end=False) for p in [5998724,22917619,67736589,67842269]]
 m=hl.filter_intervals(hl.read_matrix_table(uri),ints)
 m=m.select_rows('filters').select_cols().select_entries('GT','GQ','RGQ','FT','AD')
 rows=m.rows().collect();print('ROW_KEYS',[(str(x.locus),x.alleles) for x in rows],flush=True);assert len(rows)==4
 (R/'hail_rows.json').write_text(json.dumps([{'locus':str(x.locus),'alleles':x.alleles,'filters':list(x.filters)} for x in rows],indent=2))
 m.entries().export(str(R/'four_sites_entries.private.tsv.bgz'))
 (R/'hail_status.json').write_text(json.dumps({'status':'complete','rows':len(rows),'source':uri,'hail':hl.__version__}))
except Exception:
 (R/'hail_status.json').write_text(json.dumps({'status':'failed','error':traceback.format_exc()}));raise
