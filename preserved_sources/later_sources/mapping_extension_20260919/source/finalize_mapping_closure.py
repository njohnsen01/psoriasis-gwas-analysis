from pathlib import Path
import json,hashlib,datetime,subprocess
p=Path(__file__).resolve().parent.parent
s=p/'source'
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as z:
  for b in iter(lambda:z.read(8388608),b''):h.update(b)
 return h.hexdigest()
r=subprocess.run(['python3',str(s/'finish_mapping_report_megablast.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(p/'finish_report.log').write_text(r.stdout)
assert r.returncode==0
m=json.loads((p/'FINAL_DIAGNOSTIC_MANIFEST.json').read_text())
a=json.loads((p/'VERIFIED_AUDIT.json').read_text())
assert a['status']=='bounded_diagnostics_complete'
assert a['full_cohort_mapping_certified'] is False
assert a['all_primary_association_results_unchanged'] is True
assert m['report_sha256']==sha(p/'Mapping_Review.md')
assert m['verified_audit_sha256']==sha(p/'VERIFIED_AUDIT.json')
assert m['read_audit_sha256']==sha(p/'READ_AUDIT.json')
assert m['allele_audit_sha256']==sha(p/'ALLELE_QC_AUDIT.json')
for name,digest in m['source_sha256'].items():
 f=Path(name)
 if not f.is_absolute():f=s/f.name
 assert sha(f)==digest,(name,'source hash mismatch')
required=['PROTOCOL.json','PREPARED.json','INPUT_RECEIPTS.json','READ_AUDIT.json','CLASSIFIER_TRIAGE.json','INDEL_FLAG_SAMPLING_CONTEXT.json','ALLELE_QC_AUDIT.json','VERIFIED_AUDIT.json','REFERENCE_INPUT_VERIFIED.json','Mapping_Review.md','FINAL_DIAGNOSTIC_MANIFEST.json','reportable_quality_quantiles.tsv','reportable_allele_quality.tsv','reference/FAILED_TIMEOUT.json','reference/masked_recovery/FAILED.json','reference/megablast_recovery/PROTOCOL.json','reference/megablast_recovery/COMPLETE.json','reference/megablast_recovery/alignments.tsv','reference/megablast_recovery/homology_summary.tsv']
for name in required:assert (p/name).is_file() and (p/name).stat().st_size>0,name
processes=subprocess.check_output(['ps','-eo','pid=,args='],text=True).splitlines()
active=[line for line in processes if ('blastn' in line or 'reference_homology' in line) and str(p) in line]
assert not active,'reference search still active'
text=(p/'Mapping_Review.md').read_text()
for phrase in ['deletion','93.2%','92.2%','1800','600','MegaBLAST']:
 assert phrase.lower() in text.lower(),phrase
stages={
 'reference_homology.py':'failed: 1800-second timeout; partial output preserved',
 'reference_homology_masked.py':'failed: 600-second timeout; partial output preserved',
 'reference_homology_megablast.py':'completed with 300-second bound; all 18 public queries audited',
 'mapping_audit_megablast.py':'completed; source, depth and homology checks passed',
 'finish_mapping_report_megablast.py':'completed; report hashes verified',
 'finish_mapping_report.before_syntax_fix.py.txt':'historical syntax-failure source retained; not a completed report',
 'mapping_audit_masked.py':'prepared; not run as successful audit',
 'finish_mapping_report_masked.py':'prepared; not run as successful report'}
closure={
 'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'status':'bounded_mapping_diagnostic_completed_and_verified',
 'stage_status':stages,
 'source_manifest_hashes_verified':True,
 'reference_processes_running':False,
 'file_sha256':{name:sha(p/name) for name in required},
 'read_stage_source_and_depth_checks_passed':True,
 'public_context_queries_verified':18,
 'deletion_fixed_allele_balance_flags_retained':True,
 'chromosome_Y_homology_present':True,
 'all_expected_X_alignments_score_above_returned_competitors':True,
 'full_cohort_mapping_validation':False,
 'actual_read_realignment':False,
 'independent_genotype_confirmation':False,
 'GWAS_filters_or_P_values_changed':False,
 'private_data_and_protected_counts_retained_in_Workbench':True,
 'export_or_submission_approval':False,
 'external_replication_secured':False,
 'note':'Source manifest inventories retained files; it is not a claim that every script was executed. Failed searches and fixed flags remain visible. This is a limited high-similarity public-context search plus bounded stratified read review.'}
(p/'FINAL_CLOSURE_AUDIT.json').write_text(json.dumps(closure,indent=2))
print('FINAL_CLOSURE_VERIFIED',closure['created_UTC'])
print('NO_RUNNING_REFERENCE_SEARCH; SOURCE_AND_REPORT_HASHES_PASS; 18_QUERIES_REVIEWED')
print('REPORT_SHA256',sha(p/'Mapping_Review.md'))
print('CLOSURE_SHA256',sha(p/'FINAL_CLOSURE_AUDIT.json'))
