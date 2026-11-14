from pathlib import Path
import json,hashlib,datetime
import numpy as np
import pandas as pd
from mapping_extension import O,sha
S=json.loads((O/'STATUS.json').read_text());assert S['status']=='read_stage_complete' and S['all_selected_retrieved']
assert not json.loads((O/'failures.private.json').read_text())
P=json.loads((O/'PREPARED.json').read_text());assert sha(O/'selection.private.tsv')==P['selection_sha256'];assert sha(O/'PROTOCOL.json')==P['protocol_sha256']
R=json.loads((O/'results.private.json').read_text());sel=pd.read_csv(O/'selection.private.tsv',sep='\t',dtype={'IID':str})
assert len(R)==len(sel) and len({r['IID'] for r in R})==len(R) and {r['IID'] for r in R}==set(sel.IID)
sel=sel.set_index('IID');rows=[]
for r in R:
 assert len(r['targets'])==3 and {t['target'] for t in r['targets']}=={'g0','g1','g2'}
 assert r['sex']==int(sel.loc[r['IID'],'SEX']) and r['case_status']==int(sel.loc[r['IID'],'psoriasis_case'])
 p=O/'reads_private'/r['review_id'];assert json.loads((p/'COMPLETE.private.json').read_text())==r
 for t in r['targets']:
  assert t['depth_crosscheck_passed'] and t['expected_target_fraction']==float(sel.loc[r['IID'],t['target']])
  rows.append(dict(IID=r['IID'],sex=r['sex'],ancestry=r['ancestry'],case_status=r['case_status'],center=r['center'],**t))
d=pd.DataFrame(rows);d.to_csv(O/'all_read_metrics.private.tsv',sep='\t',index=False)
summaries=[]
metrics=['GQ','aligned_base_read_depth','allele_informative_read_depth','HQ_template_depth','HQ_target_fraction','MQ_median','low_MQ_fraction','mean_soft_clip']
for family,keys in [('sex_genotype',['target','sex','expected_target_fraction']),('ancestry_sex_status_genotype',['target','ancestry','sex','case_status','expected_target_fraction']),('center_sex_genotype',['target','center','sex','expected_target_fraction'])]:
 for key,z in d.groupby(keys,dropna=False):
  if len(z)<=20:continue
  record=dict(family=family,**dict(zip(keys,key)));record['group_size']='>20; exact denominator withheld'
  for m in metrics:
   vals=pd.to_numeric(z[m],errors='coerce').dropna()
   if len(vals)>20:
    record[m+'_q25']=float(vals.quantile(.25));record[m+'_median']=float(vals.median());record[m+'_q75']=float(vals.quantile(.75))
  summaries.append(record)
pd.DataFrame(summaries).to_csv(O/'reportable_quality_quantiles.tsv',sep='\t',index=False)
flags=d.loc[~d.adequate_HQ_depth|~d.diagnostic_genotype_concordance|d.GQ.lt(20)].copy();flags.to_csv(O/'flagged_diagnostics.private.tsv',sep='\t',index=False)
# No small participant counts or recoverable complements in the external-facing diagnostic.
summary=dict(status='read_audit_complete',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_selected_source_records_verified=True,all_three_targets_per_selected_participant=True,all_depth_crosschecks_passed=True,read_retrieval_failures=False,diagnostic_flags_present=not flags.empty,any_chrY_mate_or_alternative_tag=bool((d.mate_chrY_reads+d.alternative_chrY_reads).gt(0).any()),all_primary_association_results_unchanged=True,full_cohort_mapping_certified=False,independent_genotype_confirmation=False,limits='Bounded stratified sample; primary-aligner MQ and tags cannot exclude every mapping artifact. Indel classification uses a complete fixed haplotype span. Exact participant and flag counts withheld pending coordinated disclosure review.',source_sha256={n:sha(O/n) for n in ['PROTOCOL.json','PREPARED.json','selection.private.tsv','results.private.json','reportable_quality_quantiles.tsv']})
(O/'READ_AUDIT.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
ref=O/'reference/megablast_recovery';assert (ref/'COMPLETE.json').exists(),'Read audit saved; public reference alignment still pending'
C=json.loads((ref/'COMPLETE.json').read_text());assert C['queries']==18 and C['returncode']==0 and sha(ref/'alignments.tsv')==C['output_sha256'];assert sha(ref/'PROTOCOL.json')==C['protocol_sha256']
meta=json.loads((ref/'queries.json').read_text());cols=['query','subject','identity','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen']
a=pd.read_csv(ref/'alignments.tsv',sep='\t',header=None,names=cols);out=[]
for m in meta:
 z=a.loc[a['query'].eq(m['query'])].copy();z=z.loc[z.qstart.le(m['center'])&z.qend.ge(m['center'])]
 z['query_coverage']=(z.qend-z.qstart+1)/z.qlen
 lo=z[['sstart','send']].min(axis=1);hi=z[['sstart','send']].max(axis=1)
 expected=z.subject.eq('chrX')&lo.le(m['start']+m['center']-1)&hi.ge(m['start']+m['center']-1)
 on=z.loc[expected].sort_values('bitscore',ascending=False);assert not on.empty
 row=dict(query=m['query'],target=m['target'],allele=m['allele'],flank=m['flank'],on_target_bitscore=float(on.iloc[0].bitscore),on_target_coverage=float(on.iloc[0].query_coverage))
 for label,mask in [('chrY',z.subject.eq('chrY')),('other_locus',~expected)]:
  zz=z.loc[mask].sort_values('bitscore',ascending=False)
  if len(zz):
   v=zz.iloc[0];row.update({label+'_subject':str(v.subject),label+'_bitscore':float(v.bitscore),label+'_query_coverage':float(v.query_coverage),label+'_identity':float(v.identity),label+'_score_gap':float(on.iloc[0].bitscore-v.bitscore)})
  else:row[label+'_hit']='none covering the variant anchor under this search'
 out.append(row)
pd.DataFrame(out).to_csv(ref/'homology_summary.tsv',sep='\t',index=False)
summary['public_reference_queries_verified']=len(out);summary['whole_reference_high_similarity_context_search_complete']=True
summary['reference_context_is_not_actual_read_realignment']=True
summary['search_method']='MegaBLAST high-similarity recovery, word size 28, DUST seed masking';summary['original_unmasked_search_complete']=False;summary['seed_masked_blastn_search_complete']=False;summary['status']='bounded_diagnostics_complete';summary['reference_summary_sha256']=sha(ref/'homology_summary.tsv')
(O/'VERIFIED_AUDIT.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print(json.dumps(summary,indent=2))

