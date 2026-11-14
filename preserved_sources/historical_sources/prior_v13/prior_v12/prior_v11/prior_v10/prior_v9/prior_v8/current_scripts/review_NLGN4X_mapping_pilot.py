from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';O=I/'read_mapping_QC';F=O/'CRAM_pilot'
rr=json.loads((F/'results.private.json').read_text());assert all(x['depth_crosscheck_passed'] for x in rr)
r=pd.DataFrame(rr);r['sex_label']=r.sex.map({1:'male',2:'female'})
metrics=['raw_primary_nonduplicate_read_depth','MQ_median','HQ_template_depth','HQ_target_fraction']
summary=[]
for (sex,carrier),d in r.groupby(['sex_label','carrier']):
 row={'sex':sex,'carrier':bool(carrier),'sampled_n':len(d),'scope':'preselected pilot only'}
 for m in metrics:
  for q,label in [(0,'min'),(.25,'q25'),(.5,'median'),(.75,'q75'),(1,'max')]:row[m+'_'+label]=float(d[m].quantile(q))
 summary.append(row)
pd.DataFrame(summary).to_csv(F/'stratum_summary.private.tsv',sep='	',index=False)
# Reconcile observed high-quality template allele fractions to source genotypes without a significance claim.
concordant=[]
for x in rr:
 f=x['HQ_target_fraction'];g=x['expected_NLGN4X_dosage'];ok=(f==0 if g==0 else f==1 if g==2 else .2<=f<=.8)
 concordant.append(dict(pilot_id=x['pilot_id'],expected_dosage=g,HQ_target_fraction=f,diagnostic_concordance=ok,rule='0=all reference;2=all target;heterozygote target fraction 0.2-0.8; post-hoc descriptive pilot check, specified after read review'))
pd.DataFrame(concordant).to_csv(F/'genotype_read_concordance.private.tsv',sep='	',index=False)
cols=['query','chrom','identity','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen','qseq','sseq']
xy=pd.read_csv(O/'broad_reference/XY_alignments.tsv',sep='	',header=None,names=cols)
centers=json.loads((O/'broad_reference/query_centers.json').read_text())
xy['includes_variant']=xy.apply(lambda z:z.qstart<=centers[z['query']]<=z.qend,axis=1)
xy.to_csv(O/'broad_reference/alignments_annotated.tsv',sep='	',index=False)
z=xy.loc[xy.chrom.eq('chrY')&xy.includes_variant&(xy.qend-xy.qstart+1).ge(xy.qlen*.5)].sort_values(['query','bitscore'],ascending=[True,False]).groupby('query',as_index=False).head(1).copy()
annotations=[]
for _,a in z.iterrows():
 qi=int(a.qstart);si=int(a.sstart);step=1 if a.send>=a.sstart else -1;at=None
 for qb,sb in zip(a.qseq,a.sseq):
  if qb!='-' and qi==centers[a['query']]:at=dict(query_center_base=qb,aligned_chrY_base=sb,aligned_chrY_coordinate=si if sb!='-' else None,subject_strand='+' if step==1 else '-');break
  if qb!='-':qi+=1
  if sb!='-':si+=step
 assert at is not None
 annotations.append(at)
z=pd.concat([z.reset_index(drop=True),pd.DataFrame(annotations)],axis=1);z.to_csv(O/'broad_reference/best_long_chrY_alignments.private.tsv',sep='	',index=False)
rg=[]
for x in rr:
 for group in x['read_groups']:rg.append(dict(pilot_id=x['pilot_id'],**group))
pd.DataFrame(rg).to_csv(F/'read_group_metadata.private.tsv',sep='	',index=False)
review={'status':'pilot_reviewed','review_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sampled_participants':len(rr),'all_depth_crosschecks_passed':True,'posthoc_descriptive_genotype_read_concordance':all(x['diagnostic_concordance'] for x in concordant),'read_quality_gate_is_diagnostic_not_primary_filter':True,'full_cohort_mapping_certified':False,'alternative_mapping_tags_are_not_exhaustive':True,'primary_XY_reference_only_not_whole_genome_mappability':True,'read_group_fields_available':sorted(set(k for g in rg for k in g)-{'pilot_id'}),'batch_interpretation':'CRAM headers provide finer library/platform-unit fields in sampled participants; cohort-wide extraction and adjustment were not performed. Center-adjusted Step1/Step2 refit is separate.','variant_source_schema':'Stored MatrixTable and both VDS component schemas were checked; GQ/AD-family fields are present, dedicated DP/MQ are absent. AD sum is not full DP, RGQ is not MQ.','release':'HELD','source_hashes':{str(p.relative_to(I)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [F/'results.private.json',O/'broad_reference/XY_alignments.tsv',O/'genotype_quality_by_genotype.private.tsv',Path(__file__)]}}
(F/'REVIEW.private.json').write_text(json.dumps(review,indent=2))
print('READ_MAPPING_REVIEW_RECORDED')
