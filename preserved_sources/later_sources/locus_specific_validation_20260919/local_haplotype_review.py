from pathlib import Path
import sys,json,collections,hashlib,statistics,re
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=B/'locus_specific_validation_20260919';R=B/'technical_followup_20260919/paired_reads';sys.path.insert(0,str(O/'vendor'));import pysam,edlib
POS=5998724;FA=pysam.FastaFile(str(B/'mapping_extension_20260919/reference/Homo_sapiens_assembly38.fasta'));start=POS-501;ref=FA.fetch('chrX',start,POS+500).upper();assert ref[500]=='A';alt=ref[:500]+'G'+ref[501:];rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
rows=json.loads((O/'read_transitions.private.json').read_text());chosen=[r for r in rows if r['flagged'] and (r['bwa']['status'].startswith('HQ_') or r['minimap']['status'].startswith('HQ_'))];keys={(r['q'],r['end']) for r in chosen};meta=json.loads((R/'query_metadata.private.json').read_text());seqs={}
for name,e in [('R1.private.fastq','1'),('R2.private.fastq','2'),('orphans.private.fastq',None)]:
 p=R/name
 if not p.exists():continue
 with p.open() as f:
  while True:
   h=f.readline().rstrip()
   if not h:break
   s=f.readline().rstrip().upper();f.readline();qual=f.readline().rstrip();q=h[1:].split()[0];end=e if e else next(iter(meta[q]['read_hashes']))
   if (q,end) in keys:assert hashlib.sha256(s.encode()).hexdigest()==meta[q]['read_hashes'][end];seqs[(q,end)]=(s,qual)
assert set(seqs)==keys
out=[]
def gap_near(a):
 p=a['pos']
 for ns,op in re.findall('('+chr(92)+'d+)([MIDNSHP=X])',a['CIGAR']):
  n=int(ns)
  if op in 'ID' and abs(p-POS)<=50:return True
  if op in 'MDN=X':p+=n
 return False
for r in chosen:
 s,qual=seqs[(r['q'],r['end'])];metrics=[]
 for orientation,t in [('+',s),('-',rc(s))]:
  a=edlib.align(t,ref,mode='HW',task='locations');b=edlib.align(t,alt,mode='HW',task='locations')
  metrics.append(dict(orientation=orientation,REF_edits=a['editDistance'],ALT_edits=b['editDistance'],REF_best_overlaps_site=any(lo is not None and lo<=500<=hi for lo,hi in a['locations']),ALT_best_overlaps_site=any(lo is not None and lo<=500<=hi for lo,hi in b['locations'])))
 best=min(metrics,key=lambda x:min(x['REF_edits'],x['ALT_edits']));a,b=best['REF_edits'],best['ALT_edits'];best['preference']='ALT' if b<a else 'REF' if a<b else 'ambiguous'
 out.append(dict(q=r['q'],end=r['end'],review_id=r['review_id'],minimap_status=r['minimap']['status'],bwa_status=r['bwa']['status'],minimap_gap_near=gap_near(r['minimap']),bwa_gap_near=gap_near(r['bwa']),same_target_cycle=bool(r['minimap']['obs'] and r['bwa']['obs'] and r['minimap']['obs']['cycle']==r['bwa']['obs']['cycle']),**best))
(O/'local_alignment.private.json').write_text(json.dumps(out,indent=2));conflict=[r for r in out if {r['minimap_status'],r['bwa_status']}=={'HQ_REF','HQ_ALT'}];assert conflict
agg=dict(any_conflicting_calls=True,all_conflicts_change_target_read_cycle=all(not x['same_target_cycle'] for x in conflict),all_conflicts_have_nearby_gap_in_at_least_one_alignment=all(x['minimap_gap_near'] or x['bwa_gap_near'] for x in conflict),any_conflicts_have_nearby_gap=any(x['minimap_gap_near'] or x['bwa_gap_near'] for x in conflict),edit_preferences_among_conflicts=sorted({x['preference'] for x in conflict}),all_conflicts_edit_best_overlaps_site=all(x['REF_best_overlaps_site'] and x['ALT_best_overlaps_site'] for x in conflict),all_conflicts_HAVE_additional_edits_beyond_target=all(min(x['REF_edits'],x['ALT_edits'])>0 for x in conflict),limitation='Unit-cost edit-distance comparison of two reference haplotypes; does not model other variants or provide genotype likelihoods')
(O/'LOCAL_ALIGNMENT_AGGREGATE.json').write_text(json.dumps(agg,indent=2));print('LOCAL_ALIGNMENT_COMPLETE',flush=True)
