from pathlib import Path
import collections,json,sys,hashlib
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');T=B/'technical_followup_20260919';R=T/'paired_reads';O=B/'locus_specific_validation_20260919';M=B/'mapping_extension_20260919'
sys.path.insert(0,str(M/'source'));import read_qc_core as core
meta=json.loads((R/'query_metadata.private.json').read_text());rev=json.loads((R/'paired_genotype_review.private.json').read_text());flagged={x['review_id'] for x in rev if x['model']=='minimap_paired' and x['target']=='g2' and x['adequate'] and not x['balance_pass']};assert flagged
expected={(q,e) for q,m in meta.items() for e in m['read_hashes']}
def collect(prefix):
 out={};secondaries=collections.defaultdict(list)
 for name in ['paired_'+prefix+'.private.sam','orphan_'+prefix+'.private.sam']:
  if not (R/name).exists():continue
  for line in (R/name).open():
   if line.startswith('@'):continue
   a=line.rstrip().split(chr(9));flag=int(a[1]);q=a[0];end=('1' if flag&64 else '2') if meta[q]['paired'] else next(iter(meta[q]['read_hashes']));key=(q,end)
   tags={x.split(':',2)[0]:x.split(':',2)[2] for x in a[11:] if len(x.split(':',2))==3}
   if flag&2304:
    if a[2]=='chrY':secondaries[key].append(dict(chr=a[2],pos=int(a[3]),AS=int(tags['AS']) if 'AS' in tags else None,flag=flag))
    continue
   assert key not in out
   obs=core.observe(a,core.TARGETS[2]);status='unmapped' if flag&4 else 'off_X' if a[2]!='chrX' else 'no_site_observation' if obs is None else 'low_MQ' if not 30<=obs['MQ']<255 else 'low_BQ' if obs['BQ']<20 else 'HQ_'+obs['call']
   out[key]=dict(q=q,end=end,review_id=meta[q]['review_id'],status=status,chr=a[2],pos=int(a[3]),MQ=int(a[4]),CIGAR=a[5],flag=flag,AS=int(tags['AS']) if 'AS' in tags else None,obs=obs)
 assert set(out)==expected
 for k,v in out.items():v['Y_alternatives']=secondaries[k]
 return out
mm=collect('minimap');bw=collect('bwa');rows=[]
for key in expected:
 a,b=mm[key],bw[key]
 if a['obs'] is not None or b['obs'] is not None:
  rows.append(dict(q=key[0],end=key[1],review_id=a['review_id'],flagged=a['review_id'] in flagged,minimap=a,bwa=b))
(O/'read_transitions.private.json').write_text(json.dumps(rows))
checks={}
for label,predicate in [('all_selected',lambda r:True),('flagged_observations',lambda r:r['flagged'])]:
 z=[r for r in rows if predicate(r)];lost=[r for r in z if r['bwa']['status']=='HQ_ALT' and r['minimap']['status']!='HQ_ALT'];lostref=[r for r in z if r['bwa']['status']=='HQ_REF' and r['minimap']['status']!='HQ_REF']
 checks[label]=dict(any_BWA_ALT_lost_in_minimap=bool(lost),ALT_loss_statuses=sorted({r['minimap']['status'] for r in lost}),any_BWA_REF_lost_in_minimap=bool(lostref),REF_loss_statuses=sorted({r['minimap']['status'] for r in lostref}),any_high_quality_REF_ALT_disagreement=any({r['bwa']['status'],r['minimap']['status']}=={'HQ_REF','HQ_ALT'} for r in z),all_lost_ALT_still_chrX=bool(lost) and all(r['minimap']['chr']=='chrX' for r in lost),any_lost_ALT_has_emitted_Y_alternative=any(r['minimap']['Y_alternatives'] for r in lost),all_lost_ALT_have_emitted_Y_alternative=bool(lost) and all(r['minimap']['Y_alternatives'] for r in lost),all_lost_ALT_have_same_CIGAR=bool(lost) and all(r['minimap']['CIGAR']==r['bwa']['CIGAR'] and r['minimap']['pos']==r['bwa']['pos'] for r in lost),all_lost_ALT_have_same_allele_observation=bool(lost) and all(r['minimap']['obs'] and r['minimap']['obs']['call']=='ALT' for r in lost))
(O/'READ_TRANSITION_AGGREGATE.json').write_text(json.dumps(checks,indent=2));(O/'READ_TRANSITIONS_VERIFIED.json').write_text(json.dumps(dict(query_accounting=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),phenotype_used=False),indent=2));print('READ_TRANSITIONS_COMPLETE',flush=True)
