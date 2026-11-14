from pathlib import Path
import json,sys,subprocess,collections,statistics,hashlib
O=Path('../../deeper_analysis_20260919').resolve();M=O.parent/'mapping_extension_20260919';R=O/'read_remapping'
assert (R/'VERIFIED.json').exists()
sys.path.insert(0,str(M/'source'));from read_qc_core import TARGETS,observe,self_test
assert self_test()
ref=M/'reference/Homo_sapiens_assembly38.fasta'
s=subprocess.check_output(['/opt/workbench-tools/samtools/bin/samtools','faidx',str(ref),'chrX:5996489-5996506'],text=True);seq=''.join(s.splitlines()[1:]).upper();assert seq[5:8]=='GCT';window=(5996489,5996506,seq,seq[:5]+'G'+seq[8:])
meta=json.loads((R/'query_metadata.private.json').read_text());old=json.loads((M/'results.private.json').read_text())
primary={};yhits=collections.defaultdict(list)
for line in (R/'remapped.private.sam').open():
 if line.startswith('@'):continue
 a=line.rstrip().split(chr(9));flag=int(a[1]);tags={z.split(':',2)[0]:z.split(':',2)[2] for z in a[11:]};score=int(tags['AS']) if 'AS' in tags else None
 if not flag&2304:assert a[0] not in primary;primary[a[0]]=(a,score)
 elif a[2]=='chrY' and score is not None:yhits[a[0]].append(score)
assert set(primary)==set(meta)
by=collections.defaultdict(list)
for q,m in meta.items():
 a,score=primary[q];expected=a[2]=='chrX' and 5995000<=int(a[3])<=6000000
 for ann in m['annotations']:
  j=int(ann['target'][1:]);x=observe(a,TARGETS[j],window if j==1 else None) if expected else None
  hq=x is not None and 30<=int(a[4])<255 and x['BQ']>=20
  by[(m['review_id'],ann['target'])].append(dict(oldcall=ann['old_call'],newcall=x['call'] if hq else None,expected=expected,hq=hq,ygap=score-max(yhits[q]) if yhits[q] and score is not None else None))
def concordant(f,e):return f is not None and (f<=.1 if e==0 else f>=.9 if e==1 else .25<=f<=.75)
review=[]
for person in old:
 for t in person['targets']:
  rows=by.get((person['review_id'],t['target']),[]);assert len(rows)==t['HQ_template_depth']
  oc=collections.Counter(x['oldcall'] for x in rows);nc=collections.Counter(x['newcall'] for x in rows if x['newcall'] is not None)
  assert dict(oc)==t['high_quality_allele_counts'];e=t['expected_target_fraction'];nd=sum(nc.values());od=sum(oc.values())
  of=oc['ALT']/od if od else None;nf=nc['ALT']/nd if nd else None;ox=oc['ALT']/(oc['REF']+oc['ALT']) if oc['REF']+oc['ALT'] else None;nx=nc['ALT']/(nc['REF']+nc['ALT']) if nc['REF']+nc['ALT'] else None
  assert concordant(of,e)==t['diagnostic_genotype_concordance']
  gaps=[x['ygap'] for x in rows if x['ygap'] is not None]
  review.append(dict(review_id=person['review_id'],target=t['target'],sex=person['sex'],old_depth=od,new_depth=nd,old_fraction=of,new_fraction=nf,old_exact_fraction=ox,new_exact_fraction=nx,old_balance_flag=not concordant(of,e),new_balance_flag=not concordant(nf,e),old_exact_flag=not concordant(ox,e),new_exact_flag=not concordant(nx,e),old_depth_flag=od<10,new_depth_flag=nd<10,retention=nd/od if od else None,expected_fraction=sum(x['expected'] for x in rows)/od if od else None,median_X_minus_Y_score=statistics.median(gaps) if gaps else None,has_Y_tie_or_better=any(g<=0 for g in gaps)))
(R/'genotype_review.private.json').write_text(json.dumps(review))
summary=[]
for t in ('g0','g1','g2'):
 v=[x for x in review if x['target']==t];assert len(v)>20
 row=dict(target=t,any_old_balance_flag=any(x['old_balance_flag'] for x in v),any_new_balance_flag=any(x['new_balance_flag'] for x in v),any_old_exact_flag=any(x['old_exact_flag'] for x in v),any_new_exact_flag=any(x['new_exact_flag'] for x in v),any_new_depth_flag=any(x['new_depth_flag'] for x in v),any_Y_tie_or_better=any(x['has_Y_tie_or_better'] for x in v))
 for k in ('retention','expected_fraction','median_X_minus_Y_score'):
  vals=[x[k] for x in v if x[k] is not None];row['median_'+k]=round(statistics.median(vals),4) if len(vals)>20 else None
 summary.append(row)
(R/'GENOTYPE_REVIEW_AGGREGATE.json').write_text(json.dumps(summary,indent=2))
(R/'GENOTYPE_REVIEW_VERIFIED.json').write_text(json.dumps(dict(original_depth_and_alleles_reproduced=True,original_flag_rules_reproduced=True,new_primary_HQ_calls_recomputed=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
print('GENOTYPE_REVIEW_VERIFIED')
for row in summary:print(row)
