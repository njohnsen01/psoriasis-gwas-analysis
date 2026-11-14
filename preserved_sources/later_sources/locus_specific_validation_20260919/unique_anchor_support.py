from pathlib import Path
import sys,json,collections,hashlib,statistics
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');T=B/'technical_followup_20260919';R=T/'paired_reads';O=B/'locus_specific_validation_20260919';REF=B/'mapping_extension_20260919/reference/Homo_sapiens_assembly38.fasta'
sys.path.insert(0,str(O/'vendor'));import ahocorasick,pysam
W=[31,51,71,101];POS=5998724;rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
fa=pysam.FastaFile(str(REF));assert fa.fetch('chrX',POS-1,POS)=='A';patterns=collections.defaultdict(list);public=[]
for w in W:
 s=fa.fetch('chrX',POS-1-w//2,POS+w//2).upper();assert len(s)==w and set(s)<=set('ACGT')
 for base in 'ACGT':
  motif=s[:w//2]+base+s[w//2+1:];patterns[motif].append((w,base,'+'));patterns[rc(motif)].append((w,base,'-'))
 public.append(dict(width=w,reference_motif=s,variant_index_zero_based=w//2))
auto=ahocorasick.Automaton()
for s,vals in patterns.items():auto.add_word(s,vals)
auto.make_automaton();hits={w:set() for w in W}
for chrom in fa.references:
 seq=fa.fetch(chrom).upper()
 for end,vals in auto.iter(seq):
  for w,base,strand in vals:hits[w].add((chrom,end-w+1))
print('FULL_REFERENCE_EXACT_SEARCH_COMPLETE',flush=True)
unique={w:hits[w]=={('chrX',POS-1-w//2)} for w in W}
(O/'PUBLIC_ANCHOR_SPECIFICITY.json').write_text(json.dumps([dict(**p,reference_unique=unique[p['width']],occurrences=[dict(chrom=c,start0=start) for c,start in sorted(hits[p['width']])]) for p in public],indent=2))
meta=json.loads((R/'query_metadata.private.json').read_text());expected={(q,e) for q,m in meta.items() for e in m['read_hashes']};seen=set();obs=collections.defaultdict(list)
for name,e in [('R1.private.fastq','1'),('R2.private.fastq','2'),('orphans.private.fastq',None)]:
 p=R/name
 if not p.exists():continue
 with p.open() as f:
  while True:
   header=f.readline().rstrip()
   if not header:break
   seq=f.readline().rstrip().upper();plus=f.readline().rstrip();qual=f.readline().rstrip();q=header[1:].split()[0];end=e if e else next(iter(meta[q]['read_hashes']))
   assert plus.startswith('+') and len(seq)==len(qual);assert (q,end) in expected and (q,end) not in seen;seen.add((q,end));assert hashlib.sha256(seq.encode()).hexdigest()==meta[q]['read_hashes'][end]
   for stop,vals in auto.iter(seq):
    for w,base,strand in vals:
     if not unique[w]:continue
     start=stop-w+1;qs=[ord(c)-33 for c in qual[start:stop+1]];assert len(qs)==w
     if min(qs)<20:continue
     obs[(meta[q]['review_id'],w,q)].append(dict(call='REF' if base=='A' else 'ALT' if base=='G' else 'OTHER',BQ=min(qs),base=base,end=end,strand=strand))
assert seen==expected
prior=json.loads((R/'paired_genotype_review.private.json').read_text());people={x['review_id']:x['expected'] for x in prior if x['model']=='original' and x['target']=='g2'};flagged={x['review_id'] for x in prior if x['model']=='minimap_paired' and x['target']=='g2' and x['adequate'] and not x['balance_pass']}
counts=collections.defaultdict(collections.Counter);conflicts=collections.Counter()
for (rid,w,q),vals in obs.items():
 top=max(x['BQ'] for x in vals);calls={x['call'] for x in vals if x['BQ']==top}
 if len(calls)!=1:conflicts[(rid,w)]+=1;continue
 counts[(rid,w)][next(iter(calls))]+=1
rows=[]
for rid,expected_fraction in people.items():
 for w in W:
  c=counts[(rid,w)];depth=sum(c.values());frac=c['ALT']/depth if depth else None;balance=frac is not None and (frac<=.1 if expected_fraction==0 else frac>=.9 if expected_fraction==1 else .25<=frac<=.75)
  rows.append(dict(review_id=rid,width=w,reference_unique=unique[w],expected_fraction=expected_fraction,counts=dict(c),HQ_depth=depth,ALT_fraction=frac,adequate=depth>=10,balance_pass=balance,conflicts=conflicts[(rid,w)],previously_flagged=rid in flagged))
(O/'unique_anchor_genotypes.private.json').write_text(json.dumps(rows,indent=2));agg=[]
for w in W:
 z=[x for x in rows if x['width']==w];f=[x for x in z if x['previously_flagged']]
 agg.append(dict(width=w,reference_unique=unique[w],any_adequate_depth=any(x['adequate'] for x in z),any_adequate_depth_balance_flag=any(x['adequate'] and not x['balance_pass'] for x in z),any_inadequate_depth=any(not x['adequate'] for x in z),median_HQ_depth=statistics.median(x['HQ_depth'] for x in z),flagged_all_adequate=all(x['adequate'] for x in f),flagged_all_balance_pass=all(x['balance_pass'] for x in f),flagged_any_adequate_balance_flag=any(x['adequate'] and not x['balance_pass'] for x in f)))
(O/'UNIQUE_ANCHOR_AGGREGATE.json').write_text(json.dumps(agg,indent=2));(O/'UNIQUE_ANCHOR_VERIFIED.json').write_text(json.dumps(dict(all_query_sequences_accounted=True,all_input_sequence_hashes_matched=True,all_reference_contigs_scanned=True,mask='all four target bases and both orientations',minimum_quality='every motif base Phred20',one_observation_per_template=True,phenotype_used=False,not_orthogonal_genotyping=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2));print('UNIQUE_ANCHOR_COMPLETE',flush=True)
