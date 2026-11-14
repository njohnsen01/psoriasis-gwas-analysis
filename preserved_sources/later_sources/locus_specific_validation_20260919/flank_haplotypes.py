from pathlib import Path
import sys,json,collections,hashlib,statistics
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=B/'locus_specific_validation_20260919';R=B/'technical_followup_20260919/paired_reads';sys.path.insert(0,str(O/'vendor'));import ahocorasick,pysam,edlib
P=5998724;fa=pysam.FastaFile(str(B/'mapping_extension_20260919/reference/Homo_sapiens_assembly38.fasta'));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];specs={};patterns=collections.defaultdict(list)
for d in [20,30]:
 left=fa.fetch('chrX',P-1-d-25,P-1-d).upper();right=fa.fetch('chrX',P+d,P+d+25).upper();gap=fa.fetch('chrX',P-1-d,P+d).upper();assert gap[d]=='A';specs[d]=dict(left=left,right=right,gap=gap,alt=gap[:d]+'G'+gap[d+1:])
 for side,s in [('L',left),('R',right)]:
  patterns[s].append((d,side,'+'));patterns[rc(s)].append((d,side,'-'))
auto=ahocorasick.Automaton()
for s,vals in patterns.items():auto.add_word(s,vals)
auto.make_automaton();hits=collections.defaultdict(list)
for chrom in fa.references:
 seq=fa.fetch(chrom).upper()
 for end,vals in auto.iter(seq):
  for d,side,strand in vals:hits[(d,side,strand)].append((chrom,end-24))
unique={};specificity=[]
for d in specs:
 products=set()
 for strand in ['+','-']:
  for c,a in hits[(d,'L',strand)]:
   for c2,b in hits[(d,'R',strand)]:
    if c!=c2:continue
    distance=b-a if strand=='+' else a-b
    if 25<=distance<=275:products.add((c,min(a,b),max(a,b)+25))
 intended=('chrX',P-1-d-25,P+d+25);unique[d]=products=={intended};specificity.append(dict(inner_half_width=d,reference_unique_pair=unique[d],reference_products=[dict(chrom=c,start0=a,end0=b) for c,a,b in sorted(products)]))
(O/'PUBLIC_FLANK_SPECIFICITY.json').write_text(json.dumps(specificity,indent=2));print('FLANK_REFERENCE_SEARCH_COMPLETE',flush=True)
meta=json.loads((R/'query_metadata.private.json').read_text());prior=json.loads((R/'paired_genotype_review.private.json').read_text());people={x['review_id']:x['expected'] for x in prior if x['model']=='original' and x['target']=='g2'};flagged={x['review_id'] for x in prior if x['model']=='minimap_paired' and x['target']=='g2' and x['adequate'] and not x['balance_pass']};obs=collections.defaultdict(list);seen=set()
for name,e in [('R1.private.fastq','1'),('R2.private.fastq','2'),('orphans.private.fastq',None)]:
 p=R/name
 if not p.exists():continue
 with p.open() as f:
  while True:
   h=f.readline().rstrip()
   if not h:break
   seq=f.readline().rstrip().upper();f.readline();qual=f.readline().rstrip();q=h[1:].split()[0];end=e if e else next(iter(meta[q]['read_hashes']));assert (q,end) not in seen;seen.add((q,end));assert hashlib.sha256(seq.encode()).hexdigest()==meta[q]['read_hashes'][end]
   for orient,s,qq in [('+',seq,qual),('-',rc(seq),qual[::-1])]:
    for d,sp in specs.items():
     if not unique[d]:continue
     l=s.find(sp['left']);r=s.find(sp['right'],max(0,l+25))
     if l<0 or r<l+25:continue
     bq=min(ord(c)-33 for c in qq[l:r+25])
     if bq<20:continue
     hap=s[l+25:r];obs[(meta[q]['review_id'],d,q)].append(dict(haplotype=hap,BQ=bq))
expected={(q,e) for q,m in meta.items() for e in m['read_hashes']};assert seen==expected
counts=collections.defaultdict(collections.Counter);conflicts=collections.Counter()
for (rid,d,q),v in obs.items():
 top=max(x['BQ'] for x in v);h={x['haplotype'] for x in v if x['BQ']==top}
 if len(h)!=1:conflicts[(rid,d)]+=1;continue
 counts[(rid,d)][next(iter(h))]+=1
rows=[]
for rid in people:
 for d,sp in specs.items():
  c=counts[(rid,d)];details=[]
  for hap,n in c.items():
   er=edlib.align(hap,sp['gap'],mode='NW')['editDistance'];ea=edlib.align(hap,sp['alt'],mode='NW')['editDistance'];details.append(dict(haplotype=hap,template_count=n,length_difference=len(hap)-len(sp['gap']),REF_edits=er,ALT_edits=ea,preference='ALT' if ea<er else 'REF' if er<ea else 'ambiguous'))
  depth=sum(c.values());rows.append(dict(review_id=rid,inner_half_width=d,HQ_depth=depth,adequate=depth>=10,previously_flagged=rid in flagged,expected_fraction=people[rid],haplotypes=details,conflicts=conflicts[(rid,d)]))
(O/'flank_haplotypes.private.json').write_text(json.dumps(rows,indent=2));out=[]
for d in specs:
 z=[x for x in rows if x['inner_half_width']==d];f=[x for x in z if x['previously_flagged']];hf=[h for x in f for h in x['haplotypes']]
 out.append(dict(inner_half_width=d,reference_unique_pair=unique[d],any_adequate_depth=any(x['adequate'] for x in z),median_HQ_depth=statistics.median(x['HQ_depth'] for x in z),flagged_all_adequate=all(x['adequate'] for x in f),flagged_all_have_some_support=all(x['HQ_depth']>0 for x in f),flagged_any_nonreference_length=any(h['length_difference']!=0 for h in hf),flagged_all_supported_haplotypes_nonreference_length=bool(hf) and all(h['length_difference']!=0 for h in hf),flagged_haplotype_preferences=sorted({h['preference'] for h in hf}),flagged_any_exact_SNP_only_ALT_haplotype=any(h['ALT_edits']==0 for h in hf)))
(O/'FLANK_HAPLOTYPE_AGGREGATE.json').write_text(json.dumps(out,indent=2));(O/'FLANK_HAPLOTYPES_VERIFIED.json').write_text(json.dumps(dict(all_reference_contigs_scanned=True,all_input_sequences_accounted=True,sequence_hashes_verified=True,minimum_BQ_every_base=20,one_observation_per_template=True,phenotype_used=False,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2));print('FLANK_HAPLOTYPES_COMPLETE',flush=True)
