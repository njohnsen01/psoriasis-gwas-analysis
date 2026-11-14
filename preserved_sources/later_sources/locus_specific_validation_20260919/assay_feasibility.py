from pathlib import Path
import sys,json,subprocess,hashlib,collections
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=B/'locus_specific_validation_20260919';M=B/'mapping_extension_20260919';sys.path.insert(0,str(O/'vendor'));import primer3,pysam
P=5998724;lo=P-1001;fa=pysam.FastaFile(str(M/'reference/Homo_sapiens_assembly38.fasta'));seq=fa.fetch('chrX',lo,P+1000).upper();assert seq[1000]=='A'
args={'PRIMER_OPT_SIZE':22,'PRIMER_MIN_SIZE':20,'PRIMER_MAX_SIZE':26,'PRIMER_OPT_TM':60.0,'PRIMER_MIN_TM':58.0,'PRIMER_MAX_TM':63.0,'PRIMER_PAIR_MAX_DIFF_TM':2.0,'PRIMER_MAX_POLY_X':4,'PRIMER_NUM_RETURN':10,'PRIMER_PRODUCT_SIZE_RANGE':[[400,800]],'PRIMER_EXPLAIN_FLAG':1,'PRIMER_MAX_NS_ACCEPTED':0}
result=primer3.bindings.design_primers({'SEQUENCE_ID':'GRCh38_chrX_rs183018930_reference','SEQUENCE_TEMPLATE':seq,'SEQUENCE_TARGET':[900,201]},args);(O/'PUBLIC_PRIMER3_RESULT.json').write_text(json.dumps(result,indent=2));(O/'PUBLIC_ASSAY_SETTINGS.json').write_text(json.dumps(dict(template_chr='chrX',template_start0=lo,template_end0=P+1000,target_position1=P,settings=args,primer3_version=primer3.__version__,scope='Feasibility only; not an order-ready or validated assay; no common-variant binding-site screen yet'),indent=2))
queries=[];pairs=[]
for i in range(result['PRIMER_PAIR_NUM_RETURNED']):
 l,ll=result['PRIMER_LEFT_'+str(i)];r,rl=result['PRIMER_RIGHT_'+str(i)];left=result['PRIMER_LEFT_'+str(i)+'_SEQUENCE'];right=result['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'];pairs.append(dict(pair=i,left_5to3=left,right_5to3=right,left_start1=lo+l+1,left_end1=lo+l+ll,right_start1=lo+r-rl+2,right_end1=lo+r+1,product_size=result['PRIMER_PAIR_'+str(i)+'_PRODUCT_SIZE'],left_Tm=result['PRIMER_LEFT_'+str(i)+'_TM'],right_Tm=result['PRIMER_RIGHT_'+str(i)+'_TM']))
 queries.extend(['>p'+str(i)+'L',left,'>p'+str(i)+'R',right])
assert pairs;(O/'public_primers.fa').write_text(chr(10).join(queries)+chr(10));(O/'PUBLIC_PRIMER_PAIRS.json').write_text(json.dumps(pairs,indent=2))
cmd=['blastn','-task','blastn-short','-query',str(O/'public_primers.fa'),'-db',str(M/'reference/genome_db'),'-word_size','7','-dust','no','-soft_masking','false','-evalue','1000','-max_target_seqs','100000','-max_hsps','100000','-num_threads','2','-outfmt','6 qseqid sseqid qstart qend sstart send mismatch gaps qlen length bitscore','-out',str(O/'public_primer_blast.tsv')]
(O/'PUBLIC_PRIMER_BLAST_COMMAND.json').write_text(json.dumps(cmd,indent=2));print('PRIMER_CANDIDATES_DESIGNED_START_SPECIFICITY_SEARCH',flush=True)
with (O/'primer_blast.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=600)
qlens={a[1:]:len(b) for a,b in zip(queries[::2],queries[1::2])};hits=collections.defaultdict(list)
for line in (O/'public_primer_blast.tsv').open():
 a=line.split();q,chrom=a[:2];qs,qe,ss,se,mm,gaps,n,length=map(int,a[2:10]);assert n==qlens[q]
 if qs==1 and qe==n and gaps==0 and mm<=3:hits[q].append(dict(chr=chrom,start=min(ss,se),end=max(ss,se),strand='+' if ss<se else '-',mismatches=mm))
summary=[]
for p in pairs:
 l=hits['p'+str(p['pair'])+'L'];r=hits['p'+str(p['pair'])+'R'];products=[]
 for a in l:
  for b in r:
   if a['chr']!=b['chr'] or a['strand']==b['strand']:continue
   fwd,rev=(a,b) if a['strand']=='+' else (b,a);size=rev['end']-fwd['start']+1
   if fwd['start']<rev['start'] and 50<=size<=5000:products.append(dict(chr=a['chr'],start=fwd['start'],end=rev['end'],size=size,total_mismatches=a['mismatches']+b['mismatches']))
 intended=[x for x in products if x['chr']=='chrX' and x['start']==p['left_start1'] and x['end']==p['right_end1']];other=[x for x in products if x not in intended]
 summary.append(dict(**p,intended_exact_product_found=any(x['total_mismatches']==0 for x in intended),any_other_pair_product=bool(other),any_chrY_pair_product=any(x['chr']=='chrY' for x in other),left_has_other_full_length_hits=any(x['chr']!='chrX' or x['start']!=p['left_start1'] for x in l),right_has_other_full_length_hits=any(x['chr']!='chrX' or x['end']!=p['right_end1'] for x in r),candidate_products=products))
(O/'PUBLIC_ASSAY_FEASIBILITY.json').write_text(json.dumps(summary,indent=2));(O/'ASSAY_FEASIBILITY_COMPLETE').write_text('Public-reference feasibility search completed; no assay validated or ordered. BLAST heuristic may miss matches; no binding-site population-variant screen.');print('ASSAY_FEASIBILITY_COMPLETE',flush=True)
