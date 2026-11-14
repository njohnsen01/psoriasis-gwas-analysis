from pathlib import Path
import sys,json,urllib.request,itertools,hashlib
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=B/'locus_specific_validation_20260919';sys.path.insert(0,str(O/'vendor'));import pysam
url='https://rest.ensembl.org/overlap/region/human/X:5998674-5998774?feature=variation;content-type=application/json';raw=urllib.request.urlopen(urllib.request.Request(url,headers={'Accept':'application/json'}),timeout=60).read();catalog=json.loads(raw);(O/'PUBLIC_ENSEMBL_LOCAL_VARIANTS.json').write_bytes(raw)
fa=pysam.FastaFile(str(B/'mapping_extension_20260919/reference/Homo_sapiens_assembly38.fasta'));P=5998724;rows=json.loads((O/'flank_haplotypes.private.json').read_text());models={};skipped=[]
for d in [20,30]:
 lo=P-1-d;hi=P+d;ref=fa.fetch('chrX',lo,hi).upper();vars=[]
 for x in catalog:
  a=x['start']-1;b=x['end']
  if x.get('assembly_name')!='GRCh38' or not lo<=a<=b<=hi:continue
  alleles=x.get('alleles',[])
  if not alleles:continue
  r=alleles[0].replace('-','')
  if fa.fetch('chrX',a,b).upper()!=r:skipped.append(dict(width=d,id=x['id'],reason='catalog_reference_not_exactly_matching_interval'));continue
  for alt in alleles[1:]:
   alt=alt.replace('-','')
   if set(alt)<=set('ACGT'):vars.append(dict(id=x['id'],start=a-lo,end=b-lo,ref=r,alt=alt))
 haps={ref:[[]]}
 for n in [1,2]:
  for combo in itertools.combinations(vars,n):
   v=sorted(combo,key=lambda x:(x['start'],x['end']))
   if n==2 and (v[0]['id']==v[1]['id'] or v[0]['end']>v[1]['start']):continue
   s=ref
   for x in reversed(v):s=s[:x['start']]+x['alt']+s[x['end']:]
   haps.setdefault(s,[]).append([x['id'] for x in v])
 models[d]=haps
out=[]
for x in rows:
 if not x['previously_flagged']:continue
 for h in x['haplotypes']:
  matches=models[x['inner_half_width']].get(h['haplotype'],[])
  out.append(dict(review_id=x['review_id'],inner_half_width=x['inner_half_width'],template_count=h['template_count'],exact_catalog_matches=matches,any_match=bool(matches),any_match_containing_rs386823709=any('rs386823709' in m for m in matches),all_matches_contain_rs386823709=bool(matches) and all('rs386823709' in m for m in matches)))
(O/'catalog_haplotype_matches.private.json').write_text(json.dumps(out,indent=2));summary=[]
for d in [20,30]:
 z=[x for x in out if x['inner_half_width']==d];supported=[x for x in z if x['template_count']>=3]
 summary.append(dict(inner_half_width=d,any_haplotypes_with_three_templates=bool(supported),all_supported_haplotypes_match_public_catalog=bool(supported) and all(x['any_match'] for x in supported),all_supported_haplotypes_have_MNV_match=bool(supported) and all(x['any_match_containing_rs386823709'] for x in supported),all_exact_explanations_of_supported_haplotypes_include_MNV=bool(supported) and all(x['all_matches_contain_rs386823709'] for x in supported),any_observed_haplotype_unexplained=any(not x['any_match'] for x in z)))
(O/'CATALOG_HAPLOTYPE_AGGREGATE.json').write_text(json.dumps(summary,indent=2));(O/'PUBLIC_CATALOG_REFERENCE_EXCLUSIONS.json').write_text(json.dumps(skipped,indent=2));(O/'CATALOG_CHECK_VERIFIED.json').write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(raw).hexdigest(),reference_checked=True,max_nonoverlapping_catalog_events=2,minimum_templates_for_supported_haplotype=3,limitation='Diagnostic exact matching to a limited public catalog, not a genotype or phase call; support may be insufficient per participant',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2));print('CATALOG_HAPLOTYPE_CHECK_COMPLETE',flush=True)
