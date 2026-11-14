from pathlib import Path
import collections,hashlib,json,sys,re,statistics,math
import pandas as pd
from scipy.stats import mannwhitneyu
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912')
M=B/'mapping_extension_20260919';D=B/'deeper_analysis_20260919';O=B/'technical_followup_20260919';R=O/'paired_reads';REF=M/'reference/Homo_sapiens_assembly38.fasta'
sys.path.insert(0,str(M/'source'));import read_qc_core as core
assert core.self_test()
for marker in ['MINIMAP_COMPLETE','SAME_INPUT_SINGLE_COMPLETE','BWA_COMPLETE']:assert (R/marker).exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for data in iter(lambda:f.read(8388608),b''):h.update(data)
 return h.hexdigest()
assert sha(REF)==json.loads((D/'read_remapping/PROTOCOL.json').read_text())['reference_sha256']
fai={a[0]:list(map(int,a[1:])) for a in [s.split() for s in Path(str(REF)+'.fai').read_text().splitlines()]}
length,offset,width,bytes_per_line=fai['chrX'];lo,hi=5996489,5996506
with REF.open('rb') as f:
 seq=''
 for pos in range(lo-1,hi):
  f.seek(offset+(pos//width)*bytes_per_line+pos%width);seq+=f.read(1).decode().upper()
assert seq[5:8]=='GCT';window=(lo,hi,seq,seq[:5]+'G'+seq[8:])
meta=json.loads((R/'query_metadata.private.json').read_text());prior=json.loads((M/'results.private.json').read_text());people={x['review_id']:x for x in prior}
expected_keys={(q,e) for q,m in meta.items() for e in m['read_hashes']};trans=str.maketrans('ACGTRYMKSWBDHVN','TGCAYRKMSWVHDBN')
checks={};alignment_alternatives={}
def collect(files,mode,rid=None):
 buckets=collections.defaultdict(list);seen=set();namespaces=collections.defaultdict(set);y_secondary=False;primary_y=False
 for path in files:
  with Path(path).open() as f:
   for line in f:
    if line.startswith('@'):continue
    a=line.rstrip().split(chr(9));flag=int(a[1])
    if flag&256 and a[2]=='chrY':y_secondary=True
    if flag&2304:continue
    if mode in ['paired','single']:
     if mode=='single':q,end=a[0].rsplit('_e',1)
     else:q=a[0];end=('1' if flag&64 else '2') if meta[q]['paired'] else next(iter(meta[q]['read_hashes']))
     assert (q,end) in expected_keys and (q,end) not in seen;seen.add((q,end));s=a[9].upper();assert s!='*'
     if flag&16:s=s.translate(trans)[::-1]
     assert hashlib.sha256(s.encode()).hexdigest()==meta[q]['read_hashes'][end]
     person=meta[q]['review_id'];a[0]=q
    else:
     person=rid;tags={z.split(':',2)[0]:z.split(':',2)[2] for z in a[11:] if len(z.split(':',2))==3};namespaces[a[0]].add(tags.get('RG',''))
    if not flag&4 and a[2]=='chrY':primary_y=True
    if flag&(4+256+512+1024+2048) or a[2]!='chrX' or a[9]=='*':continue
    start=int(a[3]);endpos=start+sum(int(n) for n,op in re.findall('('+chr(92)+'d+)([MIDNSHP=X])',a[5]) if op in 'MDN=X')-1
    for j,t in enumerate(core.TARGETS):
     if start<=t[0]<=endpos:
      v=core.observe(a,t,window if j==1 else None)
      if v:buckets[(person,'g'+str(j))].append(v)
 if mode in ['paired','single']:assert seen==expected_keys
 else:assert all(len(v)==1 for v in namespaces.values())
 return buckets,dict(query_accounting_verified=mode in ['paired','single'],sequence_hashes_verified=mode in ['paired','single'],any_primary_chrY=primary_y,any_secondary_chrY=y_secondary)
def balance(frac,expected):
 return frac is not None and (frac<=.1 if expected==0 else frac>=.9 if expected==1 else .25<=frac<=.75)
def summarize(model,buckets):
 out=[]
 for rid,p in people.items():
  for old in p['targets']:
   rows=buckets.get((rid,old['target']),[]);hq,conflicts=core.unique_high_quality(rows);counts=collections.Counter(x['call'] for x in hq);depth=len(hq);exact=counts['REF']+counts['ALT'];frac=counts['ALT']/depth if depth else None;efrac=counts['ALT']/exact if exact else None
   x=dict(model=model,review_id=rid,target=old['target'],sex=p['sex'],ancestry=p['ancestry'],case_status=p['case_status'],expected=old['expected_target_fraction'],HQ_depth=depth,exact_depth=exact,ALT_fraction=frac,exact_ALT_fraction=efrac,adequate=depth>=10,balance_pass=balance(frac,old['expected_target_fraction']),exact_balance_pass=balance(efrac,old['expected_target_fraction']),counts=dict(counts),conflicts=conflicts,MQ_median=statistics.median(z['MQ'] for z in rows) if rows else None,low_MQ_fraction=sum(z['MQ']<30 for z in rows)/len(rows) if rows else None,mean_soft_clip=statistics.mean(z['soft_clipped'] for z in rows) if rows else None,allele_strand={c:{'forward':sum(z['call']==c and not z['reverse'] for z in hq),'reverse':sum(z['call']==c and z['reverse'] for z in hq)} for c in ['REF','ALT','OTHER']})
   if model=='original':
    assert depth==old['HQ_template_depth'] and dict(counts)==old['high_quality_allele_counts'];assert x['balance_pass']==old['diagnostic_genotype_concordance']
   out.append(x)
 return out
allrows=[]
for model,filename in [('original','region.private.sam'),('recovered_original','recovered.private.sam')]:
 b=collections.defaultdict(list)
 for rid in people:
  path=(M/'reads_private'/rid/filename) if model=='original' else (R/rid/filename)
  if not path.exists():path=M/'reads_private'/rid/'region.private.sam'
  d,c=collect([path],'original',rid)
  for k,v in d.items():b[k].extend(v)
 allrows.extend(summarize(model,b));print('VERIFIED',model,flush=True)
for model,mode,files in [('minimap_paired','paired',['paired_minimap.private.sam','orphan_minimap.private.sam']),('minimap_same_input_single','single',['same_input_single_minimap.private.sam']),('bwa_paired','paired',['paired_bwa.private.sam','orphan_bwa.private.sam'])]:
 b,c=collect([R/f for f in files if (R/f).exists()],mode);checks[model]=c;allrows.extend(summarize(model,b));print('VERIFIED',model,flush=True)
(R/'paired_genotype_review.private.json').write_text(json.dumps(allrows,indent=2))
a=pd.DataFrame(allrows);agg=[]
for (model,target),z in a.groupby(['model','target']):
 agg.append(dict(model=model,target=target,any_adequate_depth_balance_flag=bool((z.adequate & ~z.balance_pass).any()),any_adequate_exact_depth_balance_flag=bool(((z.exact_depth>=10)&~z.exact_balance_pass).any()),any_inadequate_HQ_depth=bool((~z.adequate).any()),median_HQ_depth=float(z.HQ_depth.median()),median_MQ=float(z.MQ_median.median()),median_low_MQ_fraction=float(z.low_MQ_fraction.median())))
(R/'PAIRED_REVIEW_AGGREGATE.json').write_text(json.dumps(agg,indent=2))
qc=[]
for (model,target,sex,ancestry,expected),z in a[a.model.isin(['original','minimap_paired','bwa_paired'])].groupby(['model','target','sex','ancestry','expected']):
 for metric in ['HQ_depth','MQ_median','low_MQ_fraction','mean_soft_clip']:
  case=z.loc[z.case_status==1,metric].dropna();control=z.loc[z.case_status==0,metric].dropna();row=dict(model=model,target=target,sex=int(sex),ancestry=ancestry,genotype_fraction=expected,metric=metric,status='support_gate_not_met')
  if len(case)>=21 and len(control)>=21:
   u,p=mannwhitneyu(case,control,alternative='two-sided');row.update(status='verified',case_median=float(case.median()),control_median=float(control.median()),p=float(p),rank_biserial=float(2*u/(len(case)*len(control))-1))
  qc.append(row)
q=pd.DataFrame(qc);ok=q.index[q.status=='verified'];order=q.loc[ok,'p'].sort_values().index
if len(order):
 vals=q.loc[order,'p'].values*len(order)/range(1,len(order)+1);q.loc[order,'BH_q']=__import__('numpy').minimum.accumulate(vals[::-1])[::-1].clip(0,1)
q.to_csv(R/'READ_QUALITY_CASE_CONTROL.csv',index=False)
receipt=dict(reference_unchanged=True,original_pileups_reproduced=True,read_group_namespaces_unambiguous=True,same_input_single_pair_control=True,alignments=checks,source_sha256=sha(__file__),independent_sequencing=False)
(R/'PAIRED_REVIEW_VERIFIED.json').write_text(json.dumps(receipt,indent=2));print('PAIRED_REVIEW_COMPLETE',flush=True)
