from pathlib import Path
import sys,json,hashlib,subprocess,datetime,collections,re,statistics
B=Path(__file__).resolve().parent.parent;O=B/'deeper_analysis_20260919';M=B/'mapping_extension_20260919';R=O/'read_remapping';R.mkdir(exist_ok=False)
sys.path.insert(0,str(M/'source'));from read_qc_core import TARGETS,observe,unique_high_quality,self_test
assert self_test();ref=M/'reference/Homo_sapiens_assembly38.fasta';aligner=json.loads((O/'ALIGNER_SOURCE.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def rc(s):return s.translate(str.maketrans('ACGTRYMKSWBDHVN','TGCAYRKMSWVHDBN'))[::-1]
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selection':'Exactly the retained bounded sample; one original HQ representative per template and target, original MQ >=30 and BQ >=20; both REF and ALT and OTHER retained without phenotype selection','reference':str(ref),'reference_sha256':sha(ref),'source_sha256':sha(Path(__file__)),'helper_sha256':sha(M/'source/read_qc_core.py'),'aligner':aligner,'method':'Whole-reference independent single-read minimap2 sr realignment; not paired-read recovery or orthogonal sequencing','options':['-ax','sr','--secondary=yes','-N','50','-p','0.2','-I','8G','-t','2'],'expected_region':'chrX:5995000-6000000','limitations':'Originally chrX-selected reads only; cannot discover reads excluded by original extraction; heuristic search with up to 50 secondary hits; alternate contigs retained and not automatically treated as distinct loci; not whole-cohort certification','privacy':'All sequences, query metadata and per-participant results private inside Workbench. Only screened aggregate summaries may be reported.'}
(R/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
proc=subprocess.Popen(['samtools','faidx',str(ref),'chrX'],stdout=subprocess.PIPE);h=hashlib.md5()
for line in proc.stdout:
 if not line.startswith(b'>'):h.update(line.strip().upper())
assert proc.wait()==0 and h.hexdigest()=='2b3a55ff7f58eb308420c8a9b11cac50'
s=subprocess.check_output(['samtools','faidx',str(ref),'chrX:5996489-5996506'],text=True);seq=''.join(s.splitlines()[1:]).upper();assert seq[5:8]=='GCT';window=(5996489,5996506,seq,seq[:5]+'G'+seq[8:])
prior=json.loads((M/'results.private.json').read_text());queries={};seqs={};counter=0
for person in prior:
 folder=M/'reads_private'/person['review_id'];rows=[[],[],[]]
 for line in (folder/'region.private.sam').read_text().splitlines():
  if line.startswith('@'):continue
  fields=line.split(chr(9));counter+=1;qid='q'+str(counter).zfill(9)
  for j,target in enumerate(TARGETS):
   x=observe(fields,target,window if j==1 else None)
   if x is not None:x['_fields']=fields;x['_qid']=qid;rows[j].append(x)
 for j,target in enumerate(TARGETS):
  chosen,conflicts=unique_high_quality(rows[j]);old=next(t for t in person['targets'] if t['target']=='g'+str(j));assert len(chosen)==old['HQ_template_depth']
  for x in chosen:
   qid=x['_qid'];f=x['_fields'];seq=f[9].upper();qual=f[10];assert qual!='*' and len(seq)==len(qual)
   if int(f[1])&16:seq=rc(seq);qual=qual[::-1]
   meta=queries.setdefault(qid,{'review_id':person['review_id'],'sex':person['sex'],'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'annotations':[]})
   meta['annotations'].append({'target':'g'+str(j),'old_call':x['call'],'old_MQ':x['MQ'],'old_BQ':x['BQ'],'expected_fraction':old['expected_target_fraction']});seqs[qid]=(seq,qual)
with (R/'queries.private.fastq').open('w') as f:
 for qid,(seq,qual) in seqs.items():f.write(chr(10).join(['@'+qid,seq,'+',qual])+chr(10))
(R/'query_metadata.private.json').write_text(json.dumps(queries));del seqs
print('READ_SELECTION_AND_REFERENCE_VERIFIED',flush=True)
command=[aligner['binary']]+protocol['options']+[str(ref),str(R/'queries.private.fastq')];(R/'COMMAND.json').write_text(json.dumps(command,indent=2))
with (R/'remapped.private.sam').open('w') as out,(R/'aligner.log').open('w') as err:
 subprocess.run(command,stdout=out,stderr=err,check=True,timeout=1800)
primary={};secondary=collections.defaultdict(list)
with (R/'remapped.private.sam').open() as f:
 for line in f:
  if line.startswith('@'):continue
  a=line.rstrip().split(chr(9));qid=a[0];assert qid in queries;flag=int(a[1]);tags={x.split(':',2)[0]:x.split(':',2)[2] for x in a[11:] if len(x.split(':',2))==3}
  if flag&(256|2048):secondary[qid].append({'rname':a[2],'pos':int(a[3]),'AS':int(tags.get('AS','-99999'))});continue
  assert qid not in primary;seq=a[9].upper();canonical=rc(seq) if flag&16 else seq;assert hashlib.sha256(canonical.encode()).hexdigest()==queries[qid]['sequence_sha256'];primary[qid]=a
assert set(primary)==set(queries)
records=[]
for qid,a in primary.items():
 meta=queries[qid];flag=int(a[1]);expected=not flag&4 and a[2]=='chrX' and 5995000<=int(a[3])<=6000000
 for annotation in meta['annotations']:
  j=int(annotation['target'][1:]);x=observe(a,TARGETS[j],window if j==1 else None);newcall=x['call'] if x else None
  records.append(dict(annotation,review_id=meta['review_id'],sex=meta['sex'],primary_rname=a[2],primary_expected_X=bool(expected),new_MQ=int(a[4]),new_call=newcall,allele_concordant=bool(expected and newcall==annotation['old_call']),HQ_concordant=bool(expected and newcall==annotation['old_call'] and int(a[4])>=30 and x['BQ']>=20),secondary_chrY=any(x['rname']=='chrY' for x in secondary[qid]),secondary_limit_reached=len(secondary[qid])>=50))
(R/'read_results.private.json').write_text(json.dumps(records))
groups=collections.defaultdict(list)
for x in records:groups[(x['target'],x['sex'],x['old_call'])].append(x)
summaries=[]
for key,rs in sorted(groups.items()):
 people=collections.defaultdict(list)
 for x in rs:people[x['review_id']].append(x)
 entry={'target':key[0],'sex':'male' if key[1]==1 else 'female','original_allele':key[2]}
 if len(people)<=20:entry['status']='not_reportable'
 else:
  entry['status']='reportable';entry['participant_median_percent']={metric:round(100*statistics.median([sum(bool(x[metric]) for x in vals)/len(vals) for vals in people.values()]),1) for metric in ['primary_expected_X','allele_concordant','HQ_concordant','secondary_chrY','secondary_limit_reached']}
 summaries.append(entry)
(R/'AGGREGATE_REVIEW.json').write_text(json.dumps(summaries,indent=2));audit={'status':'bounded_whole_reference_single_read_remapping_completed','reference_chrX_md5_verified':True,'original_HQ_template_selection_reproduced':True,'every_query_accounted_for_once_as_primary_or_unmapped':True,'every_read_sequence_unchanged_after_orientation_check':True,'raw_data_private':True,'full_cohort_or_orthogonal_validation':False,'reference_sha256':protocol['reference_sha256'],'source_sha256':sha(Path(__file__)),'input_fastq_sha256':sha(R/'queries.private.fastq'),'output_sam_sha256':sha(R/'remapped.private.sam')};(R/'VERIFIED.json').write_text(json.dumps(audit,indent=2));print('READ_REMAPPING_VERIFIED',flush=True)
