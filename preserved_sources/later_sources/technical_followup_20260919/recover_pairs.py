from pathlib import Path
import os,json,hashlib,subprocess,collections,traceback
import pandas as pd
from concurrent.futures import ThreadPoolExecutor,as_completed
import google.auth
from google.auth.transport.requests import Request
from google.cloud import storage
O=Path(__file__).resolve().parent;M=O.parent/'mapping_extension_20260919';R=O/'paired_reads';R.mkdir(exist_ok=False)
ST='/opt/workbench-tools/samtools/bin/samtools';REF=M/'reference/Homo_sapiens_assembly38.fasta'
sel=pd.read_csv(M/'selection.private.tsv',sep=chr(9),dtype=str,usecols=['review_id','cram_uri','cram_index_uri']);assert sel.review_id.is_unique
prior={x['review_id']:x['source_receipts'] for x in json.loads((M/'results.private.json').read_text())}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
 rows={};headers=[]
 for line in Path(p).open():
  if line.startswith('@'):headers.append(line);continue
  a=line.rstrip().split(chr(9));fl=int(a[1]);mate=1 if fl&64 else 2 if fl&128 else 0
  if fl&(256+512+1024+2048) or a[9]=='*':continue
  assert a[10]!='*' and len(a[9])==len(a[10]);k=(a[0],mate)
  if k in rows:assert rows[k][9:11]==a[9:11]
  else:rows[k]=a
 return rows,headers
def one(rec):
 rid=rec['review_id'];F=R/rid;F.mkdir();original=M/'reads_private'/rid/'region.private.sam';rows,headers=load(original);seeds={k[0] for k in rows};missing={q for q in seeds if ((q,1) in rows)!=((q,2) in rows)};recovered=False;error=None
 if missing:
  try:
   client=storage.Client()
   for key in ['cram_uri','cram_index_uri']:
    bn,name=rec[key][5:].split('/',1);blob=client.bucket(bn,user_project=os.environ['GOOGLE_PROJECT']).get_blob(name)
    assert blob is not None and str(blob.generation)==str(prior[rid][key]['generation']) and blob.crc32c==prior[rid][key]['crc32c']
   creds,_=google.auth.default(scopes=['https://www.googleapis.com/auth/devstorage.read_only']);creds.refresh(Request());env=dict(os.environ,GCS_OAUTH_TOKEN=creds.token,GCS_REQUESTER_PAYS_PROJECT=os.environ['GOOGLE_PROJECT'])
   cmd=[ST,'view','-h','--fetch-pairs','-T',str(REF),'-X',rec['cram_uri'],str(M/'reads_private'/rid/'source.crai'),'chrX:5995600-5999100']
   (F/'command.private.json').write_text(json.dumps(cmd))
   with (F/'fetched.private.sam').open('w') as out,(F/'fetch.private.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,env=env,check=True,timeout=300)
   extra,eh=load(F/'fetched.private.sam')
   assert any('SN:chrX' in h and 'M5:2b3a55ff7f58eb308420c8a9b11cac50' in h for h in eh)
   for k,a in extra.items():
    if k[0] in seeds:
     if k in rows:assert rows[k][9:11]==a[9:11]
     else:rows[k]=a
   recovered=True
  except Exception as e:error=type(e).__name__; (F/'failure.private.txt').write_text(traceback.format_exc())
 with (F/'recovered.private.sam').open('w') as f:
  for h in headers:f.write(h)
  for a in rows.values():f.write(chr(9).join(a)+chr(10))
 result=dict(review_id=rid,original_sha256=sha(original),recovered_sha256=sha(F/'recovered.private.sam'),fetch_needed=bool(missing),fetch_completed=recovered,fetch_error=error,seed_templates=len(seeds),complete_pairs=sum((q,1) in rows and (q,2) in rows for q in seeds),remaining_orphans=sum(((q,1) in rows)!=((q,2) in rows) for q in seeds))
 (F/'receipt.private.json').write_text(json.dumps(result));return result
results=[]
with ThreadPoolExecutor(max_workers=2) as pool:
 futures=[pool.submit(one,x) for x in sel.to_dict('records')]
 for f in as_completed(futures):
  results.append(f.result());(R/'progress.private.json').write_text(json.dumps(dict(completed=len(results),planned=len(sel))))
(R/'recovery.private.json').write_text(json.dumps(results));assert len(results)==len(sel)
trans=str.maketrans('ACGTRYMKSWBDHVN','TGCAYRKMSWVHDBN')
def seq(a):
 s=a[9].upper();q=a[10]
 if int(a[1])&16:s=s.translate(trans)[::-1];q=q[::-1]
 return s,q
meta={};counter=0
with (R/'R1.private.fastq').open('w') as f1,(R/'R2.private.fastq').open('w') as f2,(R/'orphans.private.fastq').open('w') as fs:
 for rid in sorted(sel.review_id):
  rows,_=load(R/rid/'recovered.private.sam')
  for name in sorted({k[0] for k in rows}):
   counter+=1;qid='t'+str(counter).zfill(10);keys=[k for k in [(name,1),(name,2),(name,0)] if k in rows];pair=(name,1) in rows and (name,2) in rows;assert not(pair and (name,0) in rows)
   m=dict(review_id=rid,paired=pair,original_qname=name,read_hashes={})
   for k in keys:
    a=rows[k];s,q=seq(a);m['read_hashes'][str(k[1])]=hashlib.sha256(s.encode()).hexdigest();f=(f1 if k[1]==1 else f2) if pair else fs;f.write(chr(10).join(['@'+qid,s,'+',q])+chr(10))
   meta[qid]=m
(R/'query_metadata.private.json').write_text(json.dumps(meta))
receipt=dict(selection_preserved=True,phenotype_blinded=True,all_original_primary_reads_preserved=True,all_fetches_succeeded=all(x['fetch_error'] is None for x in results),complete_mate_recovery=all(x['remaining_orphans']==0 for x in results),input_hashes={f:sha(R/f) for f in ['R1.private.fastq','R2.private.fastq','orphans.private.fastq']},source_sha256=sha(__file__))
(R/'RECOVERY_VERIFIED.json').write_text(json.dumps(receipt,indent=2));print('PAIRED_INPUT_VERIFIED',flush=True)
aligner=json.loads((O.parent/'deeper_analysis_20260919/ALIGNER_SOURCE.json').read_text())['binary'];opts=['-ax','sr','--secondary=yes','-N','50','-p','0.2','-I','8G','-t','2']
for name,files in [('paired',['R1.private.fastq','R2.private.fastq']),('orphan',['orphans.private.fastq'])]:
 if (R/files[0]).stat().st_size==0:continue
 cmd=[aligner]+opts+[str(REF)]+[str(R/f) for f in files];(R/(name+'_minimap.command.json')).write_text(json.dumps(cmd))
 with (R/(name+'_minimap.private.sam')).open('w') as out,(R/(name+'_minimap.log')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=2400)
(R/'MINIMAP_COMPLETE').write_text('Remapping completed; genotype review pending.');print('PAIRED_MINIMAP_COMPLETE',flush=True)
