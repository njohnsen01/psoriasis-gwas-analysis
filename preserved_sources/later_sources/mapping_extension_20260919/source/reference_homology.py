from pathlib import Path
import json,hashlib,subprocess,os,datetime
from google.cloud import storage
from read_qc_core import TARGETS
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=B/'mapping_extension_20260919/reference';O.mkdir(parents=True,exist_ok=True)
def run(args,name,timeout):
 with (O/(name+'.log')).open('w') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
 assert p.returncode==0,name+' failed; inspect saved log'
def main():
 lock=O/'RUNNING.lock';assert not lock.exists() and not (O/'COMPLETE.json').exists();lock.write_text(str(os.getpid()))
 try:
  sources=json.loads((B/'read_mapping_QC/broad_reference/reference_sources.json').read_text());src=sources[0];bn,key=src['uri'][5:].split('/',1)
  client=storage.Client();blob=client.bucket(bn).get_blob(key,generation=int(src['generation']));assert blob
  fa=O/'Homo_sapiens_assembly38.fasta';partial=O/'reference.partial'
  blob.download_to_filename(str(partial),if_generation_match=blob.generation,checksum='crc32c',timeout=900);assert partial.stat().st_size==blob.size;partial.replace(fa)
  run(['samtools','faidx',str(fa)],'faidx',300)
  xx=subprocess.run(['samtools','faidx',str(fa),'chrX'],capture_output=True,check=True).stdout.splitlines()[1:];seq=b''.join(xx);md5=hashlib.md5(seq).hexdigest();assert md5=='2b3a55ff7f58eb308420c8a9b11cac50';del xx,seq
  queries=[];meta=[]
  for j,(pos,ref,alt) in enumerate(TARGETS):
   for flank in [50,150,500]:
    start,end=pos-flank,pos+len(ref)-1+flank;p=subprocess.run(['samtools','faidx',str(fa),f'chrX:{start}-{end}'],capture_output=True,text=True,check=True);seq=''.join(p.stdout.splitlines()[1:]).upper();assert seq[flank:flank+len(ref)]==ref
    for label,allele in [('REF',ref),('ALT',alt)]:
     qid=f'g{j}_{label}_{flank}';qs=seq[:flank]+allele+seq[flank+len(ref):];queries.append('>'+qid+'\n'+qs+'\n');meta.append(dict(query=qid,target=j,allele=label,flank=flank,start=start,end=end,center=flank+1,length=len(qs)))
  (O/'queries.fa').write_text(''.join(queries));(O/'queries.json').write_text(json.dumps(meta,indent=2))
  run(['makeblastdb','-in',str(fa),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(O/'genome_db')],'makeblastdb',1800)
  command=['blastn','-task','blastn','-query',str(O/'queries.fa'),'-db',str(O/'genome_db'),'-dust','no','-soft_masking','false','-word_size','11','-evalue','10','-max_target_seqs','10000','-max_hsps','1000','-num_threads','2','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen','-out',str(O/'alignments.tsv')]
  run(command,'blastn',1800)
  (O/'COMPLETE.json').write_text(json.dumps(dict(status='complete',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),public_source=src['uri'],generation=blob.generation,crc32c=blob.crc32c,bytes=blob.size,chrX_MD5=md5,queries=len(meta),command=command,interpretation='Reference/ALT homology diagnostic, not independent read realignment, formal mappability, or proof all reads map correctly.'),indent=2));print('PUBLIC_REFERENCE_HOMOLOGY_COMPLETE',flush=True)
 finally:lock.unlink(missing_ok=True)
if __name__=='__main__':main()

