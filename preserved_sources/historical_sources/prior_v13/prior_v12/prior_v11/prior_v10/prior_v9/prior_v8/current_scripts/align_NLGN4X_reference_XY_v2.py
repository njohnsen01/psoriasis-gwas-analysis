from pathlib import Path
import json,hashlib,subprocess,traceback,datetime,os
from google.cloud import storage
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');O=N/'reconciliation_20260912/read_mapping_QC';F=O/'broad_reference';F.mkdir(exist_ok=True)
try:
 client=storage.Client();bucket=client.bucket('gcp-public-data--broad-references',user_project=os.environ['GOOGLE_PROJECT']);name='hg38/v0/Homo_sapiens_assembly38.fasta';b=bucket.get_blob(name);assert b is not None
 fai=bucket.blob(name+'.fai').download_as_text();dictionary=bucket.blob('hg38/v0/Homo_sapiens_assembly38.dict').download_as_text();(F/'source.fai').write_text(fai);(F/'source.dict').write_text(dictionary)
 idx={x[0]:list(map(int,x[1:])) for line in fai.splitlines() for x in [line.split(chr(9))]};dic={}
 for line in dictionary.splitlines():
  if line.startswith('@SQ'):
   r=dict(x.split(':',1) for x in line.split(chr(9))[1:]);dic[r['SN']]=r
 seqs={};sources=[]
 for chrom in ['chrX','chrY']:
  length,offset,linebases,linewidth=idx[chrom];end=offset+((length-1)//linebases)*linewidth+(length-1)%linebases
  seq=b.download_as_bytes(start=offset,end=end,if_generation_match=b.generation,timeout=180).decode().replace(chr(10),'').replace(chr(13),'').upper();assert len(seq)==length;md5=hashlib.md5(seq.encode()).hexdigest();assert md5==dic[chrom]['M5'];seqs[chrom]=seq
  (F/(chrom+'.fa')).write_text('>'+chrom+chr(10)+chr(10).join(seq[i:i+60] for i in range(0,len(seq),60))+chr(10));sources.append(dict(uri='gs://'+bucket.name+'/'+name,generation=b.generation,chromosome=chrom,byte_start=offset,byte_end=end,sequence_md5=md5,length=length))
 assert seqs['chrX'][5998724-1]=='A';assert sources[0]['sequence_md5']=='2b3a55ff7f58eb308420c8a9b11cac50','CRAM reference M5 concordance'
 (F/'reference_sources.json').write_text(json.dumps(sources,indent=2))
 with (F/'XY.fa').open('w') as out:
  for chrom in ['chrX','chrY']:out.write((F/(chrom+'.fa')).read_text())
 queries=[]
 for length in [101,151,251,501]:
  for offset in ([0] if length!=151 else [-50,-25,0,25,50]):
   start=5998724-1-length//2+offset;sequence=seqs['chrX'][start:start+length];queries.append((f'NLGN4X_len{length}_offset{offset}',sequence,5998724-start))
 (F/'lead_queries.fa').write_text(''.join('>'+name+chr(10)+seq+chr(10) for name,seq,center in queries));(F/'query_centers.json').write_text(json.dumps({name:center for name,seq,center in queries},indent=2))
 commands=[['samtools','faidx',str(F/'chrX.fa')],['makeblastdb','-in',str(F/'XY.fa'),'-dbtype','nucl','-out',str(F/'XY_db')],['blastn','-task','blastn','-query',str(F/'lead_queries.fa'),'-db',str(F/'XY_db'),'-dust','no','-soft_masking','false','-word_size','11','-evalue','10','-max_target_seqs','1000','-max_hsps','1000','-num_threads','2','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq','-out',str(F/'XY_alignments.tsv')]]
 for k,cmd in enumerate(commands):
  with (F/f'command_{k}.log').open('w') as out:r=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
  assert r.returncode==0,(k,r.returncode)
 (F/'commands.json').write_text(json.dumps(commands,indent=2));(F/'COMPLETE.json').write_text(json.dumps(dict(status='reference_alignment_completed',queries=len(queries),exact_reference_M5_matches_CRAM=True,scope='GRCh38 primary chrX and chrY reference sequence only; not whole-genome mappability, read realignment or proof against an artifact.',date_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2));print('XY_REFERENCE_ALIGNMENT_COMPLETE',flush=True)
except Exception as e:(F/'FAILED.json').write_text(json.dumps(dict(error=repr(e),traceback=traceback.format_exc()),indent=2));raise
