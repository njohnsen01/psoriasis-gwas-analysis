from pathlib import Path
import time,json,hashlib,subprocess,os
O=Path(__file__).resolve().parent;R=O/'paired_reads';REF=O.parent/'mapping_extension_20260919/reference/Homo_sapiens_assembly38.fasta'
def wait_file(p):
 until=time.monotonic()+7200
 while not p.exists():
  if time.monotonic()>until:raise RuntimeError('Dependency did not complete within the bounded wait')
  time.sleep(20)
wait_file(R/'RECOVERY_VERIFIED.json');meta=json.loads((R/'query_metadata.private.json').read_text());single=R/'same_input_single.private.fastq'
with single.open('w') as out:
 for file,end in [('R1.private.fastq','1'),('R2.private.fastq','2'),('orphans.private.fastq',None)]:
  with (R/file).open() as f:
   while True:
    name=f.readline().rstrip()
    if not name:break
    seq=f.readline().rstrip();plus=f.readline().rstrip();qual=f.readline().rstrip();assert name.startswith('@') and plus=='+' and len(seq)==len(qual);qid=name[1:];m=meta[qid];read_end=end if end else next(iter(m['read_hashes']));assert hashlib.sha256(seq.encode()).hexdigest()==m['read_hashes'][read_end]
    out.write(chr(10).join(['@'+qid+'_e'+read_end,seq,'+',qual])+chr(10))
mm=json.loads((O.parent/'deeper_analysis_20260919/ALIGNER_SOURCE.json').read_text())['binary'];cmd=[mm,'-ax','sr','--secondary=yes','-N','50','-p','0.2','-I','8G','-t','2',str(REF),str(single)]
(R/'same_input_single.command.json').write_text(json.dumps(cmd));print('START_SAME_INPUT_SINGLE_READ_CONTROL',flush=True)
with (R/'same_input_single_minimap.private.sam').open('w') as out,(R/'same_input_single_minimap.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=2400)
(R/'SAME_INPUT_SINGLE_COMPLETE').write_text('Complete; review pending.')
wait_file(O/'BWA_INDEX_READY');bwa=json.loads((O/'BWA_SOURCE.json').read_text())['binary'];prefix=(O/'BWA_INDEX_READY').read_text().strip()
for name,files in [('paired',['R1.private.fastq','R2.private.fastq']),('orphan',['orphans.private.fastq'])]:
 if (R/files[0]).stat().st_size==0:continue
 cmd=[bwa,'mem','-a','-t','2','-K','10000000',prefix]+[str(R/x) for x in files];(R/(name+'_bwa.command.json')).write_text(json.dumps(cmd));print('START_BWA_'+name,flush=True)
 with (R/(name+'_bwa.private.sam')).open('w') as out,(R/(name+'_bwa.log')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,env=dict(os.environ,OMP_NUM_THREADS='2'),check=True,timeout=2400)
(R/'BWA_COMPLETE').write_text('Complete; review pending.');print('BWA_COMPLETE',flush=True)
