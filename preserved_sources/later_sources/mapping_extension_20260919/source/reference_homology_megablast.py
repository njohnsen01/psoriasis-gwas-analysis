from pathlib import Path
import json,subprocess,os,datetime,hashlib,collections,traceback
from mapping_extension import O,sha
R=O/'reference';assert (R/'masked_recovery/FAILED.json').exists() and not (R/'masked_recovery/RUNNING.lock').exists(),'Masked search has not failed or is still active';D=R/'megablast_recovery';D.mkdir(exist_ok=True)
for pid in [3746935,3747226]:
 try:os.kill(pid,0)
 except ProcessLookupError:continue
 raise RuntimeError('Prior reference process still exists; do not duplicate')
assert not (D/'RUNNING.lock').exists() and not (D/'COMPLETE.json').exists()
meta=json.loads((R/'queries.json').read_text());assert len(meta)==18 and len({x['query'] for x in meta})==18
counts=collections.Counter()
for line in (R/'alignments.tsv').open():
 t=line.rstrip().split(chr(9));assert len(t)==13;counts[t[0]]+=1
(R/'FAILED_TIMEOUT.json').write_text(json.dumps(dict(status='failed_timeout',timeout_seconds=1800,partial_output_sha256=sha(R/'alignments.tsv'),represented_queries=dict(counts),absent_queries=[m['query'] for m in meta if m['query'] not in counts],complete=False),indent=2))
(D/'queries.json').write_bytes((R/'queries.json').read_bytes())
command=['blastn','-task','megablast','-query',str(R/'queries.fa'),'-db',str(R/'genome_db'),'-dust','yes','-soft_masking','true','-word_size','28','-evalue','10','-max_target_seqs','10000','-max_hsps','1000','-num_threads','2','-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen','-out',str(D/'alignments.tsv')]
plan=dict(method='MegaBLAST high-similarity context diagnostic, word size 28, DUST seed masking',reason='Unmasked run exceeded its fixed 1800-second bound; preserve partial output and both slower BLASTN runs were incomplete; run a bounded high-similarity MegaBLAST diagnostic',command=command,timeout_seconds=300,queries_sha256=sha(R/'queries.fa'),query_metadata_sha256=sha(R/'queries.json'),reference_database_reused=True,original_unmasked_search_incomplete=True,limitations='Not exhaustive low-complexity homology or actual-read realignment; absence of a hit is conditional on this search.')
(D/'PROTOCOL.json').write_text(json.dumps(plan,indent=2));(D/'RUNNING.lock').write_text(str(os.getpid()))
try:
 with (D/'blastn.log').open('w') as f:r=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,timeout=300)
 assert r.returncode==0
 seen=set();rows=0
 for line in (D/'alignments.tsv').open():
  t=line.rstrip().split(chr(9));assert len(t)==13;seen.add(t[0]);rows+=1
 assert seen=={m['query'] for m in meta}
 (D/'COMPLETE.json').write_text(json.dumps(dict(status='complete',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=command,returncode=r.returncode,queries=18,public_alignment_rows=rows,output_sha256=sha(D/'alignments.tsv'),protocol_sha256=sha(D/'PROTOCOL.json'),original_unmasked_search_incomplete=True),indent=2));print('MEGABLAST_REFERENCE_SEARCH_COMPLETE',flush=True)
except Exception as e:
 (D/'FAILED.json').write_text(json.dumps(dict(status='failed',error_type=type(e).__name__,error=str(e)),indent=2));raise
finally:(D/'RUNNING.lock').unlink(missing_ok=True)

