from pathlib import Path
import sys,json,subprocess,collections,hashlib
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=B/'locus_specific_validation_20260919';R=B/'technical_followup_20260919/paired_reads';M=B/'mapping_extension_20260919';sys.path.insert(0,str(O/'vendor'));sys.path.insert(0,str(M/'source'));import pysam,read_qc_core as core
fa=pysam.FastaFile(str(M/'reference/Homo_sapiens_assembly38.fasta'));P=5998724;lo=P-1001;ref=fa.fetch('chrX',lo,P+1000).upper();assert ref[1000:1003]=='AAG';haps={'reference':ref,'SNP_only':ref[:1000]+'G'+ref[1001:],'MNV_GAT':ref[:1000]+'GAT'+ref[1003:],'neighbor_only':ref[:1002]+'T'+ref[1003:]};rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
plan=dict(purpose='Diagnostic simulation motivated by observed discordance and public rs386823709; not independent replication',reference='same complete GRCh38 as previous alignments',read_length=150,fragment_length=350,target_cycles_zero_based=[15,135],orientations=['forward','reverse'],quality=40,sequencing_errors=0,haplotypes=['reference AAG','SNP-only GAG','MNV GAT','neighbor-only AAT'],parameters='identical to prior paired minimap2 and BWA-MEM2 runs',participant_sequences_used=False)
(O/'SYNTHETIC_CONTROL_PLAN.json').write_text(json.dumps(plan,indent=2));truth={};f1=(O/'synthetic_R1.fastq').open('w');f2=(O/'synthetic_R2.fastq').open('w')
for h,s in haps.items():
 for orient in ['forward','reverse']:
  for cycle in range(15,136):
   target_read_start=1000-cycle if orient=='forward' else 1000-(149-cycle);start=target_read_start if orient=='forward' else target_read_start-200;frag=s[start:start+350];assert len(frag)==350;q=h+'_'+orient+'_'+str(cycle);truth[q]=dict(haplotype=h,expected='ALT' if s[1000]=='G' else 'REF',orientation=orient,cycle=cycle)
   for f,r in [(f1,frag[:150]),(f2,rc(frag[-150:]))]:f.write('@'+q+chr(10)+r+chr(10)+'+'+chr(10)+'I'*150+chr(10))
f1.close();f2.close();(O/'synthetic_truth.json').write_text(json.dumps(truth,indent=2));summary=[]
for model,name in [('minimap2','paired_minimap.command.json'),('bwa_mem2','paired_bwa.command.json')]:
 cmd=json.loads((R/name).read_text());assert str(R/'R1.private.fastq') in cmd and str(R/'R2.private.fastq') in cmd;cmd=[str(O/'synthetic_R1.fastq') if x==str(R/'R1.private.fastq') else str(O/'synthetic_R2.fastq') if x==str(R/'R2.private.fastq') else x for x in cmd];(O/(model+'_synthetic.command.json')).write_text(json.dumps(cmd,indent=2));print('START',model,flush=True)
 with (O/(model+'_synthetic.sam')).open('w') as out,(O/(model+'_synthetic.log')).open('w') as log:subprocess.run(cmd,stdout=out,stderr=log,check=True,timeout=1200)
 obs=collections.defaultdict(list);seen=set()
 for line in (O/(model+'_synthetic.sam')).open():
  if line.startswith('@'):continue
  a=line.rstrip().split(chr(9));flag=int(a[1]);q=a[0]
  if flag&2304:continue
  end='1' if flag&64 else '2';assert (q,end) not in seen;seen.add((q,end));r=core.observe(a,core.TARGETS[2])
  if r:obs[q].append(r)
 assert seen=={(q,e) for q in truth for e in ['1','2']}
 details=[]
 for q,t in truth.items():
  v,conf=core.unique_high_quality(obs[q]);assert len(v)<=1;call=v[0]['call'] if v else 'NO_HQ_CALL';details.append(dict(q=q,**t,observed=call,correct=call==t['expected']))
 (O/(model+'_synthetic_results.json')).write_text(json.dumps(details,indent=2))
 for h in haps:
  v=[x for x in details if x['haplotype']==h];counts=collections.Counter(x['observed'] for x in v);summary.append(dict(model=model,haplotype=h,synthetic_templates=len(v),correct=sum(x['correct'] for x in v),REF=counts['REF'],ALT=counts['ALT'],OTHER=counts['OTHER'],NO_HQ_CALL=counts['NO_HQ_CALL']))
 print('DONE',model,flush=True)
(O/'SYNTHETIC_CONTROL_RESULTS.json').write_text(json.dumps(summary,indent=2));(O/'SYNTHETIC_CONTROLS_VERIFIED.json').write_text(json.dumps(dict(complete_query_accounting=True,identical_previous_parameters=True,public_reference_only=True,independent_genotyping=False,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2));print('SYNTHETIC_CONTROLS_COMPLETE',flush=True)
