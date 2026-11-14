from pathlib import Path
import os,sys,json,hashlib,subprocess,datetime,traceback,re
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import Counter
import numpy as np
import pandas as pd
from google.cloud import storage
import google.auth
from google.auth.transport.requests import Request
from read_qc_core import TARGETS,EXCLUDE,cigar_map,observe,unique_high_quality,self_test
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
B=N/'reconciliation_20260912';O=B/'mapping_extension_20260919';OLD=B/'read_mapping_QC';Q=N.parent/'sample_lead_qc_20260910';PQ=B/'joint_G_sex_20260914/provider_quality';REF=OLD/'broad_reference/chrX.fa';ENTRY=PQ/'three_sites_entries.private.tsv.bgz'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
def write(name,obj):
 p=O/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2,default=str,allow_nan=False));t.replace(p)
def cmd(args,logfile,**kw):
 with Path(logfile).open('w') as err:p=subprocess.run(args,stderr=err,**kw)
 if p.returncode:raise RuntimeError('command failed; inspect private log '+Path(logfile).name)
 return p
def prepare():
 O.mkdir(exist_ok=True);assert not (O/'selection.private.tsv').exists();assert self_test()
 protocol=dict(created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),targets=TARGETS,purpose='Read/mapping diagnostic, not association reanalysis or independent genotype confirmation',selection='Strict-eligible source cohort with called target genotypes, EUR/AFR/AMR; up to 25 per ancestry/sex/case-control/any-target-carrier stratum, ordered by SHA256 seed 20260919. Preserve eligible original pilot participants. No GQ/read-quality/association-P selection.',maximum_selected=648,max_workers=2,region='chrX:5995600-5999100',read_gate='Primary nonduplicate non-QC-failed nonsupplementary mapped reads; HQ MQ>=30 and <255; BQ>=20; overlapping templates counted once by BQ then MQ; tied conflicting calls removed.',indel='Span X:5996489-5996506. Compare reconstructed local haplotype against reference and exact deletion; equivalent shifted CIGAR deletions accepted. Reads lacking full span excluded from indel fraction.',diagnostic_flags='HQ template depth <10; target fraction >0.10 for noncarrier, <0.90 for target homozygous/hemizygous, outside 0.25-0.75 for heterozygote; GQ<20. Flags do not change GWAS samples/genotypes.',X_Y='Primary read MQ, mate/alternative chrY tags, reference/ALT context homology. Absence of tags cannot exclude alternative alignment.',limitations='Bounded stratified sample, not full-cohort read validation or proof of absence of mapping artifacts.',release='All participant/read records private; only reviewed disclosure-safe aggregates may leave Workbench')
 write('PROTOCOL.json',protocol)
 d=pd.read_csv(Q/'sample_audit.private.tsv',sep='\t',dtype={'IID':str});d=d.loc[d.strict_eligible.eq(True)].copy();assert d.IID.is_unique
 valid=set(d.IID);chunks=[]
 for e in pd.read_csv(ENTRY,sep='\t',compression='gzip',dtype={'s':str,'GT':str},chunksize=100000):
  e=e.loc[e.s.isin(valid)].copy()
  for j,(pos,ref,alt) in enumerate(TARGETS):
   z=e.loc[e.locus.eq('chrX:'+str(pos))].copy()
   if z.empty:continue
   def parse(r):
    alleles=json.loads(r.alleles);assert alleles[0]==ref and alt in alleles
    if not isinstance(r.GT,str) or not re.fullmatch(r'\d+([/|]\d+)?',r.GT):return (np.nan,np.nan)
    gt=list(map(int,re.split(r'[/|]',r.GT)));return (sum(x==alleles.index(alt) for x in gt)/len(gt),len(set(gt)))
   z[['fraction','distinct']]=pd.DataFrame([parse(r) for r in z.itertuples()],index=z.index);z['target']='g'+str(j);chunks.append(z[['s','target','fraction','distinct','GQ']])
 e=pd.concat(chunks,ignore_index=True);assert not e.duplicated(['s','target']).any()
 f=e.pivot(index='s',columns='target',values='fraction');ng=e.pivot(index='s',columns='target',values='distinct');gq=e.pivot(index='s',columns='target',values='GQ').add_prefix('GQ_')
 z=d.merge(f,left_on='IID',right_index=True,validate='one_to_one').merge(gq,left_on='IID',right_index=True,validate='one_to_one')
 for j in range(3):
  bad=z.SEX.eq(1)&z.IID.map(ng['g'+str(j)]).gt(1);z.loc[bad,'g'+str(j)]=np.nan
 z=z.loc[z[['g0','g1','g2']].notna().all(axis=1)].copy();z['carrier']=z[['g0','g1','g2']].gt(0).any(axis=1);z['ancestry_pred']=z.ancestry_pred.str.upper()
 cm=pd.read_csv(OLD/'cram_manifest.private.csv',dtype={'person_id':str});z=z.merge(cm,left_on='IID',right_on='person_id',validate='one_to_one')
 z=z.loc[z.ancestry_pred.isin(['EUR','AFR','AMR'])&z.SEX.isin([1,2])&z.psoriasis_case.isin([0,1])].copy();z['selection_hash']=z.IID.map(lambda s:hashlib.sha256(('20260919:'+s).encode()).hexdigest())
 chosen=z.sort_values('selection_hash').groupby(['ancestry_pred','SEX','psoriasis_case','carrier'],group_keys=False).head(25)
 old=pd.read_csv(OLD/'CRAM_pilot_selection.private.tsv',sep='\t',dtype={'IID':str});chosen=pd.concat([chosen,z.loc[z.IID.isin(old.IID)]]).drop_duplicates('IID').sort_values('selection_hash').reset_index(drop=True)
 assert 20<len(chosen)<=648
 chosen['review_id']=['m'+str(i+1).zfill(4) for i in range(len(chosen))];chosen.to_csv(O/'selection.private.tsv',sep='\t',index=False)
 write('INPUT_RECEIPTS.json',{str(p):sha(p) for p in [ENTRY,Q/'sample_audit.private.tsv',OLD/'cram_manifest.private.csv',OLD/'CRAM_pilot_selection.private.tsv',Path(__file__),Path(__file__).with_name('read_qc_core.py')]})
 write('PREPARED.json',dict(status='prepared',selection_sha256=sha(O/'selection.private.tsv'),protocol_sha256=sha(O/'PROTOCOL.json'),association_models_unchanged=True));print('PREPARED_FIXED_SELECTION_AND_PROTOCOL',flush=True)
def reference_window():
 p=subprocess.run(['samtools','faidx',str(REF),'chrX:5996489-5996506'],capture_output=True,text=True,check=True);seq=''.join(p.stdout.splitlines()[1:]).upper();assert seq[5:8]=='GCT'
 return (5996489,5996506,seq,seq[:5]+'G'+seq[8:])
def one(rec,window):
 folder=O/'reads_private'/rec['review_id'];folder.mkdir(parents=True,exist_ok=True);done=folder/'COMPLETE.private.json'
 if done.exists():return json.loads(done.read_text())
 creds,_=google.auth.default(scopes=['https://www.googleapis.com/auth/devstorage.read_only']);creds.refresh(Request());env=dict(os.environ,GCS_OAUTH_TOKEN=creds.token,GCS_REQUESTER_PAYS_PROJECT=os.environ['GOOGLE_PROJECT']);client=storage.Client();receipts={}
 for key in ['cram_uri','cram_index_uri']:
  bn,name=rec[key][5:].split('/',1);blob=client.bucket(bn,user_project=os.environ['GOOGLE_PROJECT']).get_blob(name);assert blob is not None
  receipts[key]=dict(uri=rec[key],generation=blob.generation,size=blob.size,crc32c=blob.crc32c)
  if key=='cram_index_uri':blob.download_to_filename(str(folder/'source.crai'),if_generation_match=blob.generation,checksum='crc32c',timeout=90)
 command=['samtools','view','-h','-T',str(REF),'-X',rec['cram_uri'],str(folder/'source.crai'),'chrX:5995600-5999100']
 with (folder/'region.private.sam').open('w') as out:cmd(command,folder/'extract.private.log',stdout=out,env=env,timeout=180)
 rows=[[],[],[]];raw_depth=[0,0,0];rgs=[];xmd5=None
 for line in (folder/'region.private.sam').read_text().splitlines():
  t=line.split('\t')
  if line.startswith('@SQ'):
   x=dict(v.split(':',1) for v in t[1:])
   if x.get('SN')=='chrX':xmd5=x.get('M5')
  if line.startswith('@RG'):rgs.append(dict(v.split(':',1) for v in t[1:]))
  if line.startswith('@'):continue
  if not int(t[1])&EXCLUDE and t[2]=='chrX' and t[9]!='*':
   bases,_,_=cigar_map(int(t[3]),t[5],t[9],t[10])
   for j,(pos,_,_) in enumerate(TARGETS):raw_depth[j]+=int(pos in bases)
  for j,target in enumerate(TARGETS):
   x=observe(t,target,window if j==1 else None)
   if x is not None:rows[j].append(x)
 assert xmd5=='2b3a55ff7f58eb308420c8a9b11cac50','reference MD5 mismatch'
 bam=folder/'region.private.bam';cmd(['samtools','view','-b','-o',str(bam),str(folder/'region.private.sam')],folder/'bam.private.log',stdout=subprocess.DEVNULL,timeout=60);cmd(['samtools','index',str(bam)],folder/'index.private.log',stdout=subprocess.DEVNULL,timeout=60)
 result=[]
 for j,(pos,ref,alt) in enumerate(TARGETS):
  depth=cmd(['samtools','depth','-r',f'chrX:{pos}-{pos}','-q','0','-Q','0','-G',str(EXCLUDE),str(bam)],folder/f'depth_g{j}.private.log',stdout=subprocess.PIPE,text=True,timeout=30).stdout
  observed=0 if not depth.strip() else int(depth.strip().split('\t')[-1]);assert observed==raw_depth[j],'independent depth check failed'
  high,conflicts=unique_high_quality(rows[j]);counts=Counter(r['call'] for r in high);frac=counts['ALT']/len(high) if high else None;expected=float(rec['g'+str(j)]);mq=[r['MQ'] for r in rows[j] if r['MQ']<255]
  concordant=(frac is not None and (frac<=.10 if expected==0 else frac>=.90 if expected==1 else .25<=frac<=.75))
  v=dict(target='g'+str(j),position=pos,expected_target_fraction=expected,GQ=None if pd.isna(rec['GQ_g'+str(j)]) else float(rec['GQ_g'+str(j)]),aligned_base_read_depth=raw_depth[j],allele_informative_read_depth=len(rows[j]),HQ_template_depth=len(high),HQ_target_fraction=frac,MQ_median=float(np.median(mq)) if mq else None,low_MQ_fraction=sum(r['MQ']<30 for r in rows[j])/len(rows[j]) if rows[j] else None,MQ_unknown_fraction=sum(r['MQ']==255 for r in rows[j])/len(rows[j]) if rows[j] else None,mean_soft_clip=float(np.mean([r['soft_clipped'] for r in rows[j]])) if rows[j] else None,high_quality_allele_counts=dict(counts),overlap_conflicts=conflicts,mate_chrY_reads=sum(r['mate_chrY'] for r in rows[j]),alternative_chrY_reads=sum(r['alternative_chrY'] for r in rows[j]),depth_crosscheck_passed=True,adequate_HQ_depth=len(high)>=10,diagnostic_genotype_concordance=bool(concordant),allele_strand={a:dict(forward=sum(r['call']==a and not r['reverse'] for r in high),reverse=sum(r['call']==a and r['reverse'] for r in high)) for a in ['REF','ALT','OTHER']})
  result.append(v)
 obj=dict(review_id=rec['review_id'],IID=rec['IID'],sex=int(rec['SEX']),ancestry=rec['ancestry_pred'],case_status=int(rec['psoriasis_case']),center=str(rec['site_id']),targets=result,source_receipts=receipts,read_groups=rgs)
 (folder/'commands.private.json').write_text(json.dumps(command,indent=2));done.write_text(json.dumps(obj,indent=2,allow_nan=False));return obj
def run():
 assert (O/'PREPARED.json').exists() and not (O/'RUNNING.lock').exists()
 prepared=json.loads((O/'PREPARED.json').read_text());assert sha(O/'selection.private.tsv')==prepared['selection_sha256'];assert sha(O/'PROTOCOL.json')==prepared['protocol_sha256']
 lock=O/'RUNNING.lock';lock.write_text(str(os.getpid()));write('STATUS.json',dict(status='running',scope='bounded stratified three-target CRAM review'))
 selected=pd.read_csv(O/'selection.private.tsv',sep='\t',dtype={'IID':str});window=reference_window();results=[];failures=[]
 try:
  with ThreadPoolExecutor(max_workers=2) as pool:
   jobs={pool.submit(one,r,window):r['review_id'] for r in selected.to_dict('records')}
   for future in as_completed(jobs):
    try:results.append(future.result())
    except Exception as e:failures.append(dict(review_id=jobs[future],error=repr(e),traceback=traceback.format_exc()))
    write('progress.private.json',dict(completed=len(results),failed=len(failures),planned=len(selected)));write('failures.private.json',failures)
  write('results.private.json',results);write('STATUS.json',dict(status='read_stage_complete' if not failures else 'read_stage_failed',all_selected_retrieved=not failures,all_depth_crosschecks_passed=all(t['depth_crosscheck_passed'] for r in results for t in r['targets']),full_cohort_mapping_certified=False,requires_aggregate_review=True));print('READ_STAGE_ENDED; see STATUS.json',flush=True)
 finally:lock.unlink(missing_ok=True)
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='run':run()

