from pathlib import Path
import os,json,hashlib,subprocess,traceback,re,datetime
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import numpy as np,pandas as pd
from scipy.stats import fisher_exact
from google.cloud import storage
import google.auth
from google.auth.transport.requests import Request
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');O=N/'reconciliation_20260912/read_mapping_QC';F=O/'CRAM_pilot';F.mkdir(exist_ok=True);REF=O/'broad_reference/chrX.fa';POS=5998724
EXCLUDE=4|256|512|1024|2048

def at_base(cigar,start,seq,qual):
 r=start;q=0
 for size,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar):
  n=int(size)
  if op in 'M=X':
   if r<=POS<r+n:
    qi=q+POS-r;return (seq[qi].upper(),ord(qual[qi])-33 if qual!='*' else np.nan,qi)
   r+=n;q+=n
  elif op in 'DN':
   if r<=POS<r+n:return None
   r+=n
  elif op in 'IS':q+=n
 return None

def one(rec):
 ident=rec['pilot_id'];folder=F/ident;folder.mkdir(exist_ok=True);done=folder/'COMPLETE.private.json'
 if done.exists():return json.loads(done.read_text())
 creds,project=google.auth.default(scopes=['https://www.googleapis.com/auth/devstorage.read_only']);creds.refresh(Request());env=dict(os.environ,GCS_OAUTH_TOKEN=creds.token,GCS_REQUESTER_PAYS_PROJECT=os.environ['GOOGLE_PROJECT'])
 client=storage.Client();receipts={}
 for key in ['cram_uri','cram_index_uri']:
  bn,name=rec[key][5:].split('/',1);b=client.bucket(bn,user_project=os.environ['GOOGLE_PROJECT']).get_blob(name);assert b is not None;receipts[key]=dict(uri=rec[key],generation=b.generation,size=b.size,crc32c=b.crc32c)
  if key=='cram_index_uri':b.download_to_filename(str(folder/'input.cram.crai'),if_generation_match=b.generation,checksum='crc32c',timeout=60)
 commands=[['samtools','view','-h','-T',str(REF),'-X',rec['cram_uri'],str(folder/'input.cram.crai'),'chrX:5998424-5999024']]
 with (folder/'region.private.sam').open('w') as out,(folder/'extract.private.log').open('w') as err:p=subprocess.run(commands[0],stdout=out,stderr=err,env=env,timeout=90)
 if p.returncode:raise RuntimeError(ident+' targeted CRAM retrieval failed; inspect private log')
 rows=[];rgs=[];xmd5=None;record_count=0
 for line in (folder/'region.private.sam').read_text().splitlines():
  t=line.split(chr(9))
  if line.startswith('@RG'):rgs.append(dict(v.split(':',1) for v in t[1:]));continue
  if line.startswith('@SQ'):
   d=dict(v.split(':',1) for v in t[1:])
   if d['SN']=='chrX':xmd5=d.get('M5')
  if line.startswith('@'):continue
  record_count+=1;flag=int(t[1]);mq=int(t[4]);tags={x.split(':',2)[0]:x.split(':',2)[2] for x in t[11:] if len(x.split(':',2))==3}
  if flag&EXCLUDE or t[2]!='chrX' or t[9]=='*':continue
  v=at_base(t[5],int(t[3]),t[9],t[10])
  if v is None:continue
  base,bq,qi=v;rows.append(dict(template=tags.get('RG','')+'|'+t[0],base=base,BQ=bq,MQ=mq,reverse=bool(flag&16),read_length=len(t[9]),read_cycle=(len(t[9])-qi if flag&16 else qi+1),soft_clipped=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',t[5]) if op=='S'),mate_chrY=t[6]=='chrY',alternative_chrY=('chrY,' in tags.get('SA','') or 'chrY,' in tags.get('XA',''))))
 assert xmd5=='2b3a55ff7f58eb308420c8a9b11cac50','CRAM reference mismatch'
 bam=folder/'region.private.bam'
 commands += [['samtools','view','-b','-o',str(bam),str(folder/'region.private.sam')],['samtools','index',str(bam)],['samtools','depth','-r','chrX:5998724-5998724','-q','0','-Q','0','-G',str(EXCLUDE),str(bam)]]
 for k,cmd in enumerate(commands[1:]):
  p=subprocess.run(cmd,capture_output=True,text=True,timeout=30);(folder/f'check_{k}.private.log').write_text(p.stderr);assert p.returncode==0,(ident,'independent depth check',k)
 depth=0 if not p.stdout.strip() else int(p.stdout.strip().split(chr(9))[-1]);assert depth==len(rows),(ident,'CIGAR/depth mismatch',depth,len(rows))
 r=pd.DataFrame(rows);r.to_csv(folder/'read_observations.private.tsv',sep=chr(9),index=False)
 raw=Counter(x['base'] for x in rows);high=[x for x in rows if x['MQ']>=30 and x['MQ']<255 and x['BQ']>=20 and x['base'] in 'ACGT'];templates={};conflicts=0
 for x in high:templates.setdefault(x['template'],[]).append(x)
 chosen=[]
 for key,vals in templates.items():
  top=max(v['BQ'] for v in vals);best=[v for v in vals if v['BQ']==top]
  if len({v['base'] for v in best})>1:conflicts+=1;continue
  chosen.append(sorted(best,key=lambda v:(-v['MQ'],v['read_cycle']))[0])
 counts=Counter(x['base'] for x in chosen);mq=[x['MQ'] for x in rows if x['MQ']<255];g=[x for x in chosen if x['base']=='G'];a=[x for x in chosen if x['base']=='A'];strand=[[sum(not x['reverse'] for x in g),sum(x['reverse'] for x in g)],[sum(not x['reverse'] for x in a),sum(x['reverse'] for x in a)]]
 result=dict(pilot_id=ident,IID=rec['IID'],ancestry=rec['ancestry_pred'],sex=int(rec['SEX']),case_status=int(rec['psoriasis_case']),center=rec['site_id'],expected_NLGN4X_dosage=float(rec['NLGN4X_dosage']),carrier=bool(rec['carrier']),region_records=record_count,raw_primary_nonduplicate_read_depth=depth,raw_base_counts=dict(raw),MQ_median=float(np.median(mq)) if mq else None,MQ_below30_n=sum(x['MQ']<30 for x in rows),MQ_unknown255_n=sum(x['MQ']==255 for x in rows),HQ_template_depth=len(chosen),HQ_base_counts=dict(counts),HQ_target_fraction=counts['G']/len(chosen) if chosen else None,overlap_conflicts=conflicts,read_mate_chrY_n=sum(x['mate_chrY'] for x in rows),alternative_chrY_tag_n=sum(x['alternative_chrY'] for x in rows),strand_fisher_p=float(fisher_exact(strand).pvalue) if g and a else None,read_groups=rgs,source_receipts=receipts,depth_crosscheck_passed=True)
 (folder/'commands.private.json').write_text(json.dumps(commands,indent=2));done.write_text(json.dumps(result,indent=2));print('CRAM_PILOT_COMPLETE',ident,flush=True);return result
try:
 assert json.loads((O/'broad_reference/COMPLETE.json').read_text())['exact_reference_M5_matches_CRAM'];assert REF.exists();selection=pd.read_csv(O/'CRAM_pilot_selection.private.tsv',sep=chr(9),dtype={'IID':str});assert selection.IID.is_unique
 protocol=dict(created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),purpose='Stratified diagnostic pilot, not full-cohort mapping-QC certification or a re-estimated GWAS.',selection='Two available samples per EUR/AFR/AMR by female/male by case/control by target-carrier/noncarrier; deterministic pre-read selection.',region_GRCh38=[5998424,5999024],target_GRCh38=POS,reference='Broad Homo_sapiens_assembly38 with chrX M5 matched to every sampled CRAM',raw_depth='Aligned-base read depth excluding unmapped, secondary, supplementary, duplicate and QC-failed records; no MQ/BQ threshold; independent samtools depth verification.',HQ='MQ>=30 and <255, BQ>=20, A/C/G/T; overlapping templates counted once using highest BQ then MQ, conflicting equal-BQ allele observations removed.',genotype_quality='GQ taken from provider variant data separately; summed AD is not complete DP.',batch='Preserve CN/LB/PU and all RG fields privately; sampled header availability is not cohort-wide batch adjustment.',export='HELD')
 (F/'PROTOCOL.private.json').write_text(json.dumps(protocol,indent=2))
 with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,selection.to_dict('records')))
 (F/'results.private.json').write_text(json.dumps(results,indent=2));pd.DataFrame([{k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in results]).to_csv(F/'results.private.tsv',sep=chr(9),index=False)
 (F/'COMPLETE.private.json').write_text(json.dumps(dict(status='pilot_completed',participants=len(results),all_depth_crosschecks_passed=True,full_cohort_mapping_certified=False,release='HELD',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2));print('CRAM_PILOT_ALL_COMPLETE',flush=True)
except Exception as e:(F/'FAILED.private.json').write_text(json.dumps(dict(error=repr(e),traceback=traceback.format_exc()),indent=2));raise
