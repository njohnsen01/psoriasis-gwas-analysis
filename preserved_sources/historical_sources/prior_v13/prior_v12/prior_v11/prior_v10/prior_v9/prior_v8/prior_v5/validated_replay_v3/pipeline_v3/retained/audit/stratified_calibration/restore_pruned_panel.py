from pathlib import Path
import json,os,subprocess,hashlib,shutil,traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent
OLD=Path(f'/home/jupyter/psoriasis_regenie_loco_indexed_range_step1/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}')
PLAN=json.loads((OLD/'step1_assembly_plan.json').read_text())
KEEP=set((OLD/'global_pruning/psoriasis_step1_common_global_ld.prune.in').read_text().split())
PANEL=ROOT/'private_panel';PANEL.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def atomic(p,d):
 q=p.with_suffix(p.suffix+'.partial');q.write_text(json.dumps(d,indent=2));q.replace(p)
def run(cmd,log):
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(str(cmd[0])+' failed: '+log.read_text()[-3000:])
def restore(rec):
 c=rec['chromosome'];out=PANEL/('chr'+c);receipt=Path(str(out)+'.json')
 if receipt.exists():return json.loads(receipt.read_text())
 tmp=PANEL/('download_chr'+c);tmp.mkdir(exist_ok=True)
 try:
  if shutil.disk_usage(ROOT).free<150*1024**3:raise RuntimeError('Protected free disk threshold reached')
  for ext in ['bim','fam','bed']:
   p=tmp/rec[ext]['name'];item=rec[ext]
   run(['gsutil','-u',os.environ['GOOGLE_PROJECT'],'cp',item['versioned_uri'],str(p)],tmp/(ext+'.download.log'))
   if p.stat().st_size!=item['size_bytes'] or sha(p)!=item['sha256']:raise RuntimeError('Source checksum failed: chr'+c+' '+ext)
  fam=tmp/rec['fam']['name'];n=sum(1 for _ in fam.open());assert n==392155
  assert sha(fam)=='d5ddf8b5b901a7e6aa1cf093088ad1c960017ba930af2776c1c34976bfd631c9'
  lines=(tmp/rec['bim']['name']).read_text().splitlines(True);stride=(n+3)//4
  selected=[i for i,l in enumerate(lines) if l.split()[1] in KEEP]
  src=tmp/rec['bed']['name'];assert src.stat().st_size==3+stride*len(lines)
  with src.open('rb') as f,Path(str(out)+'.bed').open('wb') as g:
   assert f.read(3)==bytes([108,27,1]);g.write(bytes([108,27,1]))
   for i in selected:
    f.seek(3+i*stride);b=f.read(stride);assert len(b)==stride;g.write(b)
  Path(str(out)+'.bim').write_text(''.join(lines[i] for i in selected));shutil.copyfile(fam,Path(str(out)+'.fam'))
  d={'chromosome':c,'selected_markers':len(selected),'source_checksums_verified':True,'bed_sha256':sha(Path(str(out)+'.bed'))}
  atomic(receipt,d);print('RESTORED',c,len(selected),flush=True);return d
 finally:
  shutil.rmtree(tmp,ignore_errors=True)
try:
 atomic(ROOT/'panel_status.json',{'status':'restoring','started_utc':datetime.now(timezone.utc).isoformat()})
 with ThreadPoolExecutor(max_workers=2) as ex: receipts=list(ex.map(restore,PLAN['artifact_plan']))
 prefix=PANEL/'calibration_pruned';stride=(392155+3)//4
 with Path(str(prefix)+'.bed').open('wb') as out,Path(str(prefix)+'.bim').open('w') as bim:
  out.write(bytes([108,27,1]))
  for c in range(1,23):
   p=PANEL/('chr'+str(c))
   with Path(str(p)+'.bed').open('rb') as f:
    assert f.read(3)==bytes([108,27,1]);shutil.copyfileobj(f,out,8*1024**2)
   bim.write(Path(str(p)+'.bim').read_text())
 shutil.copyfile(PANEL/'chr1.fam',Path(str(prefix)+'.fam'))
 ids=[l.split()[1] for l in Path(str(prefix)+'.bim').read_text().splitlines()]
 assert len(ids)==247255 and len(set(ids))==247255 and set(ids)==KEEP
 assert Path(str(prefix)+'.bed').stat().st_size==3+stride*len(ids)
 atomic(ROOT/'panel_status.json',{'status':'validated','completed_utc':datetime.now(timezone.utc).isoformat(),'n_samples':392155,'n_markers':len(ids),'prefix':str(prefix),'bed_sha256':sha(Path(str(prefix)+'.bed')),'chromosomes':receipts})
 print('PANEL_VALIDATED',len(ids),flush=True)
except Exception as e:
 atomic(ROOT/'panel_status.json',{'status':'failed','error':str(e),'traceback':traceback.format_exc()});raise
