from pathlib import Path
import os,json,hashlib,subprocess,zstandard
from google.cloud import storage
N=Path(__file__).resolve().parents[1];out=N/'scratch/pilot_chr22';out.mkdir(exist_ok=True)
c=storage.Client();rec={}
for ext,uri in json.loads((N/'pilot_sources.json').read_text()).items():
 bucket,name=uri[5:].split('/',1);b=c.bucket(bucket,user_project=os.environ['GOOGLE_PROJECT']).get_blob(name);assert b is not None
 path=out/('input.'+('pvar.zst' if ext=='pvar' else ext))
 if not path.exists():b.download_to_filename(str(path),if_generation_match=b.generation,checksum='crc32c')
 assert path.stat().st_size==b.size
 h=hashlib.sha256(path.read_bytes()).hexdigest();rec[ext]={'uri':uri,'generation':b.generation,'md5':b.md5_hash,'bytes':b.size,'sha256':h}
 if ext=='pvar':
  with path.open('rb') as f,(out/'input.pvar').open('wb') as g:zstandard.ZstdDecompressor().copy_stream(f,g)
(N/'status/pilot_source_receipt.json').write_text(json.dumps(rec,indent=2))
for group in ['cases','controls','EUR.controls','AFR.controls','AMR.controls']:
 cmd=['plink2','--pfile',str(out/'input'),'--keep',str(N/'inputs'/(group+'.keep.tsv')),'--threads','4','--memory','4096','--out',str(out/group)]
 cmd+=['--freq','counts','--missing','variant-only'] if group in ['cases','controls'] else ['--hardy']
 (out/(group+'.command.json')).write_text(json.dumps(cmd))
 with (out/(group+'.console.log')).open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,(group,p.returncode)
 print('PILOT_QC',group,flush=True)
print('PILOT_QC_COMPLETE',flush=True)
