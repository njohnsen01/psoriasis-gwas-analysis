from pathlib import Path
import json,requests,hashlib,tarfile,subprocess,os
O=Path(__file__).resolve().parent;T=O/'tools';j=json.loads((T/'bwa_release.json').read_text());a=next(x for x in j['assets'] if x['name']=='bwa-mem2-2.3_x64-linux.tar.bz2');p=T/a['name']
r=requests.get(a['browser_download_url'],timeout=180);r.raise_for_status();p.write_bytes(r.content);sha=hashlib.sha256(r.content).hexdigest()
if a.get('digest'):assert a['digest']=='sha256:'+sha
with tarfile.open(p) as f:
 for m in f.getmembers():assert not m.issym() and not m.islnk() and (T/m.name).resolve().is_relative_to(T.resolve())
 f.extractall(T)
b=next(T.rglob('bwa-mem2'));v=subprocess.run([str(b),'version'],capture_output=True,text=True,check=True)
(O/'BWA_SOURCE.json').write_text(json.dumps(dict(version=v.stdout.strip(),binary=str(b),source=a['browser_download_url'],sha256=sha,digest_verified=bool(a.get('digest'))),indent=2))
I=O/'bwa_index';I.mkdir();prefix=I/'GRCh38';ref=O.parent/'mapping_extension_20260919/reference/Homo_sapiens_assembly38.fasta'
with (O/'bwa_index.log').open('w') as log:subprocess.run([str(b),'index','-p',str(prefix),str(ref)],stdout=log,stderr=log,env=dict(os.environ,OMP_NUM_THREADS='2'),check=True,timeout=5400)
(O/'BWA_INDEX_READY').write_text(str(prefix));print('BWA_INDEX_READY',flush=True)
