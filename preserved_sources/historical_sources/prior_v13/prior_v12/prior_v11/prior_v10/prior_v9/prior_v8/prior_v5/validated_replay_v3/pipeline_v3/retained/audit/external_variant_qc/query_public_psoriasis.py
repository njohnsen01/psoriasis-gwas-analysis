"""Public Pan-UKB variant QC lookup; no participant data is accessed."""
from pathlib import Path
import requests, gzip, struct, io, json, hashlib, datetime
ROOT=Path(__file__).resolve().parent
URL='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/phecode-696.4-both_sexes.tsv.bgz'
TARGETS=[('rs5970952',22935736,'C','T'),('rs183018930',5916765,'A','G'),('rs5919457',67062111,'G','T'),('rs5964607',66956431,'T','C')]
s=requests.Session()
h=s.head(URL,timeout=30);h.raise_for_status();etag=h.headers['ETag']
idx=s.get('https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/phecode-696.4-both_sexes.tsv.bgz.tbi',timeout=30);idx.raise_for_status();(ROOT/'public_psoriasis.tbi').write_bytes(idx.content)
f=io.BytesIO(gzip.decompress(idx.content));assert f.read(4)==b'TBI\x01'
def ints(n):return struct.unpack('<'+'i'*n,f.read(4*n))
nref,fmt,cseq,cbeg,cend,meta,skip,nl=ints(8);names=f.read(nl).decode().rstrip('\0').split('\0');refs={}
for name in names:
 bins={}
 for _ in range(ints(1)[0]):
  b,n=struct.unpack('<Ii',f.read(8));bins[b]=[struct.unpack('<QQ',f.read(16)) for _ in range(n)]
 ni=ints(1)[0];linear=struct.unpack('<'+'Q'*ni,f.read(8*ni));refs[name]=(bins,linear)
requests_log=[]
def fetch(a,b):
 assert b-a<30000000
 r=s.get(URL,headers={'Range':f'bytes={a}-{b}','If-Match':etag},timeout=60);r.raise_for_status();assert r.status_code==206
 requests_log.append({'start':a,'end':b,'sha256':hashlib.sha256(r.content).hexdigest()})
 return r.content
def unpack(raw):
 chunks=[];i=0
 while i+18<=len(raw):
  assert raw[i:i+3]==b'\x1f\x8b\x08'
  n=struct.unpack_from('<H',raw,i+16)[0]+1
  if i+n>len(raw):break
  chunks.append(gzip.decompress(raw[i:i+n]));i+=n
 return b''.join(chunks)
fields=unpack(fetch(0,65535)).decode().splitlines()[0].lstrip('#').split('\t');out=[]
for rs,pos,ref,alt in TARGETS:
 beg=max(0,pos-100001);end=pos+100000;rb=[0]
 for off,shift in [(1,26),(9,23),(73,20),(585,17),(4681,14)]:rb.extend(range(off+(beg>>shift),off+((end-1)>>shift)+1))
 bins,linear=refs['X'];low=linear[beg>>14] if beg>>14<len(linear) else 0
 chunks=sorted(set(c for k in rb for c in bins.get(k,[]) if c[1]>low));matches={}
 for a,b in chunks:
  text=unpack(fetch(a>>16,(b>>16)+65535))[(a&65535):].decode()
  for line in text.splitlines():
   v=line.split('\t')
   if len(v)==len(fields) and v[:4]==['X',str(pos),ref,alt]:matches[line]=dict(zip(fields,v))
 assert len(matches)==1,(rs,len(matches))
 rec=next(iter(matches.values()));rec['requested_rsid']=rs;out.append(rec)
result={'source':URL,'accessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'etag':etag,'last_modified':h.headers.get('Last-Modified'),'source_bytes':h.headers.get('Content-Length'),'index_sha256':hashlib.sha256(idx.content).hexdigest(),'query':'100-kb padded indexed intervals then exact chromosome-position-ref-alt matching','byte_ranges':requests_log,'variants':out}
(ROOT/'public_psoriasis_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'source':URL,'variants':out},indent=2))
