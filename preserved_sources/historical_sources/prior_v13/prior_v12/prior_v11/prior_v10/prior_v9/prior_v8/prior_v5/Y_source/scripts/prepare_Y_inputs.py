from pathlib import Path
import argparse,base64,hashlib,json,os,shutil,sys
import pandas as pd
from google.cloud import storage
import google_crc32c
Y=Path(__file__).resolve().parents[1];N=Y.parent
p=argparse.ArgumentParser();p.add_argument('--download-missing',action='store_true');args=p.parse_args()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def crc(p):
 c=google_crc32c.Checksum()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):c.update(b)
 return base64.b64encode(c.digest()).decode()
manifest=json.loads((Y/'provider_manifest.json').read_text());receipts=[]
client=storage.Client(project=os.environ['GOOGLE_PROJECT']);bucket=client.bucket('vwb-aou-datasets-controlled',user_project=os.environ['GOOGLE_PROJECT'])
(Y/'source').mkdir(exist_ok=True)
for r in manifest:
 dest=Y/'source'/Path(r['name']).name
 if not dest.exists():
  assert args.download_missing,'Pinned source missing; rerun with --download-missing within authorized Workbench'
  assert shutil.disk_usage(Y).free>128*1024**3+int(r['bytes'])
  blob=bucket.blob(r['name'],generation=int(r['generation']));blob.reload()
  assert blob.size==r['bytes'] and blob.crc32c==r['crc32c']
  part=dest.with_suffix(dest.suffix+'.partial');blob.download_to_filename(str(part),if_generation_match=int(r['generation']),checksum='crc32c');part.rename(dest)
 assert dest.stat().st_size==r['bytes'] and crc(dest)==r['crc32c']
 receipts.append({'name':r['name'],'generation':r['generation'],'sha256':sha(dest),'bytes':dest.stat().st_size})
ps=pd.read_csv(Y/'source/acaf_threshold.chrY.psam',sep=chr(9),dtype={'#IID':str})
pa=pd.read_csv(N/'inputs/provider_X_aligned.psam',sep=chr(9),dtype={'#FID':str,'IID':str})
ph=pd.read_csv(N/'stratified_step1/male/phenotype.tsv',sep=chr(9),dtype={'FID':str,'IID':str})
assert ps['#IID'].equals(pa.IID) and ps['#IID'].is_unique and ph.IID.is_unique
assert len(ph)==151282 and ph.psoriasis_case.eq(1).sum()==2259
assert pa.loc[pa.IID.isin(ph.IID),'SEX'].eq(1).all() and pa.IID.isin(ph.IID).sum()==len(ph)
aligned=Y/'aligned_provider.psam'
if aligned.exists():
 old=pd.read_csv(aligned,sep=chr(9),dtype={'#FID':str,'IID':str});pd.testing.assert_frame_equal(old,pa)
else:pa.to_csv(aligned,sep=chr(9),index=False)
record={'status':'validated','source_receipts':receipts,'provider_X_aligned_psam_sha256':sha(N/'inputs/provider_X_aligned.psam'),'aligned_Y_psam_sha256':sha(aligned),'male_phenotype_sha256':sha(N/'stratified_step1/male/phenotype.tsv'),'script_sha256':sha(Path(__file__)),'scope':'Pinned source and cohort alignment; phenotype construction and step-1 prerequisite pipeline remains in parent analysis'}
(Y/'PREPARATION_VERIFIED.json').write_text(json.dumps(record,indent=2)+chr(10))
print('Y_PREPARATION_VERIFIED')
