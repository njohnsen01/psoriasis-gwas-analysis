from pathlib import Path
import json,hashlib,sys,os,traceback,importlib.util
import pandas as pd,numpy as np,zstandard
from threadpoolctl import threadpool_limits
_limit=threadpool_limits(limits=1)
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';R=Path('/home/jupyter/psoriasis_reproducibility_20260910');O=I/'regional_QKI';O.mkdir(exist_ok=True)
sys.path.insert(0,str(N.parent/'X_source_reanalysis_20260910/python_dependencies'));import pgenlib
spec=importlib.util.spec_from_file_location('primary_auto',N/'scripts/run_autosomes.py');auto=importlib.util.module_from_spec(spec);spec.loader.exec_module(auto)
sha=auto.sha
save=lambda name,value:(O/name).write_text(json.dumps(value,indent=2,default=str)+chr(10))
def main():
 ass=pd.read_csv(I/'QKI_primary_audit/full_QKI_regional_associations.private.tsv.gz',sep=chr(9));assert ass.ID.is_unique and ass.TEST.eq('ADD').all();assert ass.GENPOS.between(163303171,164303171).all()
 chunks=sorted({Path(s).parent.name for s in ass.source});assert chunks==['chr6_chr6_pso_chunk_0063'];folder=O/'source_chunk_0063';folder.mkdir(exist_ok=True)
 rec=json.loads((N/'results/autosomes'/chunks[0]/'source_receipts.json').read_text());receipt={}
 for ext in ['pgen','pvar','psam']:
  path=folder/('input.'+('pvar.zst' if ext=='pvar' else ext));receipt[ext]=auto.download(rec[ext]['uri'],path)
  for k in ['generation','crc32c','bytes','sha256']:assert str(receipt[ext][k])==str(rec[ext][k]),(ext,k)
  if ext=='pvar':
   with path.open('rb') as f,(folder/'input.pvar').open('wb') as g:zstandard.ZstdDecompressor().copy_stream(f,g)
 save('source_receipts.private.json',receipt)
 pv=pd.read_csv(folder/'input.pvar',sep=chr(9));ps=pd.read_csv(folder/'input.psam',sep=chr(9),dtype={'IID':str});assert pv.ID.is_unique and ps.IID.is_unique
 q=pd.read_csv(N/'results/autosomes'/chunks[0]/'variant_QC.private.tsv.gz',sep=chr(9)).set_index('ID');pos=pv.set_index('ID').POS
 region=q.loc[pos.between(163303171,164303171)];assert set(region.index[region.PASS])==set(ass.ID)
 p=pv.set_index('ID').loc[ass.ID];assert np.array_equal(p.REF,ass.ALLELE0) and np.array_equal(p.ALT,ass.ALLELE1)
 for src,part in ass.groupby('source'):assert sha(Path(src))==part.source_sha256.iloc[0]
 D=pd.read_pickle(R/'reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');D.index=D.index.astype(str)
 select=D.loc[D.unrelated & D.psoriasis_case.eq(0) & D.ancestry_pred.str.upper().isin(['EUR','AFR','AMR'])]
 ix=np.flatnonzero(ps.IID.isin(select.index)).astype(np.uint32);d=select.reindex(ps.IID.iloc[ix]);assert len(d)==len(select) and d.index.is_unique
 groups={(anc,sex):(d.ancestry_pred.str.upper().eq(anc)&(True if sex=='combined' else d.sex.eq(1 if sex=='male' else 0))).to_numpy() for anc in ['EUR','AFR','AMR'] for sex in ['combined','female','male']}
 lookup=pd.Series(np.arange(len(pv)),index=pv.ID);ref='chr6:163803171:CTTATT:C';buf=np.empty(len(ix),np.int8)
 reader=pgenlib.PgenReader(str(folder/'input.pgen').encode(),sample_subset=ix)
 def read(vid):
  reader.read(int(lookup.loc[vid]),buf);z=buf.astype(float);z[buf<0]=np.nan;return z
 base=read(ref);assert np.allclose(base,d.QKI.to_numpy(),equal_nan=True,rtol=0,atol=0),'QKI exact ALT dosage vector'
 rows=[]
 try:
  for k,vid in enumerate(ass.ID):
   z=read(vid)
   for (anc,sex),use in groups.items():
    x=z[use];y=base[use];ok=np.isfinite(x)&np.isfinite(y);n=int(ok.sum());acx=float(x[ok].sum());acy=float(y[ok].sum());macx=min(acx,2*n-acx);macy=min(acy,2*n-acy);cx=float(np.isfinite(x).mean());cy=float(np.isfinite(y).mean());reason=[]
    if n<500:reason.append('pair_complete_n')
    if min(macx,macy)<20:reason.append('MAC')
    if min(cx,cy)<.98:reason.append('call_rate')
    r2=np.nan
    if not reason:
     xx=x[ok]-x[ok].mean();yy=y[ok]-y[ok].mean();den=float(np.dot(xx,xx)*np.dot(yy,yy))
     if den<=0:reason.append('zero_variance')
     else:r2=float(np.dot(xx,yy)**2/den);assert -1e-10<=r2<=1+1e-10;r2=min(1.,max(0.,r2))
    rows.append(dict(region='QKI',variant=vid,reference=ref,ancestry=anc,sex=sex,r2=r2,status='available' if not reason else ','.join(reason),pair_complete_n=n,variant_MAC=macx,reference_MAC=macy,variant_callrate=cx,reference_callrate=cy))
   if k%500==0:print('QKI_LD_PROGRESS',k,len(ass),flush=True)
 finally:reader.close()
 ld=pd.DataFrame(rows);assert not ld.duplicated(['variant','ancestry','sex']).any();assert len(ld)==len(ass)*9
 for key,g in ld.groupby(['ancestry','sex']):assert set(g.variant)==set(ass.ID)
 ld.to_csv(O/'ancestry_control_LD.private.tsv.gz',sep=chr(9),index=False);region.to_csv(O/'regional_QC.private.tsv.gz',sep=chr(9));ass.to_csv(O/'full_associations.private.tsv.gz',sep=chr(9),index=False)
 save('PROTOCOL.private.json',dict(build='GRCh38',window=[163303171,164303171],reference=ref,LD='Pairwise-complete Pearson dosage correlation squared in provider-pruned controls, ancestry-specific combined-sex primary display and separate sex diagnostics. Unphased autosomal dosage LD; no D-prime or fine-mapping claim.',thresholds='pair-complete N>=500; MAC>=20 for each variant in paired observations; each call rate>=98%. Estimability thresholds, not power.',association='Every finite QC-passing frozen pooled test in the declared region; all rows shown without sampling.',release='HELD'))
 save('COMPLETE.private.json',dict(status='verified',association_records=len(ass),LD_records=len(ld),all_regional_QC_passing_tests_included=True,reference_vector_exactly_matches_followup=True,input_receipt=receipt,script_sha256=sha(Path(__file__)),release='HELD'));print('QKI_REGIONAL_LD_COMPLETE',flush=True)
try:main()
except Exception as e:save('FAILED.private.json',dict(error=repr(e),traceback=traceback.format_exc()));raise
