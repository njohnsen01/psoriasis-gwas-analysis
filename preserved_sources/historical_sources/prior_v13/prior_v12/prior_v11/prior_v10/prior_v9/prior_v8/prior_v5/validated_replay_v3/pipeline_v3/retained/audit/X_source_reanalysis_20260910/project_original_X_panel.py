from pathlib import Path
import sys,json,time,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'python_dependencies'))
import pgenlib
m=pd.read_csv(R/'original_X_normalized_projection_candidates.tsv',sep='\t').sort_values(['variant_index','allele_index'])
assert len(m)==376171 and m.ID.is_unique
samples=np.load(R/'source_sample_subset.npy'); offsets=np.load(R/'source_allele_offsets.npy')
ps=pd.read_csv(R/'source_order_frozen.psam',sep='\t',dtype={'FID':str,'IID':str})
assert len(ps)==len(samples)==392155
prefix=R/'private_corrected_X_original_panel'
ps.to_csv(str(prefix)+'.psam',sep='\t',index=False)
pv=m[['CHROM','POS','ID','REF','ALT']].copy();pv.columns=['#CHROM','POS','ID','REF','ALT'];pv.to_csv(str(prefix)+'.pvar',sep='\t',index=False)
buf=np.empty(len(samples),dtype=np.int8);started=time.time()
reader=pgenlib.PgenReader(str(R/'private_provider/acaf_threshold.chrX.pgen').encode(),sample_subset=samples,allele_idx_offsets=offsets)
writer=pgenlib.PgenWriter((str(prefix)+'.pgen').encode(),sample_ct=len(samples),variant_ct=len(m),nonref_flags=False)
try:
 for j,row in enumerate(m.itertuples(index=False)):
  reader.read(int(row.variant_index),buf,allele_idx=int(row.allele_index))
  writer.append_biallelic(buf)
  if j%10000==0:
   (R/'projection_status.json').write_text(json.dumps({'status':'running','completed':j,'total':len(m),'elapsed_seconds':time.time()-started}))
   print('PROJECTED',j,flush=True)
finally:
 writer.close();reader.close()
check=pgenlib.PgenReader((str(prefix)+'.pgen').encode());assert check.get_variant_ct()==len(m) and check.get_raw_sample_ct()==len(samples);check.close()
(R/'projection_status.json').write_text(json.dumps({'status':'complete','completed':len(m),'total':len(m),'elapsed_seconds':time.time()-started,'scope':'All original tested X variants; unaffected provider FT-filtered PGEN; unique minimal-allele representation mapping. Does not establish completeness of a newly QC-derived source-wide scan.'},indent=2))
print('PROJECTION_COMPLETE',flush=True)
