from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--qki-raw',type=Path,required=True);a=p.parse_args()
x=pd.read_pickle(a.inputs/'X_inputs.private.pkl');b=pd.read_pickle(a.inputs/'biomarker_inputs.private.pkl')
assert x.index.is_unique and b.index.is_unique and set(x.index)==set(b.index)
r=pd.read_csv(a.qki_raw,sep='\t',dtype={'IID':str}).set_index('IID');assert r.index.is_unique and b.index.isin(r.index).all()
col='chr6:163803171:CTTATT:C_CTTATT';assert col in r
b['ancestry_pred']=x.ancestry_pred.reindex(b.index);b['QKI']=2-pd.to_numeric(r[col],errors='coerce').reindex(b.index)
assert b.QKI.dropna().between(0,2).all()
out=a.inputs/'ancestry_autosomal_inputs.private.pkl';assert not out.exists();b.to_pickle(out)
(a.inputs/'ancestry_input_provenance.json').write_text(json.dumps({'QKI':'chr6:163803171:CTTATT:C','effect_allele':'C','source_counted_allele':'CTTATT','transformation':'2-minus-raw','raw_sha256':hashlib.sha256(a.qki_raw.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source':'reconstructed matrices plus original raw PLINK export'},indent=2))
print('ANCESTRY_MATRIX_REBUILT',len(b))
