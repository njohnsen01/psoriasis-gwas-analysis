from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
from openpyxl import load_workbook
import pandas as pd
import numpy as np
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');old=I/'supplement_staging_v3';D=I/'supplement_staging_v4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert not (D/'BUILD.private.json').exists()
parent=old/'S1_reconciled_tables_HELD.xlsx';dest=D/parent.name
assert sha(parent)=='d4061158b792b74f8d7a59ace57d2aa0ecf91c23839a6cc71fde22ac0abd9e5f'
audit_path=I/'supplement_review_v4/INDEPENDENT_AUDIT.private.json';audit=json.loads(audit_path.read_text());assert audit['status']=='PASS'
for p,h in audit['source_hashes'].items():assert sha(p)==h
w=load_workbook(dest,read_only=True);a=load_workbook(parent,read_only=True)
assert [n for n in w.sheetnames if n!='Sex_heterogeneity_current']==a.sheetnames
for name in a.sheetnames:
    for oldrow,newrow in zip(a[name].iter_rows(values_only=True),w[name].iter_rows(values_only=True)):
        assert oldrow==newrow[:len(oldrow)]
sex=pd.read_csv(I/'sex_heterogeneity_current/sex_heterogeneity.private.tsv',sep='	')
rows=list(w['Sex_heterogeneity_current'].values);assert list(rows[0])==list(sex.columns) and len(rows)==len(sex)+1
for actual,expected in zip(rows[1:],sex.itertuples(index=False,name=None)):
    for v,x in zip(actual,expected):
        if pd.isna(x):assert v is None
        elif isinstance(x,(float,int)) and not isinstance(x,bool):np.testing.assert_allclose(v,x,rtol=2e-14,atol=1e-300)
        else:assert v==x
manifest_rows=list(w['Source_manifest'].values);new_sources=[dict(zip(manifest_rows[0],r)) for r in manifest_rows[-3:]]
for r in new_sources:assert sha(r['source'])==r['sha256']
s2path=D/'S2_reconciliation_integration.private.md';s2=s2path.read_text();parent_s2=(old/s2path.name).read_text()
inserted=s2.split('## Independent review of center-refit results and dosage units')[1].split('## Reproducibility and release')[0]
assert s2.replace('## Independent review of center-refit results and dosage units'+inserted,'',1)==parent_s2
ms=json.loads((old/'S2_source_manifest.private.json').read_text());ms.update(parent_manifest={'path':str(old/'S2_source_manifest.private.json'),'sha256':sha(old/'S2_source_manifest.private.json')},v4_added_sources=new_sources+[{'source':str(I/'sex_heterogeneity_current/PROTOCOL.private.json'),'sha256':sha(I/'sex_heterogeneity_current/PROTOCOL.private.json')}],created_UTC=datetime.now(timezone.utc).isoformat(),script_sha256=sha(I/'scripts/build_supplement_v4.py'),release='HELD',not_final_S2=True,S2_sha256=sha(s2path))
(D/'S2_source_manifest.private.json').write_text(json.dumps(ms,indent=2)+chr(10))
receipt={'status':'SOURCE_AND_PRESERVATION_CHECKS_PASS','created_UTC':datetime.now(timezone.utc).isoformat(),'release':'HELD','not_final_S1':True,'not_final_S2':True,'parent_workbook':str(parent),'parent_workbook_sha256':sha(parent),'S1_sha256':sha(dest),'S2_sha256':sha(s2path),'script_sha256':sha(I/'scripts/build_supplement_v4.py'),'finalization_script_sha256':sha(__file__),'independent_audit_sha256':sha(audit_path),'new_sources':new_sources,'original_sheet_count':57,'draft_sheet_count':len(w.sheetnames),'all_original_cells_and_styles_preserved':True,'all_new_source_cells_verified':True,'native_Excel_Word_review':False,'new_scientific_files_exported':False,'execution_note':'Builder reached all workbook source/style and preservation checks, then python -c invocation lacked __file__ during manifest writing. This finalizer rechecked original cells, new sex cells, source hashes and S2 insertion preservation and completed receipts. Invoke builder as a Python file or with runpy, not bare exec.'}
(D/'BUILD.private.json').write_text(json.dumps(receipt,indent=2)+chr(10));json.loads((D/'BUILD.private.json').read_text());w.close();a.close()
print('V4_SOURCE_AND_PRESERVATION_CHECKS_PASS; complete valid receipts saved; export HELD')
