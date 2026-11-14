from pathlib import Path
import ast,csv,hashlib,json,datetime,zipfile,xml.etree.ElementTree as E
import numpy as np,pandas as pd,openpyxl
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');V=I/'supplement_staging_v5';D=V/'methods_sources';src=I.parent.parent/'SII_timing_sensitivity_20260910/run_sii.py'
rows=list(csv.DictReader((D/'SII_unit_mapping.private.tsv').open(),delimiter='\t'));pairs=[r for r in rows if r['record_type']=='observed_pair']
c=pd.DataFrame(pairs);env={'np':np,'pd':pd,'c':c};selected=[]
for n in ast.walk(ast.parse(src.read_text())):
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ['standard','source','a','b','fac'] for t in n.targets): selected.append(n)
assert len(selected)==5
selected.sort(key=lambda n:(n.lineno,n.col_offset));exec(compile(ast.Module(body=selected,type_ignores=[]),str(src),'exec'),env)
want=pd.to_numeric(c['factor_to_thousand_cells_per_uL'],errors='coerce').to_numpy();got=env['fac'].to_numpy();assert np.allclose(want,got,equal_nan=True)
w=openpyxl.load_workbook(V/'S1_reconciled_tables_HELD.xlsx',read_only=True,data_only=False);sh=w['SII_unit_mapping'];head=list(next(sh.values));assert head==list(rows[0]);vals=list(sh.values)[1:];assert len(vals)==len(rows)
for r,v in zip(rows,vals):
 for k,x in zip(head,v):
  z=r[k]
  if k=='factor_to_thousand_cells_per_uL': assert (x is None and z=='') or (isinstance(x,(int,float)) and float(z)==x)
  else: assert (x or '')==z
w.close();ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(I/'supplement_staging_v4/S1_reconciled_tables_HELD.xlsx') as old,zipfile.ZipFile(V/'S1_reconciled_tables_HELD.xlsx') as new:
 assert new.testzip() is None
 for j in range(1,35):
  name='xl/worksheets/sheet%d.xml'%j;a=old.read(name);b=new.read(name)
  if j not in [1,34]: assert a==b,name
  else:
   a=E.fromstring(a);b=E.fromstring(b);nb={c.attrib['r']:E.tostring(c) for c in b.findall('.//s:sheetData/s:row/s:c',ns)}
   assert all(nb[c.attrib['r']]==E.tostring(c) for c in a.findall('.//s:sheetData/s:row/s:c',ns)),name
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
reused=[]
for pattern in ['Lab_methods_*.png','Phenotype_concepts_*.png']:
 for p in (V/'review_proofs_final').glob(pattern): assert sha(p)==sha(V/'review_proofs'/p.name);reused.append(p.name)
receipt={'verified_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'unit_rule_source':str(src),'source_sha256':sha(src),'actual_AST_expressions_evaluated':5,'observed_unit_pairs_verified':len(pairs),'typed_unit_table_verified':True,'all_original_v4_cells_preserved':True,'original_32_unappended_sheets_byte_identical':True,'reviewed_initial_images_identical_to_final':sorted(reused),'scope':'No association or phenotype model rerun; new Methods source/table verification only.'}
(V/'INDEPENDENT_METHODS_AUDIT.private.json').write_text(json.dumps(receipt,indent=2));print('V5_INDEPENDENT_METHODS_AUDIT_PASSED',len(pairs),'unit pairs;',len(reused),'proofs identical')
