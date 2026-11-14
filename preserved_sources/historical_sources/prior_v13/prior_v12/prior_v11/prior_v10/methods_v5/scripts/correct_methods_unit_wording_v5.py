from pathlib import Path
import datetime,hashlib,json,math,xml.etree.ElementTree as E,zipfile
import pandas as pd
I=Path(__file__).resolve().parents[1];V=I/'supplement_staging_v5';D=V/'methods_sources';ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};uri=ns['x']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
oldkey='factor_to_10e3_cells_per_uL';newkey='factor_to_thousand_cells_per_uL'
oldrule='Use known factor; if both labels map, require agreement; reject unknown/conflicting pair.'
newrule='Require a recognized standard-unit name. If the source-unit label also maps, its factor must agree. An unrecognized source label is allowed when the standard name is recognized.'
oldtext='Under strict mapping, a known factor from the standard unit name or source unit label could be used; when both mapped, the factors had to agree. Unknown or conflicting pairs were excluded.'
newtext='Under strict mapping, the standard unit name had to match the explicit dictionary. If the source unit label also matched its dictionary, its factor had to agree with the standard-unit factor. An unrecognized source label was allowed when the standard unit name was recognized; an unknown standard unit name or conflicting known factors led to exclusion.'
p=I/'scripts/build_methods_integration_v5.py';s=p.read_text();assert oldkey in s and oldrule in s and oldtext in s;s=s.replace(oldkey,newkey).replace(oldrule,newrule).replace(oldtext,newtext)
s=s.replace("'factor_to_thousand_cells_per_uL':r['factor']","'factor_to_thousand_cells_per_uL':float(r['factor']) if r['factor']!='' else None")
compile(s,str(p),'exec');p.write_text(s)
f=D/'SII_unit_mapping.private.tsv';df=pd.read_csv(f,sep=chr(9),keep_default_na=False);df=df.rename(columns={oldkey:newkey});df.loc[df.record_type.isin(['standard_dictionary','source_dictionary']),'handling']=newrule;df[newkey]=pd.to_numeric(df[newkey],errors='coerce');df.to_csv(f,sep=chr(9),index=False)
xlsx=V/'S1_reconciled_tables_HELD.xlsx';before=sha(xlsx)
with zipfile.ZipFile(xlsx) as z:entries={x.filename:(x,z.read(x.filename)) for x in z.infolist()}
unitname='xl/worksheets/sheet37.xml';manifestname='xl/worksheets/sheet34.xml';tree=E.fromstring(entries[unitname][1]);cells={c.attrib['r']:c for c in tree.findall('.//x:sheetData/x:row/x:c',ns)}
def setcell(c,value):
 for child in list(c):c.remove(child)
 if value is None or (isinstance(value,float) and math.isnan(value)):c.attrib.pop('t',None);return
 if isinstance(value,(float,int)):
  c.attrib['t']='n';E.SubElement(c,'{'+uri+'}v').text=repr(value)
 else:
  c.attrib['t']='inlineStr';E.SubElement(E.SubElement(c,'{'+uri+'}is'),'{'+uri+'}t').text=str(value)
setcell(cells['D1'],newkey)
for j,r in enumerate(df.to_dict('records'),2):
 setcell(cells['D'+str(j)],None if pd.isna(r[newkey]) else float(r[newkey]));setcell(cells['E'+str(j)],r['handling'])
 row=tree.find('.//x:sheetData/x:row[@r="'+str(j)+'"]',ns);row.set('ht',str(max(float(row.get('ht','30')),61 if r['record_type']!='observed_pair' else 47)));row.set('customHeight','1')
mt=E.fromstring(entries[manifestname][1]);mc={c.attrib['r']:c for c in mt.findall('.//x:sheetData/x:row/x:c',ns)};assert ''.join(mc['A37'].itertext())=='SII_unit_mapping';setcell(mc['C37'],sha(f))
changes={unitname:E.tostring(tree),manifestname:E.tostring(mt)};tmp=xlsx.with_suffix('.tmp.xlsx')
with zipfile.ZipFile(tmp,'w') as z:
 for name,(info,data) in entries.items():z.writestr(info,changes.get(name,data))
with zipfile.ZipFile(tmp) as z:
 assert z.testzip() is None
 assert all(z.read(k)==v[1] for k,v in entries.items() if k not in changes)
 # Only the newly appended Source_manifest hash cell is changed.
 aa=E.fromstring(entries[manifestname][1]);bb=E.fromstring(z.read(manifestname));ac={c.attrib['r']:E.tostring(c) for c in aa.findall('.//x:sheetData/x:row/x:c',ns)};bc={c.attrib['r']:E.tostring(c) for c in bb.findall('.//x:sheetData/x:row/x:c',ns)};assert all(v==bc[k] for k,v in ac.items() if k!='C37')
tmp.replace(xlsx)
p=V/'S2_reconciliation_integration.private.md';s=p.read_text();assert oldtext in s;p.write_text(s.replace(oldtext,newtext))
p=I/'scripts/render_methods_v5_proofs.py';s=p.read_text();assert oldkey in s;s=s.replace(oldkey,newkey).replace("O=V/'review_proofs'","O=V/'review_proofs_final'");compile(s,str(p),'exec');p.write_text(s)
r=json.loads((V/'BUILD.private.json').read_text());r['builder_sha256']=sha(I/'scripts/build_methods_integration_v5.py');r['unit_wording_correction']={'verified_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_expression':'fac=a.where(a.notna()&(b.isna()|np.isclose(a,b,equal_nan=False)))','source_file':str(N if False else I.parent.parent/'SII_timing_sensitivity_20260910/run_sii.py'),'old_workbook_sha256':before,'new_workbook_sha256':sha(xlsx),'scope':'Only new unit-map header, handling text, numeric typing/row heights, its appended source hash and S2 explanatory sentence. Original association sheets remain byte-identical.'};r['files']={str(p.relative_to(V)):sha(p) for p in V.rglob('*') if p.is_file() and p.name!='BUILD.private.json' and 'review_proofs' not in str(p)};(V/'BUILD.private.json').write_text(json.dumps(r,indent=2)+chr(10));print('STRICT_UNIT_WORDING_AND_NUMERIC_TYPING_VERIFIED')
