from pathlib import Path
import hashlib,json,xml.etree.ElementTree as E,zipfile
import openpyxl
I=Path(__file__).resolve().parents[1];V=I/'supplement_staging_v5';builder=I/'scripts/build_methods_integration_v5.py'
s=builder.read_text();old="for s,rows in expected.items():assert list(z[s].values)==rows,s";new="for s,rows in expected.items():assert [tuple(None if v=='' else v for v in r) for r in z[s].values]==[tuple(None if v=='' else v for v in r) for r in rows],s"
assert old in s;s=s.replace(old,new);compile(s,str(builder),'exec');builder.write_text(s)
ns={'__name__':'__main__','__file__':str(builder)}
prefix=s[:s.index('oldpath=')].replace(";assert not (V/'S1_reconciled_tables_HELD.xlsx').exists()",'')
exec(compile(prefix,str(builder),'exec'),ns)
oldpath=I/'supplement_staging_v4/S1_reconciled_tables_HELD.xlsx';out=V/'S1_reconciled_tables_HELD.xlsx';oldhash=ns['sha'](oldpath)
a=zipfile.ZipFile(oldpath);b=zipfile.ZipFile(out);xn={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
files=[n for n in a.namelist() if n.startswith('xl/worksheets/sheet') and n.endswith('.xml')];diff=[]
for n in files:
 if a.read(n)==b.read(n):continue
 diff.append(n);assert n in ['xl/worksheets/sheet1.xml','xl/worksheets/sheet34.xml']
 aa=E.fromstring(a.read(n));bb=E.fromstring(b.read(n));ac={c.attrib['r']:E.tostring(c) for c in aa.findall('.//x:sheetData/x:row/x:c',xn)};bc={c.attrib['r']:E.tostring(c) for c in bb.findall('.//x:sheetData/x:row/x:c',xn)}
 assert all(bc.get(k)==v for k,v in ac.items())
# Existing style definitions also remain unchanged.
aa=E.fromstring(a.read('xl/styles.xml'));bb=E.fromstring(b.read('xl/styles.xml'))
for kind in ['fonts','fills','borders','cellStyleXfs','cellXfs','cellStyles','numFmts']:
 x=aa.find('x:'+kind,xn);y=bb.find('x:'+kind,xn)
 if x is not None:assert y is not None and all(E.tostring(c)==E.tostring(list(y)[j]) for j,c in enumerate(x)),kind
z=openpyxl.load_workbook(out,read_only=True,data_only=False)
for name,df in ns['tables'].items():
 expected=[tuple(df.columns)]+[tuple(None if (v=='' or ns['pd'].isna(v)) else (v.item() if hasattr(v,'item') else v) for v in row) for row in df.itertuples(index=False,name=None)]
 assert list(z[name].values)==expected,name
ns.update({'z':z,'out':out,'oldpath':oldpath,'oldhash':oldhash})
exec(compile(s[s.index('cohort_text='):],str(builder),'exec'),ns)
p=V/'BUILD.private.json';r=json.loads(p.read_text());r['finalizer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();r['checks']['all_prior_worksheet_XML_unchanged_except_documented_appends']=True;r['checks']['old_style_definitions_preserved']=True;r['finalization_note']='A strict comparison stopped at empty-string versus XLSX blank serialization. Finalization normalized only empty string to blank, rechecked new typed values against source tables, confirmed every prior worksheet XML except documented README/provenance appends byte-for-byte, preserved each old cell XML in the appended sheets, and preserved all prior style definitions. No data/model rerun.';p.write_text(json.dumps(r,indent=2)+chr(10));print('V5_FINALIZATION_VERIFIED')
