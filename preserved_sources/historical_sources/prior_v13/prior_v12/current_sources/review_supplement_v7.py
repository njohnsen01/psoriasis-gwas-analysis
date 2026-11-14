from pathlib import Path
import datetime,hashlib,json,math
import openpyxl
from PIL import Image,ImageDraw,ImageFont
from matplotlib import font_manager
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');V=I/'supplement_staging_v7';O=V/'review_proofs'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
b=json.loads((V/'BUILD.private.json').read_text())
for name,h in b['outputs'].items():assert sha(V/name)==h
O.mkdir(exist_ok=True);w=openpyxl.load_workbook(V/'S1_reconciled_tables_HELD.xlsx',read_only=True,data_only=True)
font=ImageFont.truetype(font_manager.findfont('DejaVu Sans'),22);head=ImageFont.truetype(font_manager.findfont('DejaVu Sans:style=normal:weight=bold'),22);titlefont=ImageFont.truetype(font_manager.findfont('DejaVu Sans'),30);proofs=[]
def frame(name):
 it=w[name].iter_rows(values_only=True);h=next(it);return [dict(zip(h,x)) for x in it]
def wrap(value,width,f=font):
 text='' if value is None else str(value);out=[]
 for para in text.split('\n'):
  line=''
  for word in para.split():
   trial=(line+' '+word).strip()
   if f.getlength(trial)<=width:line=trial
   else:
    if line:out.append(line);line=''
    while f.getlength(word)>width:
     i=max(1,int(len(word)*width/f.getlength(word)))
     while f.getlength(word[:i])>width:i-=1
     out.append(word[:i]);word=word[i:]
    line=word
  out.append(line)
 return out or ['']
def render(stem,title,rows,cols,widths,note):
 pages=[];page=[];used=160
 for row in rows:
  cells=[wrap(row.get(k),width-22) for k,width in zip(cols,widths)];height=max(52,28*max(map(len,cells))+18)
  if page and used+height>1250:pages.append(page);page=[];used=160
  page.append((cells,height));used+=height
 if page:pages.append(page)
 for j,page in enumerate(pages,1):
  width=sum(widths)+60;foot=wrap(note,width-70);height=160+sum(x[1] for x in page)+len(foot)*27+35;im=Image.new('RGB',(width,height),'white');d=ImageDraw.Draw(im);d.text((30,20),f'{title} ({j}/{len(pages)})',font=titlefont,fill='#28495C');x=30
  for col,cw in zip(cols,widths):
   d.rectangle((x,75,x+cw,150),fill='#304F63')
   for k,t in enumerate(wrap(col,cw-20,head)):d.text((x+10,83+28*k),t,font=head,fill='white')
   x+=cw
  y=150
  for rid,(cells,rh) in enumerate(page):
   x=30
   for cell,cw in zip(cells,widths):
    d.rectangle((x,y,x+cw,y+rh),fill='#F1F5F8' if rid%2==0 else 'white',outline='#D9E1E5')
    for k,t in enumerate(cell):d.text((x+10,y+9+28*k),t,font=font,fill='#22313A')
    x+=cw
   y+=rh
  for k,t in enumerate(foot):d.text((30,y+18+27*k),t,font=font,fill='#4D5960')
  path=O/f'{stem}_{j:02}.png';im.save(path);proofs.append({'file':str(path),'sha256':sha(path),'rows':len(page),'columns':cols,'width':width,'height':height})
rows=[]
for r in frame('Conditional_current'):
 rows.append({'Region':r['region'],'Family':r['family'].replace('_',' '),'Tested variant':r['ID'],'Conditioned on':r['conditioning_IDs'].replace(';',';\n'),'OR (95% Wald CI)':f"{r['OR']:.3f} ({r['CI95_Wald_low']:.3f}-{r['CI95_Wald_high']:.3f})",'REGENIE P':f"{r['P_REGENIE']:.3e}"})
render('conditional','Current conditional estimates',rows,list(rows[0]),[230,260,320,450,270,190],'Exploratory, overlapping discovery analyses. Separate family thresholds: strongest-lead 0.05/9; mutual 0.05/13. Effects per ALT-coded dosage unit. Selected columns; exact output provenance remains in S1.')
rows=[]
for r in frame('Conditional_model_QC'):
 rows.append({'Model':r['analysis'],'Requested / used':f"{r['conditional_variants_requested']} / {r['conditional_variants_used']}",'Numeric covariates':r['original_numeric_covariates'],'Failed tests':r['failed_tests'],'Low-MAC ignored':r['low_MAC_ignored'],'Output rows':r['output_rows']})
render('diagnostics','Conditional model diagnostics',rows,list(rows[0]),[460,200,230,190,230,190],'All conditioning variants were used. Exact test membership verified. The expected extract/exclude notice is preserved in the source table. All 13 unconditioned beta/SE/LOG10P values reproduced exactly.')
rows=[]
for r in frame('Conditional_LD'):
 rows.append({'Variant 1':r['ID1'],'Variant 2':r['ID2'],'Ancestry':r['ancestry'],'Dosage r2':f"{r['r2']:.4f}" if r['available'] else 'Unavailable','Availability / reason':'Available' if r['available'] else r['reason']})
render('LD','Ancestry-specific control LD',rows,list(rows[0]),[350,350,180,210,450],'Current-primary provider-pruned controls. Pair N >=500, both MAC >=20, both call rates >=0.98 and nonzero variance. Exact private denominators are retained in S1; unavailable values are not zero.')
text=(V/'Manuscript_sensitivity_insertions.private.md').read_text()
paragraphs=[p for p in text.split('\n\n') if p.strip()]
merged=[];k=0
while k<len(paragraphs):
 p=paragraphs[k]
 if p.strip().endswith('insertion') and k+1<len(paragraphs):
  p+='\n\n'+paragraphs[k+1];k+=1
 merged.append(p);k+=1
rows=[{'Manuscript insertion candidate':p} for p in merged]
render('prose_final','Manuscript integration text',rows,['Manuscript insertion candidate'],[1500],'Secure review candidate. Native manuscript integration and author approval remain open.')
manifest={'status':'rendered_pending_visual_review','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'workbook_sha256':sha(V/'S1_reconciled_tables_HELD.xlsx'),'build_sha256':sha(V/'BUILD.private.json'),'script_sha256':sha(__file__),'proofs':proofs,'scope':'Selected columns of all new rows and all manuscript insertion text; not native Excel/Word rendering or all-column layout certification.'};(O/'RENDER_MANIFEST.private.json').write_text(json.dumps(manifest,indent=2));print('V7_PROOFS',len(proofs));print('\n'.join(Path(x['file']).name for x in proofs))
