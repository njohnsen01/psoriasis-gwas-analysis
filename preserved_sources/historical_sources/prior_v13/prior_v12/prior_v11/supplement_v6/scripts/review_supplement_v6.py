from pathlib import Path
from datetime import datetime, timezone
import json,csv,hashlib,math,zipfile,sys
from PIL import Image,ImageDraw,ImageFont
from matplotlib import font_manager
import openpyxl, scipy
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912')
V=I/'supplement_staging_v6'; O=V/'review_proofs'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
 return h.hexdigest()
r=json.loads((V/'BUILD.private.json').read_text())
assert r['all_preexisting_cells_and_styles_preserved']
for name,digest in r['output_hashes'].items(): assert sha(V/name)==digest
assert not O.exists(); O.mkdir()
w=openpyxl.load_workbook(V/'S1_reconciled_tables_HELD.xlsx',read_only=True,data_only=True)
def frame(name):
 it=w[name].iter_rows(values_only=True); h=next(it); return [dict(zip(h,row)) for row in it]
assets={str(x['source']):x['sha256'] for x in frame('Regional_assets')}
inventory=[]
for folder,script in [('regional_X','plot_full_X_regions_v2.py'),('regional_QKI','plot_QKI_regions_v2.py')]:
 for p in sorted((I/folder/'figures_v2').glob('*.tiff')):
  png=p.with_suffix('.png'); pdf=p.with_suffix('.pdf'); assert png.is_file() and pdf.is_file()
  assert pdf.read_bytes()[:5]==b'%PDF-'
  with Image.open(p) as im:
   width,height=im.size; dpi=tuple(float(v) for v in im.info.get('dpi',(0,0))); mode=im.mode; compression=str(im.info.get('compression'))
   assert len(dpi)==2 and all(abs(v-600)<.1 for v in dpi)
   im.verify()
  prior=assets.get(str(p)); digest=sha(p)
  if prior: assert prior==digest
  inventory.append({'figure':p.stem,'TIFF_source':str(p),'TIFF_sha256':digest,'PDF_source':str(pdf),'PDF_sha256':sha(pdf),'PNG_source':str(png),'PNG_sha256':sha(png),'width_pixels':width,'height_pixels':height,'dpi_x':dpi[0],'dpi_y':dpi[1],'width_mm':width/dpi[0]*25.4,'height_mm':height/dpi[1]*25.4,'color_mode':mode,'compression':compression,'prior_S1_TIFF_hash_match':bool(prior),'plot_script':str(I/'scripts'/script),'plot_script_sha256':sha(I/'scripts'/script),'interpretation':'Pooled association plotted in every panel; control LD reference varies. Not independent sex/ancestry GWAS or fine-mapping.','release':'HELD'})
assert len(inventory)==11
p=V/'Regional_figure_inventory.private.tsv'
with p.open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(inventory[0]),delimiter='\t'); writer.writeheader(); writer.writerows(inventory)
legends='''# Full regional figure legends for integration

Chromosome X regional plots. The eight figure files cover four historical target variants in three genomic regions, displayed under the primary male 0/2 coding and the male 0/1 sensitivity. Each figure displays the same complete QC-passing pooled regional REGENIE association results across six panels; panels differ in the EUR, AFR or AMR and female or male provider-pruned control reference used to calculate pairwise dosage r-squared. These are LD-reference panels, not six independent association scans. Coordinates use GRCh38. Gray points indicate LD unavailable after the allele-count or call-rate gates. The dotted vertical line marks the reference variant position. Horizontal lines mark P = 5 x 10^-8 and P = 5 x 10^-9. A target excluded by primary QC has no primary association point and is identified as excluded in the figure. Physical reporting windows and pairwise LD do not demonstrate signal independence, fine-mapping, colocalization or a causal gene. Male 0/1 sensitivity results do not replace the primary male 0/2 analysis.

QKI regional plots. The combined-control, female-control and male-control figure files show the complete pooled association region around chr6:163803171:CTTATT:C (GRCh38). EUR, AFR and AMR panels reuse the same pooled association values and vary only the control LD reference. Point colors show pairwise dosage r-squared; unavailable LD is explicit. Horizontal lines mark P = 5 x 10^-8 and P = 5 x 10^-9. These are descriptive regional plots. Within-cohort subsets are not independent replication, and the figures do not establish fine-mapping or causal assignment to QKI.

The inventory records each existing TIFF, PDF and review PNG, source plotting script, checksum, pixel dimensions and stored resolution. TIFF metadata verification is not final journal acceptance or combined disclosure certification. Existing visual-QA records remain applicable to byte-identical figures. Final figure numbering and placement must be reconciled to the author-approved manuscript scope.
'''
(V/'Regional_figure_legends.private.md').write_text(legends)
regular=font_manager.findfont('DejaVu Sans'); bold=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
font=ImageFont.truetype(regular,18); head=ImageFont.truetype(bold,18); titlefont=ImageFont.truetype(bold,28)
measure=ImageDraw.Draw(Image.new('RGB',(1,1)))
def wrap(value,width,f=font):
 lines=[]
 for paragraph in str(value if value is not None else 'Unavailable').split('\n'):
  line=''
  for word in paragraph.split(' '):
   if measure.textlength(line+(' ' if line else '')+word,font=f)<=width:
    line+=(' ' if line else '')+word
   else:
    if line: lines.append(line)
    line=''
    for char in word:
     if measure.textlength(line+char,font=f)>width: lines.append(line); line=''
     line+=char
  lines.append(line)
 return lines or ['']
proofs=[]
def render(rows,cols,widths,title,stem,note):
 pages=[]; page=[]; used=150
 for row in rows:
  wrapped=[wrap(row.get(k),width-20) for k,width in zip(cols,widths)]
  height=max(48,26*max(map(len,wrapped))+16)
  if used+height>1030 and page: pages.append(page);page=[];used=150
  page.append((wrapped,height));used+=height
 if page: pages.append(page)
 for i,page in enumerate(pages,1):
  h=max(650,150+sum(h for _,h in page)+115); im=Image.new('RGB',(sum(widths)+60,h),'white');d=ImageDraw.Draw(im)
  d.text((30,20),f'{title} ({i}/{len(pages)})',font=titlefont,fill='#28495C');x=30;y=72
  for col,width in zip(cols,widths):
   d.rectangle((x,y,x+width,y+72),fill='#304F63')
   for j,line in enumerate(wrap(col,width-18,head)): d.text((x+9,y+8+j*22),line,font=head,fill='white')
   x+=width
  y=144
  for idx,(cells,height) in enumerate(page):
   x=30
   for text,width in zip(cells,widths):
    d.rectangle((x,y,x+width,y+height),fill='#F1F5F8' if idx%2==0 else 'white',outline='#D9E1E5')
    for j,line in enumerate(text): d.text((x+9,y+8+j*26),line,font=font,fill='#22313A')
    x+=width
   y+=height
  for j,line in enumerate(wrap(note,im.width-60)): d.text((30,y+18+j*23),line,font=font,fill='#4D5960')
  path=O/f'{stem}_{i:02}.png';im.save(path);proofs.append({'file':str(path),'sha256':sha(path),'rows':len(page),'columns':cols})
loo=frame('Ancestry_LOO_current'); show=[]
for x in loo:
 d={k:x[k] for k in ['subset','variant','ancestry_left_out','status','primary_QC_status']}
 d['OR (95% CI)']=f"{x['OR']:.3f} ({x['CI_low']:.3f}, {x['CI_high']:.3f})" if x['status']=='estimated' else 'Unavailable'
 d['Wald P']=f"{x['p_Wald']:.2e}" if x['status']=='estimated' else 'Unavailable'
 d['Bonferroni P']=f"{x['p_Bonferroni']:.2e}" if x['status']=='estimated' else 'Unavailable'; show.append(d)
render(show,list(show[0]),[100,295,150,160,160,240,160,170],'Corrected autosomal leave-one-group-out sensitivity','LOO','Post-hoc, within-cohort sensitivity. All 132 attempted contrasts define the correction family. Failed-QC targets remain diagnostic; unavailable rows are not negative results. No participant counts shown in these proof columns.')
ann=frame('Functional_lead_audit'); cols=['variant','current_role','original_sheet','exact_ID_match','VEP_gene_symbols','source_annotation_status']
render(ann,cols,[300,280,265,130,250,300],'Exact-identifier annotation source review','Annotation','Matching preserved source text is not a fresh annotation query. No exact match means unavailable in the two reviewed tables. Complete source text and hashes remain in the workbook.')
cond=frame('Conditional_source_audit');cols=['analysis','command_options','review_status','interpretation']
render(cond,cols,[280,650,330,425],'Historical conditional-model source audit','Conditional','These historical commands do not establish current relatedness-aware conditional independence. No participant data appear in this proof.')
figrows=[{'figure':x['figure'],'size (pixels)':f"{x['width_pixels']} x {x['height_pixels']}",'DPI':f"{x['dpi_x']:.0f} x {x['dpi_y']:.0f}",'prior TIFF hash':x['prior_S1_TIFF_hash_match'],'formats':'TIFF, PDF, PNG'} for x in inventory]
render(figrows,list(figrows[0]),[640,260,150,240,220],'Existing regional figure inventory','Figures','Existing source files were inspected and hashed, not redrawn. Final numbering, manuscript placement and combined disclosure review remain open.')
text=(V/'S2_reconciliation_integration.private.md').read_text().split('## Corrected leave-one-group-out sensitivity and historical source reconciliation',1)[1]
paragraphs=[{'text':s.strip()} for s in (text+'\n\n'+legends).split('\n\n') if s.strip()]
render(paragraphs,['text'],[1650],'Supplement and legend text additions','Text','Text proofs only. Native Word/Excel layout and final manuscript integration remain separate checks.')
record={'created_UTC':datetime.now(timezone.utc).isoformat(),'build_receipt_sha256':sha(V/'BUILD.private.json'),'workbook_sha256':sha(V/'S1_reconciled_tables_HELD.xlsx'),'figure_inventory_sha256':sha(V/'Regional_figure_inventory.private.tsv'),'figure_count':len(inventory),'TIFF_dpi_checked':600,'existing_S1_TIFF_hash_matches':sum(x['prior_S1_TIFF_hash_match'] for x in inventory),'proofs':proofs,'status':'rendered_visual_inspection_pending','Python':sys.version,'scipy':scipy.__version__,'script_sha256':sha(__file__)}
(O/'RENDER_MANIFEST.private.json').write_text(json.dumps(record,indent=2))
print('RENDERED',len(proofs),'proofs; FIGURES',len(inventory),'TIFF HASH MATCHES',record['existing_S1_TIFF_hash_matches'])
print('PROOF_FILES',[Path(x['file']).name for x in proofs])

