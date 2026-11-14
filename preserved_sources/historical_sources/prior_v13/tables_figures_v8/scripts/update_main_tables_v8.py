from update_forests_v8 import *
import textwrap,copy
from openpyxl.utils import get_column_letter

def cohort_table():
 flow=D['Cohort_flow'];rows=[]
 labels={'original_analysis_after_PC_source_join':'Original analyzed cohort','provider_coverage_failure':'After provider coverage exclusion','noncanonical_or_unavailable_ploidy':'After noncanonical/unavailable ploidy exclusion','provider_QC_flag':'Corrected primary cohort'}
 for r in flow.itertuples(index=False):
  if r.stage in labels:rows.append(['Cohort stage',labels[r.stage],int(r.N_retained),int(r.cases_retained),int(r.controls_retained)])
 a=D['Age'].set_index('phenotype');e=D['EHR_observation'].set_index('phenotype')
 fmt=lambda r:f'{r["median"]:.1f} ({r["q1"]:.1f}–{r["q3"]:.1f})'
 rows.append(['Characteristics','Age, median (IQR), years',None,fmt(a.loc['case']),fmt(a.loc['control'])])
 rows.append(['Characteristics','Recorded EHR span available',None,int(e.loc['case','N_available']),int(e.loc['control','N_available'])])
 rows.append(['Characteristics','Recorded EHR span unavailable',None,int(e.loc['case','N_missing']),int(e.loc['control','N_missing'])])
 rows.append(['Characteristics','Recorded EHR span, median (IQR), years',None,fmt(e.loc['case']),fmt(e.loc['control'])])
 rows.append(['Characteristics','Recorded EHR span equals zero',None,int(e.loc['case','zero']),int(e.loc['control','zero'])])
 a=D['Ancestry_groups'];a=a[a.population=='strict_complete'];assert a.N.sum()==386069 and a.cases.sum()==5705 and a.controls.sum()==380364
 for group in ['EUR','AFR','AMR']:
  r=a[a.group==group];assert len(r)==1;r=r.iloc[0];rows.append(['Genetic-similarity group',group,int(r.N),int(r.cases),int(r.controls)])
 z=a[~a.group.isin(['EUR','AFR','AMR'])];rows.append(['Genetic-similarity group','Other provider-assigned groups combined',int(z.N.sum()),int(z.cases.sum()),int(z.controls.sum())])
 df=pd.DataFrame(rows,columns=['Section','Measure','Eligible total N','Cases','Controls'])
 df.loc[df.Measure.eq('After noncanonical/unavailable ploidy exclusion'),['Cases','Controls']]='Withheld'
 table_file('Table_1','Corrected cohort flow and characteristics',df,'Cohort stages are nested. Intermediate case/control cells are withheld to prevent small-count disclosure by subtraction. Canonical sex-discordance removed no additional participants. The full upstream exclusion flow and all individual genetic-similarity groups are retained in S1. Group assignment uses the complete echo_v4_r2 ancestry_pred exact labels without a probability cutoff; the follow-up group gate is at least 200 cases and 500 controls, with variant-level QC separate. EHR span summaries use participants with available records and do not establish sufficient observation for every control. Eligible counts are not variant-specific analyzed N.','Cohort_flow; Age; EHR_observation; Ancestry_groups; Followup_groups')

def display_frame(t):
 d=t['data'].copy()
 if t['id']=='Table_1':d=d.drop(columns='Section')
 elif t['id']=='Table_2':
  d['OR (95% Wald CI)']=[ci(v,l,h) if pd.notna(v) else 'Not tested' for v,l,h in d[['OR','95% Wald CI lower','95% Wald CI upper']].itertuples(index=False,name=None)]
  d['P']=[fp(x) if pd.notna(x) else 'QC excluded' for x in d['REGENIE P']]
  d=d[['Variant (GRCh38)','Effect allele','OR (95% Wald CI)','P','Primary status']]
 elif t['id']=='Table_3':
  d['Beta (95% t CI)']=[ci(v,l,h) for v,l,h in d[['Beta','95% t CI lower','95% t CI upper']].itertuples(index=False,name=None)]
  d['P']=[fp(x) for x in d.P];d['FDR q']=[fp(x) for x in d['FDR q']]
  d=d[['Variant (GRCh38)','Effect allele','Biomarker','Beta (95% t CI)','P','FDR q']]
 return d

def render_table(t):
 d=display_frame(t);n=len(d);per=10 if t['id']=='Table_2' else 8 if t['id']=='Table_1' else 10
 pdfpages=[]
 from matplotlib.backends.backend_pdf import PdfPages
 with PdfPages(T/(t['id']+'.pdf')) as pdf:
  for page,start in enumerate(range(0,n,per),1):
   z=d.iloc[start:start+per].copy().fillna('—')
   if t['id']=='Table_1':widths=[.50,.16,.17,.17];wraps=[40,14,20,20]
   elif t['id']=='Table_2':widths=[.28,.09,.26,.13,.24];wraps=[27,9,25,12,22]
   else:widths=[.25,.09,.12,.30,.12,.12];wraps=[24,9,12,26,12,12]
   vals=[]
   for row in z.itertuples(index=False,name=None):
    vals.append(['\n'.join(textwrap.wrap(f'{int(v):,}' if isinstance(v,(int,float,np.integer,np.floating)) and float(v).is_integer() else str(v),w,break_long_words=True,break_on_hyphens=False)) for v,w in zip(row,wraps)])
   labels=['\n'.join(textwrap.wrap(x,w,break_long_words=False)) for x,w in zip(d.columns,wraps)]
   fig,ax=plt.subplots(figsize=(8.5,7.6));ax.axis('off');fig.subplots_adjust(left=.07,right=.96,bottom=.12,top=.88)
   ax.set_title(t['title']+(f' (continued {page})' if page>1 else ''),loc='left',fontsize=11,fontweight='bold',pad=20)
   tbl=ax.table(cellText=vals,colLabels=labels,colWidths=widths,cellLoc='left',colLoc='left',loc='upper left');tbl.auto_set_font_size(False);tbl.set_fontsize(8)
   for (r,c),cell in tbl.get_celld().items():
    cell.set_edgecolor('#d7dde3');cell.set_linewidth(.4);cell.PAD=.08
    if r==0:cell.set_facecolor('#e7edf3');cell.set_text_props(fontweight='bold');cell.set_height(.12)
    else:
     maxlines=max(v.count('\n')+1 for v in vals[r-1]);cell.set_height(.065+.024*(maxlines-1));cell.set_facecolor('white' if r%2 else '#f6f8fa')
   foot=' '.join(t['footnote'].split());fig.text(.08,.07,textwrap.fill(foot,108),fontsize=7,va='bottom',linespacing=1.5)
   fig.savefig(P/f'{t["id"]}_{page:02}.png',dpi=150);pdf.savefig(fig);plt.close(fig);pdfpages.append(str(P/f'{t["id"]}_{page:02}.png'))
 t['proofs']=pdfpages

def build_workbook():
 src=SRC/'S1_reconciled_tables_HELD.xlsx';out=O/'S1_reconciled_tables_HELD.xlsx';book=load_workbook(src)
 existing=book.sheetnames[:]
 for t in sorted(tables,key=lambda t:t['id']):
  assert t['id'] not in book
  ws=book.create_sheet(t['id']);d=t['data'];cols=list(d.columns)
  ws.append([t['title']]);ws.append([t['footnote']]);ws.append(['Source: '+t['source']]);ws.append(cols)
  for row in d.itertuples(index=False,name=None):ws.append([None if pd.isna(x) else x.item() if isinstance(x,np.generic) else x for x in row])
  ws.freeze_panes='C5';ws.auto_filter.ref=f'A4:{get_column_letter(len(cols))}{ws.max_row}';ws.sheet_view.showGridLines=False
  for row in [1,2,3]:ws.merge_cells(start_row=row,start_column=1,end_row=row,end_column=len(cols))
  ws['A1'].font=Font(name='Calibri',size=14,bold=True,color='17365D');ws['A2'].alignment=Alignment(wrap_text=True,vertical='top');ws['A3'].alignment=Alignment(wrap_text=True,vertical='top');ws.row_dimensions[1].height=24;ws.row_dimensions[2].height=78;ws.row_dimensions[3].height=28;ws.row_dimensions[4].height=45
  for c in ws[4]:c.fill=PatternFill('solid',fgColor='DCE6F1');c.font=Font(bold=True,size=10);c.alignment=Alignment(wrap_text=True,vertical='center')
  for j,col in enumerate(cols,1):
   width=40 if ('Variant' in col or 'Measure' in col or 'status' in col.lower()) else 32 if ('group' in col or col=='Section') else 19
   ws.column_dimensions[get_column_letter(j)].width=width
   for r in range(5,ws.max_row+1):
    cell=ws.cell(r,j);cell.alignment=Alignment(wrap_text=True,vertical='top');cell.font=Font(name='Calibri',size=10)
    if isinstance(cell.value,(int,float)):cell.number_format='0.00E+00' if col in ['P','REGENIE P','FDR q'] else '#,##0' if col in ['Eligible total N','Cases','Controls'] else '0.000'
  for r in range(5,ws.max_row+1):ws.row_dimensions[r].height=34
  ws.print_title_rows='1:4';ws.page_setup.orientation='landscape';ws.page_setup.paperSize=ws.PAPERSIZE_A4;ws.page_setup.fitToWidth=1;ws.page_setup.fitToHeight=0;ws.sheet_properties.pageSetUpPr.fitToPage=True
 book.save(out);book.close()
 # Exact all-cell preservation of the preceding 44 worksheets, including styles.
 a=load_workbook(src,read_only=True);b=load_workbook(out,read_only=True);cells=0
 assert b.sheetnames[:len(existing)]==existing
 for sn in existing:
  x=a[sn];y=b[sn];assert (x.max_row,x.max_column)==(y.max_row,y.max_column)
  for xr,yr in zip(x.rows,y.rows):
   for xc,yc in zip(xr,yr):
    assert xc.value==yc.value and getattr(xc,'_style_id',0)==getattr(yc,'_style_id',0),(sn,xc.coordinate);cells+=1
 a.close();b.close()
 savejson(O/'WORKBOOK_VERIFICATION.private.json',{'prior_sheets_preserved':len(existing),'cells_preserved':cells,'new_main_tables':[t['id'] for t in tables],'source_sha256':sha(src),'output_sha256':sha(out),'sheets_total':len(existing)+len(tables),'scope':'Static derived scientific tables; existing 44 worksheets preserved. Native application print rendering and coordinated release disclosure remain separate.'})
 print('WORKBOOK_VERIFIED',len(existing)+len(tables),cells,flush=True)

if __name__=='__main__':
 try:
  checkpoint=json.loads((O/'FORESTS_COMPLETE.private.json').read_text())
  for t in checkpoint['tables']:
   t['data']=pd.DataFrame(t['data']);tables.append(t)
  cohort_table()
  for t in tables:render_table(t)
  build_workbook()
  savejson(O/'TABLES_COMPLETE.private.json',{'tables':[{**t,'data':t['data'].to_dict('records')} for t in tables],'source_workbook':str(SRC/'S1_reconciled_tables_HELD.xlsx'),'code_sha256':sha(Path(__file__))})
 except Exception as e:savejson(O/'TABLES_FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise

