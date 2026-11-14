from pathlib import Path
import json,hashlib,shutil,sys,math,traceback,importlib.metadata as im
import numpy as np,pandas as pd
from scipy import stats
from openpyxl import load_workbook
from openpyxl.styles import Font,PatternFill,Alignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from PIL import Image
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
I=N/'reconciliation_20260912'; SRC=I/'supplement_staging_v7'
O=I/'tables_figures_v8'; O.mkdir(exist_ok=True)
F=O/'figures'; T=O/'tables'; P=O/'review_proofs'
for d in [F,T,P]:d.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def savejson(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
w=load_workbook(SRC/'S1_reconciled_tables_HELD.xlsx',read_only=True,data_only=True)
D={}
for s in w:
 rr=list(s.values);D[s.title]=pd.DataFrame(rr[1:],columns=rr[0])
w.close()
BLUE='#315d82'; TEAL='#007c83'; ORANGE='#b76a32'; PURPLE='#825a94'; GREY='#69727b'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
figures=[];tables=[];sources=[]
def remember(p,kind):
 p=Path(p);sources.append({'file':str(p),'sha256':sha(p),'role':kind});return p
remember(SRC/'S1_reconciled_tables_HELD.xlsx','reconciled source workbook')
def fp(p):
 if p is None or not np.isfinite(float(p)):return 'Unavailable'
 p=float(p)
 if p==0:return '<1e-300'
 return f'{p:.2e}' if p<.001 else f'{p:.3f}'
def ci(v,lo,hi):return f'{float(v):.3f} ({float(lo):.3f}, {float(hi):.3f})'
def variant(v):return str(v).replace('chr','')
def panel(ax,label):ax.text(-min(.10,.065/ax.get_position().width),1.055,label,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def export_fig(fig,name,title,legend,origin):
 assert len(legend.split())<=125,(name,len(legend.split()))
 fig.savefig(F/(name+'.pdf'),facecolor='white')
 fig.savefig(F/(name+'.png'),dpi=180,facecolor='white')
 high=F/(name+'_600dpi.png');fig.savefig(high,dpi=600,facecolor='white')
 with Image.open(high) as z:z.convert('RGB').save(F/(name+'.tiff'),compression='tiff_lzw',dpi=(600,600))
 high.unlink();plt.close(fig)
 figures.append({'id':name,'title':title,'legend':legend,'source':origin,'action':'rebuilt from current source','files':{e:{'file':str(F/(name+'.'+e)),'sha256':sha(F/(name+'.'+e))} for e in ['pdf','png','tiff']}})
 print('FIGURE',name,flush=True)
def table_file(name,title,d,footnote,source):
 d.to_csv(T/(name+'.tsv'),sep='\t',index=False)
 tables.append({'id':name,'title':title,'data':d,'footnote':footnote,'source':source})

def verify_display_sources():
 seen={}
 for sn,fc,hc in [('Exact_primary_estimates','source_file','source_sha256'),('X_followup','source_output','source_sha256'),('Ancestry_followup','source_output','source_sha256'),('Biomarkers','source_output','source_sha256'),('SII_sensitivity','source_output','source_sha256'),('Center_refit','source_file','source_sha256'),('Conditional_current','output','output_sha256')]:
  for f,h in D[sn][[fc,hc]].drop_duplicates().itertuples(index=False,name=None):
   if pd.isna(f) or pd.isna(h):continue
   if f not in seen:
    assert sha(f)==h,('Changed display source',sn,f);seen[f]=h
 savejson(O/'DISPLAY_SOURCE_VERIFICATION.private.json',{'verified_files':len(seen),'files':seen})
 print('DISPLAY_SOURCES_VERIFIED',len(seen),flush=True)

def manhattan_data():
 sys.path.insert(0,str(N/'scripts'));from finalize_pooled_calibration import inputs,sha as prior_sha
 records=[];counts={};receipts=[]
 for group in ['autosomes','male02']:
  total=0;valid=0;sig=0
  for file,h,n,e,u in inputs(group):
   assert prior_sha(file)==h
   d=pd.read_csv(file,sep=r'\s+',usecols=['CHROM','GENPOS','ID','BETA','SE','LOG10P'])
   vals=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce')
   ok=np.isfinite(vals).all(axis=1)&vals.SE.gt(0)&vals.LOG10P.ge(0)
   assert len(d)==n and int(ok.sum())==e
   hit=ok&d.LOG10P.gt(-np.log10(5e-8))
   select=ok&(((np.arange(len(d))+total)%40==0)|d.LOG10P.ge(5))
   records.append(d.loc[select,['CHROM','GENPOS','ID','LOG10P']].copy())
   total+=n;valid+=e;sig+=int(hit.sum());receipts.append({'file':str(file),'sha256':h,'valid':e})
  counts[group]={'tested':total,'valid':valid,'gws':sig}
  print('MANHATTAN_SOURCE',group,counts[group],flush=True)
 d=pd.concat(records,ignore_index=True);assert d.ID.is_unique
 d.to_csv(O/'Manhattan_display.private.tsv.gz',sep='\t',index=False)
 savejson(O/'MANHATTAN_SOURCE.private.json',{'counts':counts,'files':receipts,'display':'Every 40th source row plus every finite P <= 1e-5; no change to full-rank QQ or inflation'})
 return d,counts

def qq(ax,qfile,calfile,title,label):
 q=pd.read_csv(remember(qfile,'full-output QQ ranks'),sep='\t');c=json.loads(remember(calfile,'full-output calibration').read_text())
 n=int(q.N_valid_tests.iloc[0]);assert q.N_valid_tests.nunique()==1
 assert np.allclose(q.expected_log10P,-np.log10((q.full_output_rank-.5)/n),rtol=0,atol=1e-12)
 assert q.full_output_rank.is_monotonic_increasing and q.observed_log10P.is_monotonic_decreasing
 lam=c.get('lambda_GC_from_full_output_P',c.get('lambda_GC',c.get('lambda_gc')))
 if lam is None:raise ValueError(('Unknown calibration schema',c.keys()))
 ax.fill_between(q.expected_log10P,q.null_95_low,q.null_95_high,color='#dce2e8',linewidth=0)
 ax.scatter(q.expected_log10P,q.observed_log10P,s=2,c=TEAL if '0/2' in title else BLUE,rasterized=True,linewidths=0)
 mx=q.expected_log10P.max()*1.025;my=max(mx,q.observed_log10P.max()*1.025)
 ax.plot([0,mx],[0,mx],color=GREY,lw=.7);ax.set(xlim=(0,mx),ylim=(0,my),xlabel='Expected −log₁₀(P)',ylabel='Observed −log₁₀(P)',title=title)
 ax.text(.035,.955,f'λGC = {lam:.3f}\nTests = {n:,}',transform=ax.transAxes,ha='left',va='top',fontsize=6.5,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
 panel(ax,label)
 return c

def figure1(d,counts):
 fig=plt.figure(figsize=(180/25.4,6.1));g=fig.add_gridspec(2,3,height_ratios=[1.1,1],hspace=.52,wspace=.57,left=.09,right=.98,bottom=.09,top=.92)
 ax=fig.add_subplot(g[0,:]);offset=0;ticks=[];labs=[]
 for ch in range(1,24):
  z=d[d.CHROM.astype(str).isin([str(ch),'X' if ch==23 else str(ch)])]
  if z.empty:continue
  ax.scatter(offset+z.GENPOS,z.LOG10P,s=1.7,color=TEAL if ch==23 else [BLUE,'#96afc5'][ch%2],rasterized=True,linewidths=0)
  length=max(float(z.GENPOS.max()),1);ticks.append(offset+length/2);labs.append('X' if ch==23 else str(ch));offset+=length+8e6
 ax.axhline(-np.log10(5e-8),color=ORANGE,lw=.65,ls='--');ax.axhline(-np.log10(5e-9),color=GREY,lw=.6,ls=':')
 ax.set(xlim=(0,offset),xticks=ticks,xticklabels=labs,xlabel='Chromosome (GRCh38)',ylabel='−log₁₀(P)',title='Corrected primary association scan')
 plt.setp(ax.get_xticklabels(),rotation=90,ha='center');panel(ax,'a');ax.text(.995,.96,'X: male 0/2',transform=ax.transAxes,ha='right',va='top',color=TEAL)
 qroot=N/'diagnostics/pooled_complete'
 for col,(grp,label,title) in enumerate([('autosomes','b','Autosomes'),('male02','c','X: male 0/2'),('male01','d','X: male 0/1 sensitivity')]):qq(fig.add_subplot(g[1,col]),qroot/(grp+'_QQ_full_ranks.tsv.gz'),qroot/(grp+'_calibration.json'),title,label)
 export_fig(fig,'Figure_1','Corrected association landscape and full-output calibration','(a) Primary autosomal and chromosome X associations use corrected source genotypes and male 0/2 X dosage. Display points are a systematic subset plus all P ≤ 10⁻⁵; the full vertical range is retained. Dashed and dotted lines mark P = 5 × 10⁻⁸ and 5 × 10⁻⁹. (b–d) QQ coordinates and λGC use every valid association with unit weight. Display thinning follows complete-output ranking; points are never reranked. Shading denotes pointwise 95% beta order-statistic bands. Male 0/1 is a sensitivity analysis. No QC-passing X association meets P < 5 × 10⁻⁸ in the primary male 0/2 model.','Frozen primary output receipts and complete-output QQ ranks')

if __name__=='__main__':
 try:
  verify_display_sources();data,counts=manhattan_data();figure1(data,counts)
  savejson(O/'FIGURE1_COMPLETE.private.json',{'figures':figures,'sources':sources,'code_sha256':sha(Path(__file__))})
 except Exception as e:savejson(O/'BUILD_FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise

