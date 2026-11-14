from pathlib import Path
import json,hashlib,traceback
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
from PIL import Image
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';D=I/'regional_X';O=D/'figures';O.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':7.5,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def main():
 checkpoint=json.loads((D/'COMPLETE.private.json').read_text());assert checkpoint['all_QC_passing_regional_tests_present'];r=pd.read_csv(D/'full_regional_associations.private.tsv.gz',sep=chr(9));l=pd.read_csv(D/'ancestry_sex_control_LD.private.tsv.gz',sep=chr(9));q=pd.read_csv(D/'all_regional_QC.private.tsv.gz',sep=chr(9));records=[]
 names={'chrX:5998724:A:G':'rs183018930','chrX:22917619:C:T':'rs5970952','chrX:67736589:T:C':'rs5964607','chrX:67842269:G:T':'rs5919457'}
 for region,(lo,hi,refs) in checkpoint['windows'].items():
  ymax=max(9.,float(r.loc[r.region.eq(region),'LOG10P'].max())*1.12)
  for ref in refs:
   refq=q.loc[q.ID.eq(ref)];assert len(refq)==1;passed=bool(refq.iloc[0].PASS);refx=int(ref.split(':')[1])/1e6
   for model in ['male02','male01']:
    dat=r.loc[r.region.eq(region)&r.model.eq(model)].sort_values('GENPOS');assert not dat.ID.duplicated().any();assert dat.LOG10P.notna().all();fig,axs=plt.subplots(2,3,figsize=(180/25.4,6.1),sharex=True,sharey=True);fig.subplots_adjust(left=.09,right=.985,bottom=.23,top=.80,wspace=.17,hspace=.29);title=region+' regional association';scale='primary male 0/2' if model=='male02' else 'male 0/1 sensitivity'
    fig.suptitle(title,fontsize=11,fontweight='bold',y=.975);fig.text(.5,.935,'Pooled REGENIE P values | '+scale,ha='center',fontsize=8);fig.text(.5,.899,'LD reference: '+names[ref]+' | '+('passed primary QC' if passed else 'excluded by primary HWE gate'),ha='center',fontsize=8)
    for j,sex in enumerate(['female','male']):
     for k,anc in enumerate(['EUR','AFR','AMR']):
      ax=axs[j,k];z=l.loc[l.region.eq(region)&l.reference.eq(ref)&l.ancestry.eq(anc)&l.sex.eq(sex)];assert not z.variant.duplicated().any();a=dat.merge(z[['variant','r2','status']],left_on='ID',right_on='variant',how='left',validate='one_to_one');assert len(a)==len(dat) and a.status.notna().all();ok=a.status.eq('available');ax.scatter(a.loc[~ok,'GENPOS']/1e6,a.loc[~ok,'LOG10P'],s=5,c='#b3b3b3',linewidths=0,zorder=2);ax.scatter(a.loc[ok,'GENPOS']/1e6,a.loc[ok,'LOG10P'],s=6,c=a.loc[ok,'r2'],norm=Normalize(0,1),cmap='viridis',linewidths=0,zorder=3);ax.axvline(refx,color='#333333',lw=.7,ls=':',zorder=1);ax.axhline(-np.log10(5e-8),color='#aa3939',lw=.8,ls='--');ax.axhline(-np.log10(5e-9),color='#444444',lw=.65,ls='--');v=a.loc[a.ID.eq(ref)]
      if passed:assert len(v)==1;ax.scatter([refx],v.LOG10P,s=25,marker='D',facecolors='none',edgecolors='black',linewidths=.8,zorder=5)
      else:assert len(v)==0
      ax.set_xlim(lo/1e6,hi/1e6);ax.set_ylim(0,ymax);ax.set_title(anc+' | '+sex+' control LD',pad=7);ax.text(-.04,1.08,'abcdef'[j*3+k],transform=ax.transAxes,fontweight='bold',fontsize=9);ax.xaxis.set_major_locator(MaxNLocator(3));ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'));ax.yaxis.set_major_locator(MaxNLocator(integer=True,nbins=5));ax.tick_params(length=2)
      if j==1:ax.set_xlabel('Chromosome X (Mb, GRCh38)')
      if k==0:ax.set_ylabel(r'$-\log_{10}(P)$')
      records.append({'figure':region+'_'+names[ref]+'_'+model,'ancestry':anc,'sex':sex,'points':len(a),'source_records':len(dat),'all_passing_points_present':True,'reference_passed':passed,'ld_unavailable_points':int((~ok).sum())})
    cbax=fig.add_axes([.30,.15,.40,.02]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=cbax,orientation='horizontal');cb.set_label('Within-sex control dosage correlation ($r^2$)',labelpad=2);cb.set_ticks([0,.2,.4,.6,.8,1]);fig.text(.09,.087,'Gray: LD unavailable after count/call-rate gates. Dotted vertical line: reference position.',fontsize=6.8);fig.text(.09,.06,'Dashed lines: P = 5 × 10⁻⁸ (red), 5 × 10⁻⁹ (gray). Each panel contains every passing regional test.',fontsize=6.8);fig.text(.09,.033,'LD uses provider-pruned controls. Failed-QC leads have no primary P point. No fine-mapping claim.',fontsize=6.8)
    stem=region+'_'+names[ref]+'_'+model;fig.savefig(O/(stem+'.png'),dpi=180);fig.savefig(O/(stem+'.pdf'));fig.savefig(O/(stem+'.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
    with Image.open(O/(stem+'.tiff')) as im:assert abs(im.size[0]/600*25.4-180)<.1 and abs(im.info['dpi'][0]-600)<.1
 pd.DataFrame(records).to_csv(O/'plot_point_QA.private.tsv',sep=chr(9),index=False);files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in O.iterdir() if p.suffix in ['.png','.pdf','.tiff']};(O/'MANIFEST.private.json').write_text(json.dumps({'status':'rendered; visual review pending','files':files,'source_hashes':{p.name:sha(p) for p in [D/'full_regional_associations.private.tsv.gz',D/'ancestry_sex_control_LD.private.tsv.gz',D/'all_regional_QC.private.tsv.gz']},'script_sha256':sha(Path(__file__)),'release':'HELD'},indent=2));print('REGIONAL_FIGURES_RENDERED',flush=True)
try:main()
except Exception as e:(O/'FAILED.private.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc()}));raise
