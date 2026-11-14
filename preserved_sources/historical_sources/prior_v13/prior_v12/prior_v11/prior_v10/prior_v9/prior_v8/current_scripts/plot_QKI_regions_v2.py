from pathlib import Path
import json,hashlib,math
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.lines import Line2D
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');O=N/'reconciliation_20260912/regional_QKI';F=O/'figures_v2';F.mkdir(exist_ok=True)
assert json.loads((O/'COMPLETE.private.json').read_text())['status']=='verified'
a=pd.read_csv(O/'full_associations.private.tsv.gz',sep=chr(9));ld=pd.read_csv(O/'ancestry_control_LD.private.tsv.gz',sep=chr(9));assert np.isfinite(a.LOG10P).all() and a.ID.is_unique
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'pdf.fonttype':42,'ps.fonttype':42})
colors=['#173c75','#1f9bbb','#46a965','#efa82c','#cb3335'];bounds=[0,.2,.4,.6,.8,1.000001];cmap=ListedColormap(colors);norm=BoundaryNorm(bounds,len(colors));receipts=[];ref='chr6:163803171:CTTATT:C'
for sex in ['combined','female','male']:
 fig,axes=plt.subplots(3,1,figsize=(180/25.4,175/25.4),sharex=True,sharey=True);fig.subplots_adjust(left=.105,right=.98,bottom=.22,top=.88,hspace=.30)
 ymax=max(9,math.ceil(a.LOG10P.max()+.8))
 for j,(ax,anc) in enumerate(zip(axes,['EUR','AFR','AMR'])):
  l=ld.loc[ld.sex.eq(sex)&ld.ancestry.eq(anc)];assert l.variant.is_unique and set(l.variant)==set(a.ID)
  z=a.merge(l,left_on='ID',right_on='variant',validate='one_to_one');ok=z.status.eq('available');assert len(z)==len(a);x=z.GENPOS/1e6
  ax.scatter(x[~ok],z.LOG10P[~ok],s=5,c='#c7c7c7',edgecolors='none',rasterized=True,zorder=2)
  ax.scatter(x[ok],z.LOG10P[ok],s=5,c=z.loc[ok,'r2'],cmap=cmap,norm=norm,edgecolors='none',rasterized=True,zorder=3)
  lead=z.loc[z.ID.eq(ref)];assert len(lead)==1;ax.scatter(lead.GENPOS/1e6,lead.LOG10P,s=30,marker='D',facecolor='#8c4a9e',edgecolor='#222222',linewidth=.5,zorder=5)
  ax.axvline(163.803171,color='#8c4a9e',lw=.65,alpha=.8);ax.axhline(-np.log10(5e-8),color='#b43d3d',lw=.7,ls='--');ax.axhline(-np.log10(5e-9),color='#666666',lw=.6,ls=':')
  ax.set_title(chr(97+j)+') '+anc+' control LD',loc='left',pad=3,fontweight='bold');ax.set_ylabel(r'$-\log_{10}(P)$');ax.set_ylim(0,ymax);ax.set_xlim(163.303171,164.303171);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#ececec',lw=.5)
  receipts.append({'sex':sex,'ancestry':anc,'displayed_records':len(z),'all_finite_frozen_tests_shown':True,'missing_LD_points_visible':True})
 axes[-1].set_xlabel('Chromosome 6 position (Mb, GRCh38)')
 fig.suptitle('QKI region: pooled psoriasis association and ancestry-specific LD',x=.105,ha='left',y=.975,fontweight='bold',fontsize=10)
 fig.text(.105,.933,'LD reference: chr6:163803171:CTTATT:C  |  Controls: '+('both sexes' if sex=='combined' else sex),fontsize=8)
 handles=[Line2D([0],[0],marker='o',lw=0,color=c,label=l,markersize=4) for c,l in zip(colors,['0–0.2','0.2–0.4','0.4–0.6','0.6–0.8','0.8–1.0'])];handles+=[Line2D([0],[0],marker='o',lw=0,color='#c7c7c7',label='Unavailable',markersize=4),Line2D([0],[0],marker='D',lw=0,color='#8c4a9e',label='Reference',markersize=4)]
 fig.legend(handles=handles,ncol=7,loc='lower center',bbox_to_anchor=(.53,.075),frameon=False,title=r'Pairwise dosage $r^2$',fontsize=7,title_fontsize=7,handletextpad=.3,columnspacing=.7)
 fig.text(.105,.035,'Each panel shows the same complete pooled association region; only the LD reference group changes.',fontsize=7)
 fig.text(.105,.014,'Dashed: P = 5 × 10⁻⁸; dotted: P = 5 × 10⁻⁹. Descriptive regional plot; no fine-mapping claim.',fontsize=7)
 stem=F/('QKI_full_region_'+sex)
 for ext in ['png','pdf','tiff']:
  kwargs={'dpi':600}
  if ext=='tiff':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(str(stem)+'.'+ext,**kwargs)
 plt.close(fig)
(F/'PLOT_QA.private.json').write_text(json.dumps({'status':'numerical_plot_checks_passed','panels':receipts,'visual_review':'pending','font':'DejaVu Sans; final journal font check pending','release':'HELD'},indent=2));print('QKI_FIGURES_NUMERICALLY_VERIFIED',flush=True)
