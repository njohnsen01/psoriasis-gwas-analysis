import os
from pathlib import Path
import json,hashlib,numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
A=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09')
F=Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}/final_aggregate')
plot=pd.read_csv(F/'psoriasis_corrected_manhattan_plot_data.tsv.gz',sep='\t',dtype={'CHROM':str})
priority=pd.read_csv(F/'psoriasis_corrected_priority_variant_results.tsv',sep='\t',dtype={'CHROM':str})
plot.CHROM=plot.CHROM.str.replace('chr','',regex=False).replace({'23':'X'})
assert len(plot)==136989 and plot.ID.is_unique and (plot.LOG10P>=-np.log10(5e-8)).sum()==4658
S=json.loads((A/'full_calibration_summary.json').read_text());stats={r['group']:r for r in S['groups']}
q=pd.read_csv(A/'exact_qq_points.tsv',sep='\t')
C=[str(i) for i in range(1,23)]+['X'];mx=plot.groupby('CHROM').POS.max().to_dict();off={};cent={};cur=0
for c in C:off[c]=cur;cent[c]=cur+mx[c]/2;cur+=mx[c]+8000000
end=cur-8000000
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'xtick.major.width':.6,'ytick.major.width':.6})
fig=plt.figure(figsize=(7.08,7.35));gs=fig.add_gridspec(2,2,height_ratios=[1.06,1],wspace=.3,hspace=.38)
ax=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,0]);cax=fig.add_subplot(gs[1,1])
for i,c in enumerate(C):
 z=plot[plot.CHROM==c];ax.scatter(z.POS+off[c],np.minimum(z.LOG10P,21.3),s=1.7,c='#C94F70' if c=='X' else ['#A4B7CA','#4C78A8'][i%2],alpha=.72,linewidths=0,rasterized=True)
th=-np.log10(5e-8);ax.axhline(th,color='#6B7280',lw=.8,ls=(0,(4,3)));ax.text(.003,th+.22,r'$P=5\times10^{-8}$',transform=ax.get_yaxis_transform(),fontsize=7,color='#606770')
ax.set(xlim=(-.01*end,1.01*end),ylim=(0,22.5),xlabel='Chromosome',ylabel=r'$-\log_{10}(P)$')
ax.set_xticks([cent[c] for c in C]);ax.set_xticklabels(C,fontsize=6.5,rotation=45,ha='right',rotation_mode='anchor')
ax.text(-.065,1.03,'a',transform=ax.transAxes,fontweight='bold',fontsize=10)
labels=[('IL12B','chr5:159402519:A:T','#24557A',(-28,11),'center'),('TNIP1','chr5:151087129:GC:G','#24557A',(-20,18),'center'),('MHC/HLA','chr6:31299200:G:A','#24557A',(32,11),'center'),('QKI','chr6:163803171:CTTATT:C','#D97706',(18,16),'left'),('TYK2','chr19:10352442:G:C','#24557A',(0,22),'center'),('NLGN4X','chrX:5998724:A:G','#D97706',(-17,8),'center'),('PTCHD1-AS','chrX:22917619:C:T','#D97706',(0,27),'center'),('Xq12','chrX:67842269:G:T','#D97706',(13,15),'right')]
records=[]
for label,vid,color,offset,align in labels:
 z=plot[plot.ID==vid]
 if z.empty:
  z=priority[priority.ID==vid].copy()
  if 'POS' not in z:z=z.rename(columns={'GENPOS':'POS'})
 assert len(z)==1,vid
 r=z.iloc[0];chrom=str(r.CHROM).removeprefix('chr');chrom='X' if chrom=='23' else chrom;x=int(r.POS)+off[chrom];y=min(float(r.LOG10P),21.3)
 marker='o' if vid=='chr6:31299200:G:A' else ('D' if color=='#D97706' else 'o')
 ax.scatter([x],[y],s=25,marker=marker,c=color,edgecolors='white',linewidths=.7,zorder=5)
 ax.annotate(label,(x,y),xytext=offset,textcoords='offset points',ha=align,va='bottom',fontsize=7.3,fontweight='bold',fontstyle='italic' if label not in ['MHC/HLA','Xq12'] else 'normal',arrowprops={'arrowstyle':'-','lw':.7,'color':color,'shrinkA':1,'shrinkB':2})
 records.append({'label':label,'ID':vid,'CHROM':chrom,'POS':int(r.POS),'LOG10P':float(r.LOG10P)})
for axis,key,color,label in [(b,'autosomes','#4C78A8','All autosomal tests'),(b,'autosomes_no_MHC','#14866D','MHC excluded'),(cax,'X','#C94F70','Chromosome X')]:
 z=q[q.group==key];n=stats[key]['n'];v=np.load(A/(key+'_sorted_logp.npy'),mmap_mode='r')
 assert np.array_equal(z.observed.to_numpy(),v[z['rank'].to_numpy()-1]) and np.allclose(z.expected,-np.log10((z['rank']-.5)/n),rtol=0,atol=2e-14)
 axis.scatter(z.expected,z.observed,s=3,alpha=.65,c=color,linewidths=0,rasterized=True,label=label+rf', $\lambda_{{GC}}={stats[key]["lambda_GC"]:.3f}$')
for axis,key,letter,title in [(b,'autosomes','b','Autosomes'),(cax,'X','c','Chromosome X')]:
 xmax=-np.log10(.5/stats[key]['n'])+.2;ymax=max(stats[key]['max_log10p']*1.05,xmax)
 axis.plot([0,xmax],[0,xmax],c='#7C8795',lw=.8,ls=(0,(4,3)));axis.set(xlim=(0,xmax),ylim=(0,ymax),xlabel=r'Expected $-\log_{10}(P)$',ylabel=r'Observed $-\log_{10}(P)$')
 axis.set_title(title+'\n'+f'{stats[key]["n"]:,} valid tests',fontsize=8.5,pad=7)
 axis.text(-.16,1.08,letter,transform=axis.transAxes,fontweight='bold',fontsize=10)
 axis.legend(loc='upper left',frameon=False,fontsize=6.4,handletextpad=.25)
for axis in [ax,b,cax]:axis.spines[['top','right']].set_visible(False);axis.grid(axis='y',color='#D9DDE3',lw=.45,alpha=.7);axis.set_axisbelow(True)
fig.subplots_adjust(left=.09,right=.965,top=.96,bottom=.075)
fig.text(.5,.012,'QQ panels use complete association output, unit weights, and exact ranks before display thinning.',ha='center',fontsize=6.3,color='#606770')
stem=A/'Figure1_Full_Output_Calibration'
fig.savefig(str(stem)+'.png',dpi=300,facecolor='white');fig.savefig(str(stem)+'.pdf',facecolor='white');fig.savefig(str(stem)+'_600dpi.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'},facecolor='white');plt.close(fig)
record={'full_output_sha256':S['sha256'],'manhattan_input_sha256':hashlib.sha256((F/'psoriasis_corrected_manhattan_plot_data.tsv.gz').read_bytes()).hexdigest(),'qq_sampling_weights':False,'all_displayed_observations_equal_full_output_order_statistics':True,'manhattan_label_sources':records,'outputs':[str(stem)+'.png',str(stem)+'.pdf',str(stem)+'_600dpi.tiff']}
(A/'Figure1_Full_Output_Calibration.validation.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
