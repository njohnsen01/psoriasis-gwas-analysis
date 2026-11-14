from update_forests_v8 import *
import string,re

def reuse_regions():
 specs=[]
 for letter,(vid,lab) in zip('abcd',XL):
  rs=lab.split()[-1];stem=('NLGN4X_' if 'NLGN4X' in lab else 'PTCHD1-AS_' if 'PTCHD1-AS' in lab else 'Xq12_')+rs
  specs.append((I/'regional_X/figures_v2'/f'{stem}_male02','Supplementary_Figure_S2'+letter,lab+' primary-dosage regional diagnostics','Full regional QC-passing association results are shown with ancestry- and sex-specific control dosage LD. Male X dosage is 0/2. Both Xq12 target panels represent one region. A QC-excluded historical target is a reference position only; it has no primary association estimate. Unavailable LD is shown separately from low LD. These plots describe regional association and correlation, not fine-mapping, colocalization, signal independence or replication.'))
  specs.append((I/'regional_X/figures_v2'/f'{stem}_male01','Supplementary_Figure_S5'+letter,lab+' alternative-dosage regional diagnostics','Full regional association and control dosage-LD results use male 0/1 sensitivity coding. The primary male 0/2 analysis remains the primary result. Historical QC-excluded targets are reference positions, not discoveries. Unavailable LD is distinguished from low LD. The two Xq12 panels describe one region. These regional plots do not establish causal variants, conditional independence, colocalization or independent replication.'))
 specs.append((I/'regional_QKI/figures_v2/QKI_full_region_combined','Supplementary_Figure_S2e','QKI full regional association and ancestry-specific LD','Every finite frozen primary association in the QKI region is shown. Panels use EUR, AFR and AMR control dosage LD with both sexes combined, relative to the specified QKI target. The same pooled association P values are reused in all panels. Missing LD remains visible. QKI was selected from overlapping discovery data and remains exploratory. Dashed and dotted lines mark P = 5 × 10⁻⁸ and 5 × 10⁻⁹. Regional correlation does not establish causal variants, fine-mapping or independent replication.'))
 for letter,sex in [('a','female'),('b','male')]:specs.append((I/f'regional_QKI/figures_v2/QKI_full_region_{sex}','Supplementary_Figure_S6'+letter,'QKI '+sex+'-control LD sensitivity','Panels reuse the same pooled QKI regional association P values and change only the sex-specific control LD reference within EUR, AFR and AMR. These are not sex-specific GWAS panels. Every finite frozen primary regional association is included; unavailable LD remains visible. QKI follow-up is exploratory and uses overlapping participants. Gene proximity and LD do not establish causality, fine-mapping or independent replication.'))
 for src,name,title,legend in specs:
  files={}
  for ext in ['pdf','png','tiff']:
   p=src.with_suffix('.'+ext)
   if not p.exists() and ext=='tiff':p=src.with_name(src.name+'_600dpi').with_suffix('.tiff')
   assert p.exists(),p;remember(p,'verified full regional figure');dst=F/(name+'.'+ext);shutil.copy2(p,dst);assert sha(dst)==sha(p);files[ext]={'file':str(dst),'sha256':sha(dst)}
  assert len(legend.split())<=125
  figures.append({'id':name,'title':title,'legend':legend,'source':str(src),'action':'preserved verified full-region plot; assigned consistent supplement number','files':files})
 print('REGIONS_INTEGRATED',len(specs),flush=True)

def supplementary_calibration(only=None):
 for letter,group in zip('abcde',['female','male','EUR','AFR','AMR']):
  if only is not None and group!=only:continue
  fig,axes=plt.subplots(2,2,figsize=(180/25.4,5.9));fig.subplots_adjust(left=.10,right=.97,bottom=.10,top=.90,hspace=.55,wspace=.42)
  for ax,dom,lab,title in [(axes[0,0],'all','a','Autosomes'),(axes[0,1],'MHC_excluded','b','Autosomes excluding MHC')]:
   root=N/'diagnostics/autosomes_stratified_complete'/group;qq(ax,root/(dom+'_QQ_full_ranks.tsv.gz'),root/(dom+'_calibration.json'),group+': '+title,lab)
  root=N/'diagnostics/X_stratified_complete'/group
  qq(axes[1,0],root/'male02_QQ_full_ranks.tsv.gz',root/'male02_calibration.json',('Female: X dosage 0/1/2' if group=='female' else group+': X male 0/2'),'c')
  if group=='female':axes[1,1].axis('off');axes[1,1].text(.02,.7,'Female dosage is 0/1/2.\nThe identical female model is shown once.\nNo independent replication is implied.',fontsize=8,linespacing=1.8)
  else:qq(axes[1,1],root/'male01_QQ_full_ranks.tsv.gz',root/'male01_calibration.json',group+': X male 0/1 sensitivity','d')
  export_fig(fig,'Supplementary_Figure_S4'+letter,group+' complete-output calibration','Autosomal, MHC-excluded autosomal and X QQ plots use complete valid outputs within the indicated group. Each test has unit weight for QQ ranks and λGC. Display points retain their full-output ranks; significance enrichment does not alter quantiles. Shading is a pointwise 95% beta order-statistic envelope. Male 0/2 is primary and male 0/1 is a sensitivity assumption. The identical female model is displayed once. MHC exclusion uses the recorded chromosome 6 interval. These within-cohort subgroup diagnostics assess calibration and are not independent replication.','Full-output group-specific diagnostic rank and calibration files')
 # Pooled MHC-excluded view is retained separately instead of using the old weighted extract.
 root=N/'diagnostics/pooled_complete'
 paths=list((N/'diagnostics').rglob('*MHC_excluded*QQ_full_ranks*'))
 pooled=[p for p in paths if 'pooled' in str(p)]
 if pooled:
  assert len(pooled)==1
  p=pooled[0];cp=p.with_name(p.name.replace('_QQ_full_ranks.tsv.gz','_calibration.json'))
  fig,ax=plt.subplots(figsize=(87.5/25.4,3.3));fig.subplots_adjust(left=.19,right=.97,bottom=.18,top=.85);qq(ax,p,cp,'Pooled autosomes excluding MHC','a')
  export_fig(fig,'Supplementary_Figure_S4f','Pooled MHC-excluded complete-output calibration','The pooled MHC-excluded autosomal QQ plot uses every valid test with unit weight. Expected quantiles use complete-output ranks. Display thinning does not rerank observations or alter λGC. Shading denotes a pointwise 95% beta order-statistic envelope. The MHC interval follows the frozen diagnostic definition.','Pooled MHC-excluded full-output diagnostics')

def supplementary_y():
 root=N/'Y_source_20260912/results';a=D['Y_male02'];b=D['Y_male01'];a=a.sort_values('ID');b=b.sort_values('ID')
 assert list(a.ID)==list(b.ID) and np.allclose(a.LOG10P,b.LOG10P,rtol=0,atol=1e-10)
 sys.path.insert(0,str(N/'scripts'));from finalize_pooled_calibration import quantiles
 lam,q=quantiles(a.LOG10P.to_numpy());assert abs(lam-float(D['Y_calibration'].iloc[0].lambda_GC))<1e-10
 fig,axs=plt.subplots(1,2,figsize=(180/25.4,3.5));fig.subplots_adjust(left=.09,right=.98,bottom=.19,top=.82,wspace=.4)
 ax=axs[0];ax.scatter(a.GENPOS/1e6,a.LOG10P,s=4,c=TEAL,rasterized=True,linewidths=0);ax.axhline(-np.log10(5e-8),c=ORANGE,ls='--',lw=.8);ax.axhline(-np.log10(.05/len(a)),c=GREY,ls=':',lw=.8);ax.set(xlabel='Chromosome Y position (Mb, GRCh38)',ylabel='−log₁₀(P)',title='Y association scan',ylim=(0,max(8,float(a.LOG10P.max())+.4)));panel(ax,'a')
 ax=axs[1];ax.fill_between(q.expected_log10P,q.null_95_low,q.null_95_high,color='#dce2e8');ax.scatter(q.expected_log10P,q.observed_log10P,s=3,c=TEAL,rasterized=True);mx=max(q.expected_log10P.max(),q.observed_log10P.max())*1.03;ax.plot([0,mx],[0,mx],color=GREY,lw=.7);ax.set(xlim=(0,mx),ylim=(0,mx),xlabel='Expected −log₁₀(P)',ylabel='Observed −log₁₀(P)',title=f'Full-output QQ: λGC = {lam:.3f}');ax.text(.04,.94,f'{len(a):,} valid tests',transform=ax.transAxes,va='top');panel(ax,'b')
 export_fig(fig,'Supplementary_Figure_S3','Corrected chromosome Y association and calibration','(a) Every valid corrected-source chromosome Y association is shown. Dashed and dotted lines denote P = 5 × 10⁻⁸ and the Y-family Bonferroni threshold, respectively. (b) QQ quantiles and λGC use all valid tests with unit weight. Shading is a pointwise 95% beta order-statistic envelope. The male 0/2 and 0/1 fits have identical test P values under rescaling; the primary 0/2 display is shown once. No association meets either threshold. Dosage units remain explicit in S1.','Y_male02; Y_male01; Y_calibration')

def conditional_figures():
 d=D['Conditional_current'];families=list(d.family.unique());assert len(families)==2
 fig,axs=plt.subplots(2,1,figsize=(180/25.4,9.8),gridspec_kw={'height_ratios':[len(d[d.family==x]) for x in families]});fig.subplots_adjust(left=.36,right=.81,top=.94,bottom=.06,hspace=.48)
 for j,fam in enumerate(families):
  z=d[d.family==fam].sort_values(['region','GENPOS']);rows=[]
  for r in z.itertuples():rows.append({'label':variant(r.ID),'OR':r.OR,'low':r.CI95_Wald_low,'high':r.CI95_Wald_high,'p':r.P_REGENIE,'color':TEAL if r.family_significant else GREY,'marker':'s' if r.genome_wide_significant else 'o'})
  forest_rows(axs[j],rows,string.ascii_lowercase[j],fam.replace('_',' ')+' (family '+str(len(z))+')')
 export_fig(fig,'Supplementary_Figure_S7a','Current multi-lead autosomal conditional estimates','Conditional REGENIE estimates use the current primary participants, source genotypes, covariates and LOCO predictions. The unconditioned rerun reproduced all 13 source estimates exactly at printed precision. Bars are 95% Wald CIs; P values come from the corresponding conditional REGENIE LOG10P. Teal points meet the applicable family threshold and squares also meet P < 5 × 10⁻⁸. Strongest-lead-conditioned tests form a nine-test family; mutually conditioned tests form a separate 13-test family. Physical reporting groups do not establish independent signals. Historical QC-excluded Xq12 targets are not reintroduced as primary leads.','Conditional_current; Conditional_model_QC; conditional source verification')
 ld=D['Conditional_LD'];regions=sorted(ld.region.unique(),key=lambda x:tuple(int(v) for v in re.findall(r'\d+',x)));fig,axs=plt.subplots(len(regions),3,figsize=(180/25.4,9.6));fig.subplots_adjust(left=.15,right=.87,bottom=.055,top=.96,hspace=.95,wspace=.70)
 for i,reg in enumerate(regions):
  rr=ld[ld.region==reg];ids=sorted(set(rr.ID1)|set(rr.ID2),key=lambda x:int(x.split(':')[1]))
  for j,anc in enumerate(['EUR','AFR','AMR']):
   z=rr[rr.ancestry==anc];m=np.full((len(ids),len(ids)),np.nan)
   for r in z.itertuples():
    x=ids.index(r.ID1);y=ids.index(r.ID2)
    if r.available:m[x,y]=m[y,x]=r.r2;m[x,x]=m[y,y]=1
   ax=axs[i,j];cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#e2e2e2');imv=ax.imshow(m,vmin=0,vmax=1,cmap=cmap)
   for x in range(len(ids)):
    for y in range(len(ids)):ax.text(y,x,'NA' if np.isnan(m[x,y]) else f'{m[x,y]:.2f}',ha='center',va='center',fontsize=6.2,color='white' if np.isfinite(m[x,y]) and m[x,y]>.65 else '#222222')
   labels=[x.split(':')[1] for x in ids];ax.set_xticks(range(len(ids)));ax.set_xticklabels(labels,rotation=90,fontsize=6.2);ax.set_yticks(range(len(ids)));ax.set_yticklabels(labels,fontsize=6.2);ax.set_title(anc+' '+reg.replace('_',' ').replace('131133','131–133'),fontsize=7.5);panel(ax,string.ascii_lowercase[i*3+j])
 cax=fig.add_axes([.895,.30,.012,.40]);fig.colorbar(imv,cax=cax,label='Dosage r²')
 export_fig(fig,'Supplementary_Figure_S7b','Ancestry-specific LD for current conditional leads','Pair-complete ALT-dosage r² is estimated separately among provider-pruned EUR, AFR and AMR controls for the current multi-lead autosomal regions. Available pairs meet N ≥ 500, both allele counts ≥ 20, call rate ≥ 0.98 and positive dosage variance. Unavailable comparisons are grey and labelled NA, not zero. Axis labels give GRCh38 positions within each named chromosome region. Self-correlation is displayed only where an available pair supports the corresponding variant. LD is descriptive and does not establish independent association signals. Conditional estimates appear in Supplementary Figure S7a.','Conditional_LD')

def sensitivity_figures():
 for letter,sn in [('a','Center_refit'),('b','EHR_sensitivity')]:
  d=D[sn]
  if sn=='Center_refit':d=d[d.model.astype(str).str.contains('02')]
  else:d=d[(d.subset=='provider_pruned')&(d.coding=='male02')]
  fig,axs=plt.subplots(2,2,figsize=(180/25.4,5.8));fig.subplots_adjust(left=.22,right=.88,top=.88,bottom=.13,hspace=.85,wspace=1.8)
  for j,((vid,lab),ax) in enumerate(zip(XL,axs.ravel())):
   z=d[d['ID' if sn=='Center_refit' else 'variant']==vid];assert len(z)>0;rows=[]
   for _,r in z.iterrows():
    if sn=='Center_refit':label={'matched':'Matched baseline','center_refit':'Center adjusted'}[r.adjustment];lo=r.CI95_low;hi=r.CI95_high;p=r.P_from_REGENIE_LOG10P
    else:label={'A_original':'Original','B_observed':'EHR observed','C_no_drug':'EHR / no drug exclusion'}[r.definition];lo=r.CI_low;hi=r.CI_high;p=r.p
    rows.append({'label':label,'OR':r.OR,'low':lo,'high':hi,'p':p,'color':TEAL if j==0 else GREY})
   forest_rows(ax,rows,string.ascii_lowercase[j],lab,'OR per coded ALT unit')
  export_fig(fig,'Supplementary_Figure_S8'+letter,'Sequencing-center sensitivity' if sn=='Center_refit' else 'EHR-control sensitivity','Corrected-source historical X targets are compared under male 0/2 coding. QC-excluded targets remain diagnostic and cannot support a primary discovery claim. Points and bars denote ORs and 95% Wald CIs. Sequencing-center comparisons use matched participants; EHR-control comparisons use provider-pruned participants and show the recorded original, EHR-observed, and EHR-observed without blanket drug-exclusion definitions. P values and intervals are retained from each corresponding source model. These within-cohort sensitivity analyses are not external replication.','Center_refit' if sn=='Center_refit' else 'EHR_sensitivity')
 s=D['SII_sensitivity'];s=s[(s.unit_rules=='strict_units')&(s.subset=='pruned')&(s.status=='full')&(s.male_scale==2)]
 b=D['Biomarkers'];ids=b[(b.subset=='pruned')&(b.status=='full')&(b.male_scale==2)&(b.trait=='SII')&(b.q<.05)&(b.primary_QC_status=='passed')].variant.tolist();assert len(ids)==3
 fig,ax=plt.subplots(figsize=(180/25.4,4.0));fig.subplots_adjust(left=.34,right=.94,bottom=.20,top=.85)
 for i,vid in enumerate(ids):
  z=s[s.variant==vid];assert len(z)==2 and z.n.nunique()==1
  for _,r in z.iterrows():
   old=r.outcome=='original_rebuilt';dy=.12 if old else -.12;c=GREY if old else TEAL
   ax.errorbar(r.beta,i+dy,xerr=[[r.beta-r.CI_low],[r.CI_high-r.beta]],fmt='o' if old else 's',ms=4,color=c,capsize=2,label=('Separate medians' if old else 'Concurrent counts') if i==0 else None)
 ax.set_yticks(range(len(ids)));ax.set_yticklabels([variant(x) for x in ids]);ax.axvline(0,color=GREY,lw=.7);ax.set_xlabel('Standardized log-SII difference per coded target-ALT unit');ax.set_title('Matched SII construction sensitivity');ax.legend(frameon=False,fontsize=7);panel(ax,'a')
 export_fig(fig,'Supplementary_Figure_S9','Concurrent-count versus separate-median SII','For the three SII associations shown in Figure 4, this sensitivity compares the reconstructed separate-median composite with concurrent-count SII under strict unit rules in matched provider-pruned participants. Both estimates for each target use the same analyzed participants. Points and bars show standardized log-SII coefficients and 95% Student’s t CIs from the corresponding models. Outcomes can differ because forming SII before versus after component summarization is not equivalent. Target selection reflects Figure 4 and overlapping data; the complete 25-target sensitivity families and multiplicity corrections remain in S1.','SII_sensitivity; Biomarkers target selection from Figure 4')

if __name__=='__main__':
 try:
  reuse_regions();supplementary_calibration();supplementary_y();conditional_figures();sensitivity_figures()
  savejson(O/'SUPPLEMENT_FIGURES_COMPLETE.private.json',{'figures':figures,'sources':sources,'code_sha256':sha(Path(__file__))})
 except Exception as e:savejson(O/'SUPPLEMENT_FIGURES_FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise

