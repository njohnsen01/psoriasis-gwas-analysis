from pathlib import Path
import json,sys,time,traceback,hashlib,datetime,importlib.metadata as im
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from reconcile_X_site_filters import sha,save
from finalize_pooled_calibration import quantiles,inputs
N=Path(__file__).resolve().parents[1]
O=N/'submission_assets_20260911';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600})
GROUPS=['pooled','female','male','EUR','AFR','AMR']
LEADS=[('rs183018930','NLGN4X',5998724,'A','G'),('rs5970952','PTCHD1-AS',22917619,'C','T'),('rs5964607','Xq12',67736589,'T','C'),('rs5919457','Xq12',67842269,'G','T')]
RECEIPTS={};REPORTS=[];LEADROWS=[];REGIONS={}
def sources(group,domain,model='male02'):
 if group=='pooled':
  yield from inputs('autosomes' if domain=='autosomes' else model);return
 root=N/'results'/('autosomes_stratified' if domain=='autosomes' else 'X_stratified')
 assert json.loads((root/'COMPLETE.json').read_text())['status']=='complete'
 paths=sorted(root.glob('*/'+group+'/COMPLETE.json')) if domain=='autosomes' else sorted((root/group).glob('*/COMPLETE.json'))
 assert len(paths)==(1104 if domain=='autosomes' else 134)
 for p in paths:
  r=json.loads(p.read_text());assert r['status']=='validated'
  n=r['passed'] if domain=='autosomes' else r['QC_passed']
  if not n:continue
  m=r if domain=='autosomes' else r['models'][model]
  yield p.parent/('assoc_psoriasis_case.regenie' if domain=='autosomes' else model+'_psoriasis_case.regenie'),m['results_sha256'],n,m['estimated'],m['failed_or_unavailable']
def load(group,domain,model='male02'):
 parts=[];seen=set();checks=[];unavailable=0
 for p,h,n,e,u in sources(group,domain,model):
  assert sha(p)==h
  d=pd.read_csv(p,sep=r'\s+',dtype={'ID':str,'ALLELE0':str,'ALLELE1':str})
  assert len(d)==n and d.ID.is_unique and d.TEST.eq('ADD').all() and not seen.intersection(d.ID)
  seen.update(d.ID);v=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce')
  ok=np.isfinite(v).all(axis=1)&v.SE.gt(0)&v.LOG10P.ge(0)
  assert int(ok.sum())==e and int((~ok).sum())==u
  d=d.loc[ok].copy();d['CHROM']=pd.to_numeric(d.CHROM.replace({'X':23,'chrX':23}),errors='raise').astype('int16')
  parts.append(d[['CHROM','GENPOS','ID','ALLELE0','ALLELE1','N','N_CASES','N_CONTROLS','BETA','SE','LOG10P']]);unavailable+=u
  checks.append({'path':str(p),'sha256':h,'tested':n,'finite':e,'unavailable':u})
 d=pd.concat(parts,ignore_index=True);assert d.ID.is_unique
 key=group+'_'+domain+'_'+model;save(O/(key+'_source_receipts.private.json'),checks)
 RECEIPTS[key]={'n_finite':len(d),'unavailable':unavailable,'source_files':len(checks)}
 return d
def export_figure(fig,name):
 fig.savefig(O/(name+'.pdf'),dpi=600)
 fig.savefig(O/(name+'.svg'))
 fig.savefig(O/(name+'_600dpi.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 fig.savefig(O/(name+'_preview.png'),dpi=150)
 plt.close(fig)
def qq(ax,d,title):
 lam,q=quantiles(d.LOG10P.to_numpy())
 ax.fill_between(q.expected_log10P,q.null_95_low,q.null_95_high,color='#dce3e6',linewidth=0)
 ax.scatter(q.expected_log10P,q.observed_log10P,s=2,color='#006b75',rasterized=True,linewidths=0)
 xmax=float(q.expected_log10P.max())*1.04;ymax=max(xmax,float(q.observed_log10P.max())*1.04)
 ax.plot([0,xmax],[0,xmax],color='#87959a',lw=.7);ax.set(xlim=(0,xmax),ylim=(0,ymax),xlabel='Expected -log10(P)',ylabel='Observed -log10(P)',title=title)
 ax.text(.03,.96,'N = '+format(len(d),',')+'\nLambda GC = '+format(lam,'.3f'),transform=ax.transAxes,va='top',fontsize=8)
 ax.spines[['top','right']].set_visible(False)
 return lam
def manhattan(ax,d,title):
 ticks=[];labels=[];offset=0
 for i,(chrom,z) in enumerate(d.groupby('CHROM',sort=True)):
  pos=z.GENPOS.to_numpy();x=(pos+offset)/1e6
  ax.scatter(x,z.LOG10P,s=.6,color=['#416782','#5fa7a0'][i%2],linewidths=0,rasterized=True)
  ticks.append((offset+(pos.min()+pos.max())/2)/1e6);labels.append('X' if chrom==23 else str(chrom));offset+=int(pos.max())+10000000
 ax.axhline(-np.log10(5e-8),color='#a95142',lw=.8,ls='--');ax.axhline(-np.log10(5e-9),color='#7f8390',lw=.6,ls=':')
 ax.set(xticks=ticks,xticklabels=labels,xlabel='Chromosome' if len(labels)>1 else 'Chromosome X',ylabel='-log10(P)',title=title,ylim=(0,max(10,float(d.LOG10P.max())*1.04)))
 ax.tick_params(axis='x',labelsize=7);ax.spines[['top','right']].set_visible(False)
def calibration(group,domain,model,d):
 lam,q=quantiles(d.LOG10P.to_numpy());q.to_csv(O/(group+'_'+domain+'_'+model+'_QQ_ranks.private.tsv.gz'),sep='\t',index=False)
 r={'group':group,'domain':domain,'model':model,'valid_tests':len(d),'lambda_GC':lam,**{'P_lt_'+str(t):int(d.LOG10P.gt(-np.log10(t)).sum()) for t in [5e-8,3e-8,5e-9,1e-9]}}
 REPORTS.append(r);return r

def autosomes(group):
 done=O/(group+'_autosomal_assets_COMPLETE.json')
 if done.exists():return
 d=load(group,'autosomes');r=calibration(group,'autosomes','autosomal',d)
 no=d.loc[~(d.CHROM.eq(6)&d.GENPOS.between(25000000,34000000))];r2=calibration(group,'autosomes_MHC_excluded','autosomal',no)
 fig=plt.figure(figsize=(7.2,7.1));gs=fig.add_gridspec(2,2,height_ratios=[1,1.15]);manhattan(fig.add_subplot(gs[0,:]),d,'A  '+group.capitalize()+' autosomal association')
 qq(fig.add_subplot(gs[1,0]),d,'B  All autosomal tests');qq(fig.add_subplot(gs[1,1]),no,'C  MHC excluded')
 fig.subplots_adjust(left=.1,right=.98,bottom=.09,top=.94,hspace=.42,wspace=.35);export_figure(fig,'Autosomal_'+group)
 sig=d.loc[d.LOG10P.gt(-np.log10(5e-8))].copy();sig['OR']=np.exp(sig.BETA);sig['CI95_Wald_low']=np.exp(sig.BETA-1.95996398454*sig.SE);sig['CI95_Wald_high']=np.exp(sig.BETA+1.95996398454*sig.SE)
 sig.sort_values('LOG10P',ascending=False).to_csv(O/(group+'_autosomal_significant.private.tsv.gz'),sep='\t',index=False)
 selected=[];pending=sig.sort_values('LOG10P',ascending=False).copy()
 while len(pending):
  lead=pending.iloc[0];mhc=int(lead.CHROM)==6 and 25000000<=lead.GENPOS<=34000000
  take=pending.CHROM.eq(lead.CHROM)&(pending.GENPOS.between(25000000,34000000) if mhc else pending.GENPOS.between(lead.GENPOS-100000,lead.GENPOS+100000))
  row=lead.to_dict();row.update({'reporting_group_variants':int(take.sum()),'group_definition':'MHC chr6:25-34 Mb' if mhc else 'lead-centered +/-100 kb; not LD independence'});selected.append(row);pending=pending.loc[~take]
 pd.DataFrame(selected).to_csv(O/(group+'_autosomal_reporting_groups.private.tsv'),sep='\t',index=False)
 save(done,{'status':'complete','full':r,'MHC_excluded':r2,'source_scope':'historical fixed autosomal panel re-QC in strict cohort; not maximal source-wide discovery'})
 print('AUTOSOMAL_ASSETS_COMPLETE',group,flush=True)
def chromosome_x(group):
 models=['male02'] if group=='female' else ['male02','male01'];fig,axes=plt.subplots(len(models),2,figsize=(7.2,3.5*len(models)),squeeze=False)
 for row,model in enumerate(models):
  d=load(group,'X',model);r=calibration(group,'X',model,d)
  manhattan(axes[row,0],d,group.capitalize()+' X: male 0/'+model[-1]);qq(axes[row,1],d,'Full chromosome X calibration')
  for rs,gene,pos,ref,alt in LEADS:
   z=d.loc[d.ID.eq('chrX:'+str(pos)+':'+ref+':'+alt)]
   rec={'group':group,'model':model,'rsID':rs,'region':gene,'GRCh38_position':pos,'REF':ref,'ALT_effect_allele':alt,'status':'not in subgroup QC-passing test set' if z.empty else 'estimated'}
   if len(z):
    assert len(z)==1;v=z.iloc[0];rec.update({k:v[k].item() if isinstance(v[k],np.generic) else v[k] for k in ['N','N_CASES','N_CONTROLS','BETA','SE','LOG10P']});rec.update({'OR':np.exp(v.BETA),'CI95_Wald_low':np.exp(v.BETA-1.95996398454*v.SE),'CI95_Wald_high':np.exp(v.BETA+1.95996398454*v.SE)})
   LEADROWS.append(rec)
  if group=='pooled':REGIONS[model]=d.loc[d.GENPOS.between(5498724,6498724)|d.GENPOS.between(22417619,23417619)|d.GENPOS.between(67236589,68342269)].copy()
 fig.subplots_adjust(left=.1,right=.98,bottom=.09,top=.95,hspace=.4,wspace=.35);export_figure(fig,'Chromosome_X_'+group)
 print('X_ASSETS_COMPLETE',group,flush=True)
def regional():
 fig,axs=plt.subplots(3,2,figsize=(7.2,8.2));regions=[('NLGN4X',5498724,6498724),('PTCHD1-AS',22417619,23417619),('Xq12',67236589,68342269)]
 for row,(name,lo,hi) in enumerate(regions):
  ymax=max(10,*[float(REGIONS[m].loc[REGIONS[m].GENPOS.between(lo,hi),'LOG10P'].max())*1.04 for m in ['male02','male01']])
  for col,m in enumerate(['male02','male01']):
   ax=axs[row,col];z=REGIONS[m].loc[REGIONS[m].GENPOS.between(lo,hi)]
   ax.scatter(z.GENPOS/1e6,z.LOG10P,s=5,color='#3b8291',rasterized=True,linewidths=0)
   for rs,gene,pos,ref,alt in LEADS:
    if lo<=pos<=hi:
     ax.axvline(pos/1e6,color='#6e7074',ls=':',lw=.7);lead=z.loc[z.ID.eq('chrX:'+str(pos)+':'+ref+':'+alt)]
     if len(lead):ax.scatter([pos/1e6],lead.LOG10P,s=30,marker='D',color='#b45443',zorder=3)
   ax.axhline(-np.log10(5e-8),color='#a95142',lw=.7,ls='--');ax.set(xlim=(lo/1e6,hi/1e6),ylim=(0,ymax),xlabel='GRCh38 chromosome X position (Mb)',ylabel='-log10(P)',title=name+' | male 0/'+m[-1]);ax.spines[['top','right']].set_visible(False)
   if name!='NLGN4X':ax.text(.03,.96,'Original lead(s) excluded by HWE',transform=ax.transAxes,fontsize=7,va='top')
 fig.subplots_adjust(left=.1,right=.98,bottom=.07,top=.95,hspace=.42,wspace=.3);export_figure(fig,'Chromosome_X_complete_regional')
 for m,d in REGIONS.items():d.to_csv(O/('X_regions_'+m+'.private.tsv.gz'),sep='\t',index=False)
def cohort_tables():
 a=N/'cohort_audit';assert json.loads((a/'COMPLETE.json').read_text())['status']=='complete'
 flow=pd.read_csv(a/'sample_flow.private.tsv',sep='\t').set_index('stage');out=flow.loc[['original_analysis_cohort','after_canonical_sex_discordance','after_provider_QC_flag']].copy();out.index=['Original analysis cohort','After coverage, canonical ploidy and sex concordance','Strict provider-QC cohort'];assert not out.apply(lambda x:x.between(1,20)).any().any()
 out.to_csv(O/'Cohort_flow_aggregate_for_review.tsv',sep='\t')
 ar=pd.read_csv(a/'ancestry_group_counts.private.tsv',sep='\t');ar=ar.loc[ar.population.eq('strict_cohort_complete_file')].copy();small=ar[['N','cases','controls']].apply(lambda x:x.between(1,20)).any(axis=1)
 suppress=set(ar.index[small]);candidates=ar.loc[~small].sort_values('cases').index.tolist()
 if len(suppress)==1 and candidates:suppress.add(candidates[0])
 for idx in suppress:ar.loc[idx,['N','cases','controls']]='Suppressed'
 ar[['group','N','cases','controls']].to_csv(O/'Strict_cohort_all_ancestry_groups_for_review.tsv',sep='\t',index=False)
 save(O/'cohort_disclosure_review.private.json',{'source_sha256':{p.name:sha(p) for p in [a/'sample_flow.private.tsv',a/'ancestry_group_counts.private.tsv']},'rule':'Participant counts 1-20 suppressed; an additional entire ancestry row is suppressed when needed to avoid a unique complementary deduction. Sequential detailed exclusions retained privately and combined for display. Review alongside all manuscript tables before release.','release':'HELD; not program approval'})
def finish_manifest(complete):
 pd.DataFrame(REPORTS).to_csv(O/'Calibration_and_thresholds_for_review.tsv',sep='\t',index=False)
 pd.DataFrame(LEADROWS).to_csv(O/'X_lead_estimates_for_review.tsv',sep='\t',index=False)
 files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p),'release_status':'HELD: approved export process and disclosure review pending'} for p in sorted(O.iterdir()) if p.is_file() and p.name not in ['MANIFEST.json','COMPLETE.json','FAILED.json']]
 save(O/'MANIFEST.json',{'status':'complete' if complete else 'pooled_and_X_complete_waiting_autosomal_strata','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files,'code_sha256':sha(Path(__file__)),'versions':{p:im.version(p) for p in ['numpy','pandas','scipy','matplotlib']},'release_status':'HELD; no download attempted'})
def main():
 assert not (O/'COMPLETE.json').exists()
 (O/'README_FOR_REVIEW.md').write_text('Secure JID replacement assets. All files remain within Workbench; none is approved for export.\n\nAutosomal Manhattan plots contain every finite QC-passing test in the retained fixed panel. X plots contain every finite QC-passing source-X test. QQ statistics use all tests at unit weight, ranked before deterministic display thinning; bands are pointwise 95% beta order-statistic envelopes. Dotted threshold is 5e-9; dashed threshold is 5e-8. MHC exclusion is chr6:25-34 Mb, GRCh38. No LD-based independent-locus claim is made.\n\nRegional X plots include every corrected QC-passing test in each window, with the exact original lead positions marked. Excluded original leads have no primary estimate plotted. The plots do not encode LD, recombination or causal genes.\n\nX effects are per coded dosage unit: female 0/1/2, non-PAR male 0/2 primary and 0/1 sensitivity, PAR diploid. Wald CIs use the same BETA/SE; REGENIE LOG10P is retained. Group N and case/control N are variant-specific. These overlapping subgroup analyses are diagnostics, not replication.\n\nParticipant-level data and genotype QC files are excluded from figures but source paths/receipts and all result tables require release review. Cross-document disclosure review, author approval and final manuscript/S1/S2 reconciliation remain separate.\n')
 cohort_tables();autosomes('pooled')
 for g in GROUPS:chromosome_x(g)
 regional();finish_manifest(False)
 save(O/'POOLED_AND_X_COMPLETE.json',{'status':'complete','scope':'pooled autosomal, all X subgroup, regional plots and cohort accounting; autosomal subgroup figures pending'})
 while not (N/'results/autosomes_stratified/COMPLETE.json').exists():
  fail=N/'status/autosomal_stratified_FAILED.json'
  if fail.exists():raise RuntimeError('Subgroup run failure requires review')
  time.sleep(30)
 for g in GROUPS[1:]:autosomes(g)
 finish_manifest(True);save(O/'COMPLETE.json',{'status':'complete','code_sha256':sha(Path(__file__)),'release_status':'HELD pending approved export and disclosure review','figures_need_visual_review':True})
 print('ALL_SUBMISSION_ASSETS_BUILT_INSIDE_WORKBENCH',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
