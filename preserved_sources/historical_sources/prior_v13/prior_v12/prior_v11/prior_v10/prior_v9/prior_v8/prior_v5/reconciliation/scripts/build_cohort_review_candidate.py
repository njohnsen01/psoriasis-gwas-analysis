from pathlib import Path
import json,hashlib,datetime,textwrap,importlib.metadata
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';C=I/'cohort';O=I/'cohort_review_candidate';O.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def tab(n):return pd.read_csv(C/(n+'.private.tsv'),sep=chr(9))
def out(n,d):d.to_csv(O/(n+'.review.tsv'),sep=chr(9),index=False)
def fmt(n):return f'{int(n):,}'
def md(df):
 def cell(x):return '' if pd.isna(x) else str(x).replace('|','/').replace(chr(10),' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(cell(x) for x in row)+' |' for row in df.itertuples(index=False,name=None))
def main():
 rec=json.loads((C/'RECONCILED.private.json').read_text());assert all(rec['checks'].values())
 f=tab('full_sequential_flow');a=tab('ancestry_all_groups');gap=tab('historical_ancestry_gap');st=tab('five_followup_groups');age=tab('strict_age_summary').set_index('phenotype');ehr=tab('strict_EHR_summary').set_index('phenotype');life=tab('phenotype_record_and_date_lifecycle');obs=tab('strict_observation_criteria')
 # Suppress all intermediate QC phenotype counts; preserve the total and reason-specific N.
 q=f.copy();q['cases_retained']=q.cases_retained.astype(object);q['controls_retained']=q.controls_retained.astype(object)
 interm=q.stage.isin(['provider_coverage_failure','noncanonical_or_unavailable_ploidy','canonical_sex_discordance'])
 q.loc[interm,['cases_retained','controls_retained']]='suppressed; combined QC count below'
 q.loc[q.stage.eq('provider_QC_flag'),'stage']='strict_analysis_after_all_genomic_QC'
 q.loc[q.cases_retained.isna(),['cases_retained','controls_retained']]='not assigned at this stage'
 q=pd.concat([q,pd.DataFrame([{'stage':'combined_genomic_QC_exclusions','N_retained':pd.NA,'N_excluded_at_stage':rec['original']['N']-rec['strict']['N'],'cases_retained':pd.NA,'controls_retained':pd.NA,'definition':'Combined excluded participants: '+fmt(rec['original']['cases']-rec['strict']['cases'])+' cases; '+fmt(rec['original']['controls']-rec['strict']['controls'])+' controls. Do not expose sequential case/control differences.'}])],ignore_index=True)
 out('cohort_flow',q)
 # Individual N for all ancestry labels; phenotype counts for EAS/MID/SAS are combined.
 sa=a.loc[a.population.eq('strict_complete')].copy(); assert len(sa)==6
 ar=[]
 for _,r in sa.iterrows():
  smallgroup=r.group in ['EAS','MID','SAS']; ar.append({'group':r.group,'N':int(r.N),'cases':'combined below' if smallgroup else int(r.cases),'controls':'combined below' if smallgroup else int(r.controls),'followup_status':'separate ancestry model' if r.group in ['EUR','AFR','AMR'] else 'retained in pooled scan; below 200-case followup gate'})
 other=sa.loc[sa.group.isin(['EAS','MID','SAS'])];assert other.cases.lt(200).all()
 ar.append({'group':'EAS + MID + SAS (subtotal; not an additional group)','N':int(other.N.sum()),'cases':int(other.cases.sum()),'controls':int(other.controls.sum()),'followup_status':'combined for disclosure; no combined ancestry model'})
 ar.append({'group':'All strict-cohort groups','N':rec['strict']['N'],'cases':rec['strict']['cases'],'controls':rec['strict']['controls'],'followup_status':'pooled analysis'})
 ar=pd.DataFrame(ar);out('ancestry_groups',ar);out('historical_ancestry_gap',gap);out('five_followup_groups',st)
 clinical=[]
 for label in ['case','control']:
  z=age.loc[label];e=ehr.loc[label]
  clinical.append({'phenotype':label,'cohort_N':int(z.N_total),'age_median_years':round(z['median'],1),'age_Q1':round(z.q1,1),'age_Q3':round(z.q3,1),'age_mean':round(z['mean'],1),'age_SD':round(z.sd,1),'age_available_N':int(z.N_available),'EHR_span_available_N':int(e.N_available),'EHR_span_missing_N':int(e.N_missing),'EHR_span_median_years':round(e['median'],1),'EHR_span_Q1':round(e.q1,1),'EHR_span_Q3':round(e.q3,1),'EHR_span_mean':round(e['mean'],1),'EHR_span_SD':round(e.sd,1)})
 clinical=pd.DataFrame(clinical);out('age_and_EHR',clinical);out('phenotype_lifecycle',life);out('observation_criteria',obs)
 # Direct small-cell screen is necessary but does not certify cross-file disclosure safety.
 screens=[]
 for filename,cols in [('cohort_flow',['N_retained','N_excluded_at_stage','cases_retained','controls_retained']),('ancestry_groups',['N','cases','controls']),('historical_ancestry_gap',['N','cases','controls']),('five_followup_groups',['N','cases','controls']),('age_and_EHR',['cohort_N','age_available_N','EHR_span_available_N','EHR_span_missing_N']),('phenotype_lifecycle',['N','cases','controls']),('observation_criteria',['N_meets','N_does_not_meet','N_total'])]:
  data=pd.read_csv(O/(filename+'.review.tsv'),sep=chr(9));bad=[]
  for col in cols:
   v=pd.to_numeric(data[col],errors='coerce');bad.extend([(int(i),col) for i in data.index[v.between(1,20)]])
  assert not bad,(filename,bad);screens.append({'file':filename,'direct_small_cells':0})
 def stat(z):return f"{z['median']:.1f} years (IQR {z.q1:.1f}–{z.q3:.1f})"
 method=("The cohort was reconstructed from the retained Cohort Builder SQL and phenotype-processing code in CDR C2024Q3R8. Initial eligibility required WGS availability, age_at_cdr 18–125 years, no recorded death, and one of the two source-specified recorded-sex concepts. The original SQL selected the union of a case arm and a control arm. The case arm required at least two distinct (participant, entry_date, concept) records for a qualifying concept in the source event index. The control arm excluded the source-specified condition and systemic-drug concept sets; these drug exclusions were not applied to the case arm. Final case status required at least two distinct calendar dates of qualifying condition records. The retained source code labels all other selected participants as controls; the executed cohort audit found no participants with exactly one qualifying record or exactly one qualifying date in either the original or corrected analysis cohort. An inner join to the retained complete provider PC source removed no clinically selected participants. Corrected sample QC sequentially removed provider coverage failures, noncanonical or unavailable XX/XY ploidy, canonical sex discordance, and remaining provider sample-QC flags. These sets were verified to be disjoint and to reconcile to the frozen corrected phenotype and covariate inputs.")
 ancestry=("The earlier EUR/AFR/AMR follow-up groups used a smaller retained ancestry file. Participants absent from that file and participants in other groups remained in the pooled GWAS. Corrected groups use the exact uppercase ancestry_pred label in echo_v4_r2.ancestry_preds.tsv, with no additional assignment-probability threshold. All six provider labels are accounted for. EUR, AFR and AMR were analyzed separately; EAS, MID and SAS each fell below the code's minimum of 200 cases (500 controls were also required) and remain in the pooled analysis. Meeting this sample-count gate does not establish power, and variant-specific genotype/allele-count criteria may further reduce analyzed N. Female/male strata use the frozen sex covariate concordant with provider XX/XY ploidy. These are overlapping within-cohort analyses, not independent replication.")
 observation=("Age covariates were calculated from the retained 2023-10-01 UTC reference date and date of birth, divided by 365.25; this reference date is not asserted to be the verified clinical data cutoff. Eligibility used provider age_at_cdr separately. EHR span was the interval from the first to last recorded condition, procedure, drug or eligible clinical-visit date, divided by 365.25, with missing values when the EHR-availability flag was absent. This measures recorded observation span, not continuous follow-up. The primary phenotype did not require a minimum EHR observation duration. An EHR flag plus at least two distinct clinical visit dates separated by at least 365 days was tabulated as a diagnostic criterion; it has not been substituted for the frozen primary controls.")
 results=(f"The reconstructed original cohort contained {fmt(rec['original']['N'])} participants ({fmt(rec['original']['cases'])} cases and {fmt(rec['original']['controls'])} controls). Corrected sample QC retained {fmt(rec['strict']['N'])} participants ({fmt(rec['strict']['cases'])} cases and {fmt(rec['strict']['controls'])} controls). The historical ancestry-follow-up gap comprised {fmt(gap.loc[gap.group.eq('absent_from_legacy_file'),'N'].iloc[0])} participants absent from the legacy assignment file and {fmt(gap.loc[gap.group.eq('other_legacy_groups'),'N'].iloc[0])} assigned to other groups; these participants were not excluded from the original pooled scan. In the corrected cohort, age was {stat(age.loc['case'])} in cases and {stat(age.loc['control'])} in controls. EHR span was available for {fmt(ehr.loc['case','N_available'])} cases and {fmt(ehr.loc['control','N_available'])} controls, with medians of {stat(ehr.loc['case'])} and {stat(ehr.loc['control'])}, respectively. EHR span was unavailable for {fmt(ehr.loc['control','N_missing'])} controls. These differences motivate the outstanding EHR-observed control sensitivity analysis.")
 parts=['# Cohort reconciliation: manuscript integration candidate','HELD INSIDE WORKBENCH. Not approved for export. Not yet inserted in Word/S1/S2.','## Methods candidate',method,ancestry,observation,'## Results candidate',results,'## Full cohort flow',md(q),'## Corrected ancestry groups',md(ar),'## Historical ancestry gap',md(gap),'## Corrected sex and ancestry follow-up models',md(st),'## Age and EHR summaries',md(clinical),'## Phenotype record/date lifecycle',md(life),'## Observation diagnostics',md(obs),'## Disclosure and integration requirements','Exact private tables retain all counts. The candidate suppresses intermediate QC phenotype counts and combines EAS/MID/SAS case/control counts. Do not combine this candidate with historical tables that allow suppressed counts to be reconstructed. A joint review of Word, S1, S2, plots and all supplementary files remains required. Passing the direct small-cell screen is not release approval. Figure cohort N must remain distinguished from variant-specific analyzed N.','## Scientific limits','This audit verifies reconstruction and reporting; it does not provide chart validation, restore the historical X discovery claims, establish adequate observation for all controls, or complete the relatedness-aware alternative-control analysis.']
 (O/'MANUSCRIPT_COHORT_INTEGRATION_CANDIDATE.md').write_text('\n\n'.join(parts).replace('\\n','\n')+'\n')
 # Vector cohort flow; exact N and QC reason counts, no small phenotype subcells.
 by=f.set_index('stage');plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
 fig,ax=plt.subplots(figsize=(8.2,10.0));fig.subplots_adjust(left=0,right=1,top=1,bottom=0);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off');texts=[];boxes=[]
 ax.text(.04,.97,'Psoriasis GWAS cohort reconstruction',fontsize=16,fontweight='bold',va='top',color='#17384A')
 ax.text(.04,.935,'Retained source definition and corrected sample QC',fontsize=10,color='#455A64')
 def box(x,y,w,h,lines,primary=False):
  p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.007,rounding_size=0.008',facecolor='#EAF2F7' if primary else '#F4F4F4',edgecolor='#67879A' if primary else '#B2B7BA',linewidth=.8);ax.add_patch(p)
  t=ax.text(x+w/2,y+h/2,lines,ha='center',va='center',fontsize=10 if primary else 9,linespacing=1.45,color='#172F3E');texts.append(t);boxes.append((t,p));return p
 def arrow(x1,y1,x2,y2):ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'-|>','color':'#677A86','lw':1.1,'shrinkA':2,'shrinkB':2})
 w=.46;x=.04;mid=x+w/2;ys=[.815,.655,.505,.345,.185];hs=[.08,.08,.08,.10,.10]
 labels=[f"WGS availability flag present\nN = {fmt(by.loc['WGS_available','N_retained'])}",f"Demographic eligibility met\nN = {fmt(by.loc['eligible_recorded_sex','N_retained'])}",f"After control condition exclusions\nN = {fmt(by.loc['clinical_condition_exclusion','N_retained'])}",f"Original analysis cohort\nN = {fmt(rec['original']['N'])}\n{fmt(rec['original']['cases'])} cases; {fmt(rec['original']['controls'])} controls",f"Corrected primary cohort\nN = {fmt(rec['strict']['N'])}\n{fmt(rec['strict']['cases'])} cases; {fmt(rec['strict']['controls'])} controls"]
 for y,h,t in zip(ys,hs,labels):box(x,y,w,h,t,True)
 branches=[f"Excluded sequentially\nRecorded death: {fmt(by.loc['no_recorded_death','N_excluded_at_stage'])}\nOther/unavailable recorded sex: {fmt(by.loc['eligible_recorded_sex','N_excluded_at_stage'])}",f"Control condition exclusions\nN = {fmt(by.loc['clinical_condition_exclusion','N_excluded_at_stage'])}",f"Control drug exclusions\nN = {fmt(by.loc['control_drug_exclusion','N_excluded_at_stage'])}",f"Genomic QC exclusions: {fmt(rec['original']['N']-rec['strict']['N'])}\nCoverage failure: {fmt(by.loc['provider_coverage_failure','N_excluded_at_stage'])}\nOther/unavailable ploidy: {fmt(by.loc['noncanonical_or_unavailable_ploidy','N_excluded_at_stage'])}\nCanonical sex discordance: {fmt(by.loc['canonical_sex_discordance','N_excluded_at_stage'])}\nRemaining provider QC flags: {fmt(by.loc['provider_QC_flag','N_excluded_at_stage'])}"]
 for j in range(4):
  lo=ys[j+1]+hs[j+1];hi=ys[j];ym=(lo+hi)/2;arrow(mid,hi,mid,lo);height=[.092,.074,.074,.138][j];box(.57,ym-height/2,.38,height,branches[j]);arrow(mid+.015,ym,.557,ym)
 notes=['Eligibility: age_at_cdr 18–125 years, no recorded death, and source-specified recorded sex.','The WGS subset had no exclusions for age or missing person records; the PC-source join removed none.','Case and control-arm rules precede final phenotype coding. Drug exclusions apply to the control arm.','Cohort N is before variant-specific QC and genotype missingness. Related participants were retained.','Detailed group assignment, case/control accounting and age/EHR summaries accompany this flow.']
 for j,line in enumerate(notes):ax.text(.04,.135-j*.022,line,fontsize=8.1,ha='left',va='top',color='#405461')
 fig.canvas.draw();renderer=fig.canvas.get_renderer();fit=[]
 for t,p in boxes:
  tb=t.get_window_extent(renderer);pb=p.get_window_extent(renderer);ok=pb.contains(tb.x0,tb.y0) and pb.contains(tb.x1,tb.y1);assert ok,t.get_text();fit.append(ok)
 for ext in ['pdf','svg','png','tiff']:
  kwargs={'dpi':600 if ext=='tiff' else 160}
  if ext=='tiff':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(O/('cohort_flow.'+ext),**kwargs)
 plt.close(fig)
 with Image.open(O/'cohort_flow.tiff') as im:assert im.size==(4920,6000) and tuple(round(v) for v in im.info['dpi'])==(600,600)
 sources={str(p):sha(p) for p in C.glob('*.private.tsv')};sources[str(C/'RECONCILED.private.json')]=sha(C/'RECONCILED.private.json')
 report={'status':'cohort_review_candidate_built','created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_hashes':sources,'direct_count_screen':screens,'figure_text_boxes_fit':all(fit),'visual_QA':'pending','workbook_integration':'pending authoritative current S1 source; existing Workbench workbook is historical','release':'HELD; coordinated disclosure review and approved export required','script_sha256':sha(__file__),'versions':{p:importlib.metadata.version(p) for p in ['pandas','numpy','matplotlib','pillow']},'files':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in O.iterdir() if p.is_file()}}
 (O/'BUILD_MANIFEST.private.json').write_text(json.dumps(report,indent=2)+'\n');print('COHORT_REVIEW_CANDIDATE_BUILT',flush=True)
if __name__=='__main__':main()
