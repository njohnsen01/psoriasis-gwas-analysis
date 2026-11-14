from pathlib import Path
import ast,copy,datetime,hashlib,importlib.metadata,json,math,re,shutil,textwrap
import pandas as pd
import openpyxl
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter
I=Path(__file__).resolve().parents[1]; N=I.parent; B=N.parent
R=Path('/home/jupyter/psoriasis_reproducibility_20260910'); P=R/'pipeline_v3'
W=B/'phenotype_validation_20260910'; T=B/'SII_timing_sensitivity_20260910'
V=I/'supplement_staging_v5'; V.mkdir(exist_ok=True);assert not (V/'S1_reconciled_tables_HELD.xlsx').exists()
D=V/'methods_sources';D.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''):h.update(c)
 return h.hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,default=str)+chr(10))
def literals(p):
 out={}
 for n in ast.walk(ast.parse(p.read_text())):
  if isinstance(n,ast.Assign):
   for t in n.targets:
    if isinstance(t,ast.Name):
     try:out[t.id]=ast.literal_eval(n.value)
     except (ValueError,TypeError):pass
 return out
manifest=json.loads((P/'manifest.json').read_text())
assert all(sha(P/f)==h for f,h in manifest['files'].items())
source_paths=[P/'run.py',P/'manifest.json',P/'ancestry_source.json',P/'REPLAY_VALIDATION.json',P/'raw_laboratory/VALIDATION.json',P/'raw_laboratory/CLI_VALIDATION.json',P/'raw_laboratory/query.sql',P/'raw_laboratory/executed_reconstruct_raw_labs.py',R/'replay_cohort_v2/COMPLETE.json',R/'fresh_lab_processing_v3/COMPLETE.json']
source_paths += [P/'cells'/('cell_%03d.py'%n) for n in [3,4,6,11,12,13,14,15,16,95,97,98,99,100]]
source_paths += [I/'scripts'/n for n in ['prepare_cohort_reconciliation.py','reconcile_cohort_tables.py','EHR_final_cohort.py']]
source_paths += [W/n for n in ['condition_roots.tsv','expanded_case_concepts.tsv','phenotype_audit.sql','phenotype_sensitivity_protocol.json']]
source_paths += [T/n for n in ['analysis_protocol.json','lab_concept_map.tsv','unit_mapping.tsv','run_sii.py','verify_sii_and_timing.py','timing.sql']]
source_paths += [P/'retained/audit/SII_timing_sensitivity_20260910/run_sii.py',P/'retained/audit/calculate_cohort_and_ehr_summaries.py',B/'ehr_observation_span.sql',I/'cohort/RECONCILED.private.json',I/'EHR_final_cohort_sensitivity/COMPLETE.private.json']
source_paths += [R/'replay_cohort_v2/runtime_versions.json',R/'fresh_lab_processing_v3/runtime_versions.json']
source_paths += [D/n for n in ['lab_vocabulary_current.private.tsv','lab_vocabulary_query.sql','lab_vocabulary_query_receipt.private.json','lab_concept_usage.private.json']]
assert all(p.is_file() for p in source_paths)
replay=json.loads((P/'REPLAY_VALIDATION.json').read_text())
assert all(replay[k] for k in ['lab_ids_match','laboratory_values_match','psoriasis_pheno.tsv_IDs_match','psoriasis_pheno.tsv_values_match','psoriasis_covariates.tsv_IDs_match','psoriasis_covariates.tsv_values_match','SII_numbers_match'])
raw=json.loads((P/'raw_laboratory/VALIDATION.json').read_text());assert raw['status']=='validated' and raw['person_level_values_match'] and raw['classified_row_multiset_128bit_fingerprints_match']
assert sha(T/'run_sii.py')==sha(P/'retained/audit/SII_timing_sensitivity_20260910/run_sii.py')
source_rows=[]
for p in source_paths:
 source_rows.append({'file':p.name,'source':str(p),'sha256':sha(p),'lines':len(p.read_text().splitlines()),'role':'Executed source or retained execution/validation evidence; historical and corrected stages are distinguished in S2','release':'HELD'})
source_df=pd.DataFrame(source_rows)
# Extract code constants without executing the analysis cells.
c97=literals(P/'cells/cell_097.py');c98=literals(P/'cells/cell_098.py');sii=literals(T/'run_sii.py')
roots=pd.read_csv(W/'condition_roots.tsv',sep=chr(9),dtype=str,keep_default_na=False)
expanded=pd.read_csv(W/'expanded_case_concepts.tsv',sep=chr(9),dtype=str,keep_default_na=False)
case_ids={'140168','36715780','37205056','37205057','4031141','4031142','4031649','4063433','4223485'}
roots['record_type']='query_root';roots['case_query_root']=roots.concept_id.isin(case_ids);roots['control_condition_exclusion_root']=True
expanded['record_type']='expanded_case_concept';expanded['case_query_root']='not a root row';expanded['control_condition_exclusion_root']='case-set descendant'
expanded['vocabulary_id']='not supplied in retained expansion table';expanded['domain_id']='not supplied in retained expansion table';expanded['standard_concept']='S';expanded['invalid_reason']='not supplied in retained expansion table'
# Confirm expanded labels shared with root table have identical vocabulary codes/names.
m=expanded.merge(roots,on='concept_id',suffixes=('_e','_r'))
assert (m.concept_code_e==m.concept_code_r).all() and (m.concept_name_e==m.concept_name_r).all()
phen=pd.concat([roots,expanded],ignore_index=True).fillna('not applicable')
phen=phen[['record_type','concept_id','concept_name','vocabulary_id','concept_code','case_query_root','control_condition_exclusion_root','in_arthritis_set','standard_concept','invalid_reason']]
concept=pd.read_csv(T/'lab_concept_map.tsv',sep=chr(9),dtype=str,keep_default_na=False)
assert concept.concept_id.is_unique
mapping={}
for key,label in [('LDL_IDS','LDL'),('HDL_IDS','HDL'),('TG_IDS','TG'),('PLT_IDS','platelet_count'),('ANC_IDS','ANC'),('ALC_IDS','ALC')]:
 for c in c97[key]:assert c not in mapping;mapping[c]=label
assert set(mapping)-set(concept.concept_id)=={'1497147'}
current=pd.read_csv(D/'lab_vocabulary_current.private.tsv',sep=chr(9),dtype=str,keep_default_na=False)
assert set(current.concept_id)==set(concept.concept_id)
for col in ['concept_name','concept_code','vocabulary_id','standard_concept','invalid_reason']:
 assert concept.set_index('concept_id')[col].sort_index().equals(current.set_index('concept_id')[col].sort_index())
usage=json.loads((D/'lab_concept_usage.private.json').read_text());assert usage['unmatched_vocabulary_concept_rows']==0
concept['annotation_status']='CDR vocabulary returned; original annotation agrees'
concept=pd.concat([concept,pd.DataFrame([{'concept_id':'1497147','concept_name':'Not found in queried CDR v8 vocabulary','concept_code':'unavailable','vocabulary_id':'unavailable','standard_concept':'unavailable','invalid_reason':'unavailable','annotation_status':'Classifier constant only; no retained classified rows'}])],ignore_index=True)
concept['observed_in_classified_data']=concept.concept_id.isin(usage['observed_concepts'])
concept.insert(0,'biomarker',concept.concept_id.map(mapping))
concept['standardized_unit']=concept.biomarker.map(lambda x:'mg/dL' if x in ['LDL','HDL','TG'] else '10^3 cells/uL')
concept['inclusive_min']=concept.biomarker.map(lambda x:c98['filters'][x][0]);concept['inclusive_max']=concept.biomarker.map(lambda x:c98['filters'][x][1])
concept['unit_rule']=concept.biomarker.map(lambda x:'If unit text contains mmol: multiply by '+('88.57' if x=='TG' else '38.67')+'; otherwise original numeric value retained. No strict unknown-unit rejection.' if x in ['LDL','HDL','TG'] else 'Original substring-based absolute-count selection/conversion; strict exact-unit mapping sensitivity in SII_unit_mapping.')
unit=[]
for dictionary in ['standard','source']:
 for label,factor in sii[dictionary].items():unit.append({'record_type':dictionary+'_dictionary','unit_name':label if dictionary=='standard' else '', 'unit_source_value':label if dictionary=='source' else '', 'factor_to_thousand_cells_per_uL':factor,'handling':'Require a recognized standard-unit name. If the source-unit label also maps, its factor must agree. An unrecognized source label is allowed when the standard name is recognized.'})
u=pd.read_csv(T/'unit_mapping.tsv',sep=chr(9),keep_default_na=False)
for r in u.to_dict('records'):unit.append({'record_type':'observed_pair','unit_name':r['unit_name'],'unit_source_value':r['unit_source_value'],'factor_to_thousand_cells_per_uL':float(r['factor']) if r['factor']!='' else None,'handling':'Recorded pair from executed strict mapping audit; blank factor is unavailable, not zero.'})
unit=pd.DataFrame(unit)
tables={'Phenotype_concepts':phen,'Lab_methods':concept,'SII_unit_mapping':unit,'Methods_provenance':source_df}
for name,df in tables.items():df.to_csv(D/(name+'.private.tsv'),sep=chr(9),index=False)
oldpath=I/'supplement_staging_v4/S1_reconciled_tables_HELD.xlsx';oldhash=sha(oldpath)
w=openpyxl.load_workbook(oldpath);oldnames=w.sheetnames.copy();snapshot={s:[[(c.value,copy.copy(c._style),c.number_format) for c in row] for row in w[s]] for s in oldnames}
expected={}
def val(x):
 if pd.isna(x):return None
 return x.item() if hasattr(x,'item') else x
for name,df in tables.items():
 assert name not in w;ws=w.create_sheet(name);ws.append(df.columns.tolist())
 for row in df.itertuples(index=False,name=None):ws.append([val(x) for x in row])
 expected[name]=list(ws.values)
 ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.showGridLines=False
 for j,c in enumerate(ws[1],1):
  c.font=Font(name='Arial',size=10,bold=True,color='FFFFFF');c.fill=PatternFill('solid',fgColor='234A61');c.alignment=Alignment(wrap_text=True,vertical='top')
  h=str(c.value);ws.column_dimensions[get_column_letter(j)].width=70 if h in ['source','unit_rule','role','handling'] else 53 if h=='concept_name' else 24
 ws.row_dimensions[1].height=42
 for row in ws.iter_rows(min_row=2):
  needed=1
  for c in row:
   c.font=Font(name='Arial',size=10);c.alignment=Alignment(wrap_text=True,vertical='top');c.fill=PatternFill('solid',fgColor='F2F6F8' if c.row%2==0 else 'FFFFFF')
   width=ws.column_dimensions[c.column_letter].width;needed=max(needed,math.ceil(len(str(c.value or ''))/max(12,width-3)))
  ws.row_dimensions[row[0].row].height=max(30,needed*14+5)
 ws.page_setup.orientation='landscape';ws.page_setup.paperSize=ws.PAPERSIZE_A3;ws.page_setup.fitToWidth=1;ws.page_setup.fitToHeight=0;ws.sheet_properties.pageSetUpPr.fitToPage=True;ws.print_title_rows='1:1';ws.oddFooter.center.text='HELD — secure Methods integration draft'
 sf=w['Source_manifest'];f=D/(name+'.private.tsv');sf.append([name,str(f),sha(f),len(df),len(df.columns),'Executed cohort/laboratory Methods support','HELD'])
 for c in sf[sf.max_row]:c._style=copy.copy(sf.cell(sf.max_row-1,c.column)._style)
 w['Source_manifest'].auto_filter.ref=w['Source_manifest'].dimensions
for label,detail in [('Methods integration v5','Exact phenotype concepts, laboratory concept/unit/plausibility rules, concurrency limitations and source hashes are added. Prior association values are preserved.'),('Remaining Methods limits','Historical lab unit heuristics and lack of diagnosis/treatment restriction are disclosed. No chart validation, treatment-naive phenotype, common-specimen confirmation or fresh-machine replay is established.')]:
 ws=w['README'];ws.append([label,detail]);
 for c in ws[ws.max_row]:c._style=copy.copy(ws.cell(ws.max_row-1,c.column)._style)
 ws.row_dimensions[ws.max_row].height=65
out=V/'S1_reconciled_tables_HELD.xlsx';w.save(out);z=openpyxl.load_workbook(out)
for s in oldnames:
 for i,row in enumerate(snapshot[s],1):
  for j,(value,style,fmt) in enumerate(row,1):
   c=z[s].cell(i,j);assert c.value==value and c._style==style and c.number_format==fmt,(s,i,j)
for s,rows in expected.items():assert [tuple(None if v=='' else v for v in r) for r in z[s].values]==[tuple(None if v=='' else v for v in r) for r in rows],s
assert sha(oldpath)==oldhash
cohort_text='''### Executed cohort selection and case assignment

The original cohort was constructed from All of Us Controlled Tier CDR v8 (C2024Q3R8). The generated selection query required a whole-genome-variant availability flag, provider age_at_cdr from 18 through 125 years, a person record, no recorded death, and one of the two source-specified sex-at-birth concepts. It combined a case-selection arm with a control-selection arm. The case arm required at least two distinct event dates for at least one qualifying standard concept in the descendant-expanded psoriasis set. The control arm excluded any event in the broader psoriasis/psoriatic-arthritis condition set and any descendant of the source-specified medication roots. Case-arm eligibility took precedence over those control exclusions. The complete executed SQL is retained with its checksum in Methods_provenance.

After selection, case labels were assigned from condition_occurrence using at least two distinct calendar dates of qualifying diagnoses, irrespective of whether the codes differed. Participants were inner-joined to the versioned provider principal-component source. The two selection/labeling steps are therefore documented separately. The executed lifecycle audit found no participant with exactly one qualifying record or exactly one qualifying calendar date in either the original analyzed cohort or the corrected cohort; such participants were not included as controls. The phenotype was broader than plaque psoriasis: the expanded case set includes psoriatic arthritis, pustular, guttate, flexural, scalp and nail psoriasis, among other subtypes. Psoriatic-arthritis-only qualifying records can satisfy the coded case definition. Phenotype_concepts provides the roots and all expanded concepts, with the arthritis-set overlap. No independent chart-adjudication validation has been established.

The corrected primary cohort applied the verified provider coverage-failure list, removed noncanonical or unavailable XX/XY ploidy assignments, removed canonical XX/XY assignments discordant with the retained analysis sex covariate, and removed the remaining provider-flagged samples. Original sex covariates were not silently recoded from updated provider ploidy. Cohort_flow reports the sequential exclusions and case/control denominators. Ancestry_groups retains all provider genetic-similarity assignments and distinguishes the historical incomplete ancestry file from the complete provider assignments. EUR, AFR and AMR follow-up sample eligibility required at least 200 cases and 500 controls; per-variant count/call-rate gates were additional requirements and do not demonstrate statistical power.

Age was calculated as the integer day difference between 2023-10-01 UTC and date of birth, divided by 365.25. That is the implemented analysis reference date; the code comment calling it a CDR cutoff does not establish an independent data-freeze date. Sex was coded male=1 and female=0 in the covariate file. The original covariates were sex, age, age squared, age-by-sex and PC1–PC10. Later model-specific covariate modifications are reported in their respective sections.

### Clinical observation and control sensitivity

The original primary control definition did not require has_ehr_data=1 or a minimum clinical-observation period. EHR observation time was calculated as (latest minus earliest recorded date)/365.25 using condition, procedure, drug-exposure and visit dates. The visit component excluded visit_type_concept_id 44818519. Measurement dates were not included in that span query. Span was marked unavailable where the provider EHR flag was not 1; an observed zero-day span was retained as zero. This records available clinical history and does not establish continuous enrollment, observation before diagnosis, or absence of psoriasis.

The post-hoc targeted EHR sensitivity retained the same cases and compared original controls, EHR-observed original controls, and EHR-observed controls with medication-only exclusions restored. EHR-observed controls required has_ehr_data=1, at least two distinct non-study visit_start_date values spanning at least 365 days, and no qualifying psoriasis/psoriatic-arthritis event. Restored controls remained subject to the condition exclusion and the same corrected provider sample-QC gates. Both all-participant and provider-pruned analyses were performed. The restored-control targeted models did not use a newly refitted genome-wide REGENIE predictor, and they do not establish genome-wide calibration or independent phenotype validation.
'''
lab_text='''### Executed laboratory extraction, units and participant summaries

Laboratory records came from the CDR v8 measurement table for the original phenotype/covariate cohort. The extraction selected non-null numeric measurements using the candidate concept list plus standard/source concept-name searches for lipids and absolute blood-cell counts. It retained measurement concept and source metadata, unit labels, numeric values and measurement_datetime. The subsequent classifier declared 21 concept IDs, all listed in Lab_methods. A direct CDR-vocabulary lookup returned the 20 annotated LOINC concepts and confirmed their prior labels/codes. The additional classifier constant, 1497147, was not found in the queried vocabulary and had no rows in the retained classified laboratory data. It is documented as an unmatched, unused classifier entry; no inferred concept identity was assigned and the original classifier was preserved. The broader extraction candidates were not themselves the final analyzed concept set.

LDL, HDL and triglycerides were selected by concept membership. Platelet, absolute-neutrophil and absolute-lymphocyte counts also required an absolute-count substring in the combined unit labels; percentage counts were not intended as eligible components. In the original normalization, lipid unit text containing mmol triggered multiplication by 38.67 for LDL/HDL or 88.57 for triglycerides. Otherwise the numeric value was retained. Blood-cell values recognized as cells per microliter were divided by 1,000 unless a thousands-scale label was also recognized. The original heuristic was not a complete explicit rejection rule for unknown or conflicting lipid units; the stricter sensitivity below applied to SII cell counts.

Values were converted to numeric, and missing or nonpositive values were excluded before aggregation. The inclusive standardized limits were LDL 1–500 mg/dL, HDL 1–250 mg/dL, triglycerides 1–2,000 mg/dL, platelets 1–2,000 x 10^3 cells/uL, and neutrophils/lymphocytes 0.01–100 x 10^3 cells/uL. These are the implemented plausibility limits, not clinical reference ranges. Within each participant and component, the median of retained measurements was calculated separately. The original processing did not deduplicate identical rows, require concurrent components, or restrict measurements to a period before diagnosis or treatment. Different components may therefore summarize different visits. Missing measurement times did not by themselves exclude a value from the original median.

The original SII variable was median(platelets) x median(neutrophils) / median(lymphocytes), requiring positive components and a finite positive derived result. It was not the median of a same-draw SII. Triglycerides and SII were natural-log transformed. LDL, HDL, log triglycerides and log SII were standardized using the available participant-level mean and sample standard deviation in the original biomarker construction. Outcome availability and subsequent model-specific participant selection remain distinct. The retained model-specific Student-t inference and multiplicity families in the association tables are unchanged.

### Same-timestamp SII sensitivity and clinical timing

The later sensitivity rebuilt the original SII and then matched components by participant and exact measurement_datetime, parsed in UTC. Null timestamps and timestamps equal to midnight were excluded; there was no same-day rounding or nearest-neighbor match. Exact duplicate normalized component records were collapsed. A timestamp was discarded if any component had conflicting distinct normalized values, and all three components were required. SII was calculated at each qualifying timestamp, required to be finite and positive, and then summarized by the participant median. Missing/zero lymphocytes were excluded; no zero replacement or pseudocount was used.

Two unit analyses were retained: original conversion rules and strict explicit mappings. Under strict mapping, the standard unit name had to match the explicit dictionary. If the source unit label also matched its dictionary, its factor had to agree with the standard-unit factor. An unrecognized source label was allowed when the standard unit name was recognized; an unknown standard unit name or conflicting known factors led to exclusion. SII_unit_mapping gives the executed dictionaries and observed-pair audit, with factors to 10^3 cells/uL. The same plausibility limits were applied. Comparisons of original and concurrent SII used identical complete-case participants within each contrast, natural-log transformation, and sample-standard-deviation scaling within that comparison population. All-participant and provider-pruned analyses, full-cohort and cases-only models, and both X dosage assumptions remained separately labeled. The 25-target within-outcome family and historical 100-test benchmark are distinct, post-hoc correction contexts.

Exact timestamp matching is a concurrency proxy: no specimen identifier or specimen link was available to verify a common blood draw. Midnight exclusion also removes any genuine midnight draws. The timing audit compared qualifying timestamps with the earliest observed case-set diagnosis and earliest exposure in the original medication set. It did not restrict the association analyses to pretreatment values. The retained timing summary concerned the provider-pruned cases represented in that earlier biomarker audit and is not silently substituted for a new corrected-cohort summary. Recorded first diagnosis does not establish disease onset, and recorded medication exposure does not establish a psoriasis indication or treatment-naive status.
'''
old_s2=I/'supplement_staging_v4/S2_reconciliation_integration.private.md';oldtext=old_s2.read_text();newtext=oldtext
anchor='\n## Primary estimate provenance and dosage units';assert anchor in newtext
newtext=newtext.replace(anchor,chr(10)+cohort_text+anchor,1)
anchor2='\n## Functional evidence, novelty and replication';assert anchor2 in newtext
newtext=newtext.replace(anchor2,chr(10)+lab_text+anchor2,1)
oldinstruction='Read their protocol and actual covariate/phenotype source before wording the main Methods.'
assert oldinstruction in newtext;newtext=newtext.replace(oldinstruction,'The executed cohort and sensitivity definitions are detailed below and linked to source hashes.',1)
oldinstruction2='Laboratory concept/unit mappings, plausible-value limits, timestamp rules and zero-count handling must be carried from the executed processing code and original validation record.'
assert oldinstruction2 in newtext;newtext=newtext.replace(oldinstruction2,'The following subsections and S1 Methods tables document the executed concept/unit mappings, plausibility limits, timestamp rules and zero-count handling.',1)
newtext += chr(10)+'## Cohort and laboratory source audit (secure v5)'+chr(10)*2+'The original pipeline source-manifest hashes were rechecked. Retained completed cohort and fresh-laboratory reconstruction receipts report matching participants, phenotype/covariate values, classified-row fingerprints and participant-level laboratory summaries. These checks concern the documented source stages; they are not a new GWAS or a fresh-machine end-to-end replay. Methods_provenance records the reviewed files and their current hashes, including execution versions. The v4 association tables and numerical S2 tables are preserved. This v5 supplement remains a secure integration draft pending complete scope review, disclosure review, approved release, native-file checks and author approval.'+chr(10)
(V/'S2_reconciliation_integration.private.md').write_text(newtext)
# Independently preserve every old numerical markdown-table line.
oldtables=[x for x in oldtext.splitlines() if x.startswith('|')];newtables=[x for x in newtext.splitlines() if x.startswith('|')];assert oldtables==newtables
# Produce compact manuscript insertion text from saved workbook rows, retaining exact private denominators in the secure source only.
def frame(s):
 rows=list(z[s].values);return pd.DataFrame(rows[1:],columns=rows[0]).set_index('phenotype')
age=frame('Age');ehr=frame('EHR_observation');obs=list(z['Observation_criteria'].values);obs=pd.DataFrame(obs[1:],columns=obs[0]);life=list(z['Phenotype_lifecycle'].values);life=pd.DataFrame(life[1:],columns=life[0])
assert (life.loc[life.stage.eq('strict_analysis') & life.criterion.isin(['exactly_one_qualifying_record','exactly_one_qualifying_date']),'N']==0).all()
results=f"In the corrected primary cohort, median age was {age.loc['case','median']:.1f} years (IQR {age.loc['case','q1']:.1f}–{age.loc['case','q3']:.1f}) among cases and {age.loc['control','median']:.1f} years (IQR {age.loc['control','q1']:.1f}–{age.loc['control','q3']:.1f}) among controls. Among participants with an available EHR span, median recorded clinical history was {ehr.loc['case','median']:.2f} years (IQR {ehr.loc['case','q1']:.2f}–{ehr.loc['case','q3']:.2f}) for cases and {ehr.loc['control','median']:.2f} years (IQR {ehr.loc['control','q1']:.2f}–{ehr.loc['control','q3']:.2f}) for controls. Observation time was unavailable for {int(ehr.loc['control','N_missing']):,} controls; controls were not universally EHR-observed. Available zero-day spans were retained. The cohort lifecycle audit found no participant with only one qualifying record or one qualifying diagnosis date in the analyzed cohort. These checks evaluate implementation and ascertainment, not chart-adjudicated diagnostic validity."
main='''# Secure manuscript insertion candidate — cohort and laboratory Methods

HELD in Workbench. This is an insertion candidate for author review, not an exported/final manuscript. Existing corrected negative primary-X interpretation and discovery limitations remain unchanged.

## Results insertion

'''+results+'''

## Methods insertion

Psoriasis was defined by at least two distinct calendar dates of qualifying diagnoses after the source-specific cohort selection and provider sample-QC steps. The descendant-expanded concept set included psoriatic arthritis and non-plaque psoriasis; complete concepts and exclusions are supplied in Supplementary File S1. The primary control definition did not require a minimum observation period. We therefore performed a targeted sensitivity using EHR-observed controls with at least two non-study visit dates spanning 365 days, including a comparison restoring medication-only exclusions. The same cases and corrected sample-QC gates were retained. These analyses were exploratory phenotype sensitivities.

Age used the implemented 2023-10-01 reference date. Clinical observation span was the interval between the earliest and latest recorded condition, procedure, drug-exposure or qualifying visit date, with unavailable spans reported separately. A recorded span does not establish continuous clinical observation. Full cohort flow, case/control denominators and genetic-similarity assignments are reported in Supplementary File S1.

For biomarkers, retained measurements were converted and screened using the concept, unit and plausibility rules in Supplementary Files S1/S2. Component medians were calculated separately within participants. The original SII was the product of the platelet and neutrophil medians divided by the lymphocyte median. A separate sensitivity calculated SII from all three components at matching non-midnight timestamps before taking participant medians, under original and strict explicit cell-count unit rules. Matching timestamps could not confirm a common specimen. Measurements were not restricted to a period before psoriasis diagnosis or treatment; first observed diagnosis and medication exposure were assessed descriptively and do not establish onset or treatment-naive status. Full original and sensitivity processing code and its execution evidence are retained in the secure source collection.
'''
(V/'Manuscript_cohort_laboratory_insertions.private.md').write_text(main)
checks={'pipeline_manifest_verified':True,'retained_replay_receipts_checked':True,'prior_workbook_cells_and_styles_preserved':True,'new_table_cells_equal_sources':True,'old_S2_numerical_tables_preserved':True,'original_sii_source_identical':True,'original_v4_unchanged':True,'new_GWAS_or_model_fit_performed':False,'fresh_machine_replay':False,'native_Excel_Word_QA':False}
receipt={'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'SOURCE AND NUMERICAL CHECKS PASSED; VISUAL REVIEW PENDING','release':'HELD','checks':checks,'prior_S1_sha256':oldhash,'prior_S2_sha256':sha(old_s2),'builder_sha256':sha(Path(__file__)),'files':{str(p.relative_to(V)):sha(p) for p in V.rglob('*') if p.is_file()},'versions':{k:importlib.metadata.version(k) for k in ['pandas','numpy','openpyxl']},'limitations':['No independent clinical validation','No common-specimen confirmation','Original lipid-unit heuristics remain documented','No new corrected-cohort timing summary','Final scientific scope, native-file QA, approved release and author approval open']}
save(V/'BUILD.private.json',receipt)
print('V5_METHODS_BUILD_VERIFIED',len(source_rows),'sources',len(tables),'new tables',len(z.sheetnames),'total sheets')
