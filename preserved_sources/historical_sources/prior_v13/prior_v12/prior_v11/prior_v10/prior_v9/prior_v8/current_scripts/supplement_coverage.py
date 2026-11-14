from pathlib import Path
import json, hashlib, datetime
import pandas as pd
import openpyxl
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912')
O=I/'supplement_coverage'; O.mkdir(exist_ok=True)
original=Path('/home/jupyter/S1_current57_preserved_20260910.xlsm'); current=I/'supplement_staging_v2/S1_reconciled_tables_HELD.xlsx'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(original)=='41416a9c98f0d4c830f26efc87870ada606c3d3f55a7a0d5692b70b947605734'
a=openpyxl.load_workbook(original,read_only=True,keep_vba=True); b=openpyxl.load_workbook(current,read_only=True)
assert len(a.sheetnames)==57; rows=[]
mapping={1:['README'],2:['Cohort_flow','Ancestry_groups','Followup_groups','Age','EHR_observation','Observation_criteria','Phenotype_lifecycle'],3:['Autosomal_GWS'],4:['Autosomal_reporting_groups'],5:['Exact_primary_estimates','Target_QC','Autosomal_GWS'],6:['X_followup'],7:['X_followup'],8:['Ancestry_followup','X_followup'],9:['Ancestry_heterogeneity'],10:['Ancestry_followup'],11:['Ancestry_heterogeneity'],12:['Y_male02','Y_male01','Y_calibration'],14:['Target_QC','Genotype_quality','X_alleles_missingness'],15:['X_followup','Source_manifest'],16:['Biomarkers','SII_sensitivity'],17:['Biomarkers','SII_sensitivity'],20:['Autosomal_GWS'],21:['Autosomal_GWS'],22:['Autosomal_reporting_groups'],23:['Exact_primary_estimates','Target_QC','Autosomal_GWS'],24:['X_followup'],25:['Ancestry_followup','X_followup'],26:['Y_male02','Y_male01'],27:['Biomarkers'],28:['Biomarkers'],29:['Autosomal_reporting_groups','Target_QC'],30:['Calibration','Y_calibration'],31:['Exact_primary_estimates','X_followup','Biomarkers'],32:['SII_sensitivity','EHR_sensitivity'],33:['Y_male02','Y_male01','Y_calibration'],34:['Y_male02','Y_male01'],36:['Exact_primary_estimates','Target_QC'],37:['Ancestry_heterogeneity'],38:['Ancestry_heterogeneity'],39:['Ancestry_followup'],40:['Ancestry_followup','Source_manifest'],46:['Regional_assets'],47:['README','Source_manifest'],55:['Target_QC','X_alleles_missingness'],56:['Genotype_quality','NLGN4X_read_pilot'],57:['Source_manifest']}
evidence={13,18,19,35,44,45,48,49,50,51,52,53,54}; historical_conditional={41,42,43}
assert set(mapping)|evidence|historical_conditional==set(range(1,58))
extra={7:[I/'sex_heterogeneity_current/sex_heterogeneity.private.tsv',I/'sex_heterogeneity_current/PROTOCOL.private.json'],10:[I/'QKI_primary_audit/COMPLETE.private.json',I/'regional_QKI/COMPLETE.private.json'],46:[I/'regional_X/COMPLETE.private.json',I/'regional_QKI/COMPLETE.private.json'],55:[I/'read_mapping_QC/genotype_quality_COMPLETE.private.json'],56:[I/'read_mapping_QC/CRAM_pilot/REVIEW.private.json'],57:[I.parent/'code_release_candidate_20260912_v7/ARCHIVE_RECEIPT.json',I/'clean_process_check_v3/COMPLETE.json']}
for idx,name in enumerate(a.sheetnames,1):
    tabs=mapping.get(idx,[]); assert all(t in b.sheetnames for t in tabs)
    refs=extra.get(idx,[]); assert all(p.exists() for p in refs)
    if idx in evidence:
        status='PRESERVED_EVIDENCE; CORRECTED-LEAD/CURRENT-SCOPE REVIEW OPEN'; note='Original evidence retained under its exact original sheet name and workbook hash. No wholesale reuse as current-lead evidence or first-discovery support; review source dates, corrected variants and original association sources.'
    elif idx in historical_conditional:
        status='HISTORICAL CONDITIONAL SOURCE; CURRENT-SCOPE DECISION OPEN'; note='Original chromosome 5 conditional outputs preserved. Current reporting groups are not labeled independent loci. If retaining an independence claim, reconcile current source/QC, conditioning dosage and relatedness and integrate its exact estimates.'
    else:
        status='CURRENT NUMERICAL SOURCES STAGED; FINAL PRESENTATION/RELEASE OPEN'; note='Source rows available in linked current tabs; this is a mapping of required content, not a claim that original rows/labels should be copied or that final submission coverage is complete.'
    if idx==3: note+=' Derive the final MHC display from current chromosome 6 region rows; no HLA fine-mapping or historical total substitution.'
    if idx==7: note+=' Later post-hoc corrected-model sex contrasts are separate from the historical three-test family; do not relabel the correction as prespecified.'
    if idx==38: note+=' Current Q/I-squared heterogeneity does not replace historical leave-one-group-out estimates. Retain the latter only with current-source reconciliation or explicitly remove that promise.'
    if idx in {13,32,47,57}: note+=' Final S2 must include actual cohort/phenotype and laboratory concept/unit/timing rules and executable stage/access instructions, not only narrative summaries.'
    rows.append(dict(original_order=idx,original_sheet=name,original_workbook=str(original),original_sha256=sha(original),current_workbook=str(current),current_sha256=sha(current),current_tabs=' | '.join(tabs),additional_secure_sources=' | '.join(str(p) for p in refs),additional_source_hashes=' | '.join(sha(p) for p in refs),status=status,remaining=note,release='HELD'))
out=pd.DataFrame(rows); target=O/'original57_to_current_sources.private.tsv'; out.to_csv(target,sep=chr(9),index=False)
assert pd.read_csv(target,sep=chr(9))['original_sheet'].tolist()==a.sheetnames; a.close(); b.close()
receipt=dict(status='all_original_sheets_accounted_for_with_explicit_open_items',original_sheet_count=57,original_sha256=sha(original),current_workbook_sha256=sha(current),mapping_sha256=sha(target),script_sha256=sha(__file__),final_submission_coverage_certified=False,release='HELD')
(O/'COMPLETE.private.json').write_text(json.dumps(receipt,indent=2)); print(json.dumps(receipt))