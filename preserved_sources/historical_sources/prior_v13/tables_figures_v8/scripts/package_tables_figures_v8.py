from pathlib import Path
import json,hashlib,csv,shutil,zipfile,platform
from datetime import datetime,timezone
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912');O=I/'tables_figures_v8';N=I.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
save=lambda p,r:Path(p).write_text(json.dumps(r,indent=2,default=str))
manifest=json.loads((O/'TABLE_FIGURE_MANIFEST.private.json').read_text());assert len(manifest['figures'])==28
revisions=[('Original Table 1 / review cohort tables','Table_1','Reconcile nested flow, genetic groups, age and EHR distributions; withhold protected intermediate cells.'),('Original Table 2 / review historical-X table','Table_2','Use frozen current primary estimates and explicit QC exclusions; QKI and NLGN4X remain subthreshold.'),('Original biomarker Table 3','Table_3','Use QC-passing, provider-pruned combined-cohort Bonferroni results; retain original test family.'),('Cloud author-review calibration Table 3','S1 Calibration; Figures 1 and S4a–f','Calibration is supplementary; Table 3 now denotes the biomarker table.'),('Original Figure 1','Figure_1','Corrected Manhattan and unit-weight full-output QQ; separate primary X02 and X01 sensitivity.'),('Original Figure 2 / Review Figure 1','Figure_2','Current autosomal reporting groups and historical-X/QKI status; no unsupported independent-locus claim.'),('Original Figure 3','Figure_3','Corrected exploratory sex/ancestry effects and heterogeneity; total eligible N distinguished from analyzed N.'),('Original Figure 4','Figure_4','QC-passing current biomarker results; distinguish Bonferroni from FDR-only.'),('Original Supplementary Figure 1','Supplementary_Figure_S1','Cases-only nominal biomarker findings remain exploratory; corrected q values visible.'),('Original sparse Supplementary Figure 2','Supplementary_Figure_S2a–e; S5a–d; S6a–b','Replace sparse annotated extracts with full region and LD plots. Current S2a is NLGN4X, S2b PTCHD1-AS, S2c/d Xq12, S2e QKI. Alternate dosage and sex-specific control LD are separate. No fine-mapping or colocalization claim.'),('Original Supplementary Figure 3','Supplementary_Figure_S3','Rebuild Y scan and full QQ from current output.'),('Additional calibration evidence','Supplementary_Figure_S4a–f','Full sex/ancestry QQ plus pooled MHC-excluded QQ.'),('Current conditional/LD analyses','Supplementary_Figure_S7a–b','Two explicit conditional test families and ancestry-specific autosomal LD; missing LD remains NA.'),('Sequencing center and phenotype sensitivity','Supplementary_Figure_S8a–b','Matched center and provider-pruned EHR-control comparisons; excluded historical targets remain diagnostic only.'),('Biomarker timing sensitivity','Supplementary_Figure_S9','Matched separate-median and same-draw SII comparisons; no equivalence claim.')]
with (O/'Table_figure_revision_map.private.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='	');w.writerow(['previous_asset','current_asset','required_change']);w.writerows(revisions)
readme='''# Psoriasis tables and figures — secure v8

This folder is the current corrected table/figure candidate. Files remain HELD in Workbench. No scientific assets were exported in this revision. It does not replace the external Word/Drive manuscript until permitted transfer, coordinated disclosure review and native integration are completed.

## Files

- S1_reconciled_tables_HELD.xlsx: 47 sheets, preserving all 44 prior sheets and adding three editable main tables. Original 57-sheet XLSM remains separately preserved; this is not a macro conversion.
- tables/: three TSV tables and five PDF proof pages across three PDFs.
- figures/: four main figures and 24 supplementary files across S1–S9; each has PDF, PNG and 600-dpi RGB TIFF.
- Figure_legends.private.md, Table_legends.private.md and Table_figure_revision_map.private.tsv: final captions and exact replacement/numbering instructions.
- TABLE_FIGURE_MANIFEST.private.json: authoritative current hashes; older intermediate build receipts have been superseded in its stated precedence order.
- FINAL_NUMERICAL_AND_FILE_AUDIT.private.json and VISUAL_QA.private.json: numerical/file and visual checks.

## Interpretation and integration

The corrected primary male 0/2 scan has no genome-wide significant QC-passing X association. Historical HWE exclusions are distinguished from nonsignificant results. Male 0/1, sex/ancestry, conditional and biomarker findings remain sensitivity or exploratory analyses; within-cohort tests are not independent replication. Physical reporting groups do not by themselves establish independent signals.

Table 1 withholds intermediate case/control cells that would enable small-count subtraction. Exact denominator and cohort accounting remain in the secure S1; the full private workbook is not cleared for public release. Figure 3 distinguishes primary eligible cohort size from varying analyzed N. It retains current nominal heterogeneity and explicitly labels superseded historical values.

The main set has three tables plus four figures. Journal placement must be reconciled with the journal's suggested combined limit; no scientific table has been silently removed. Final editable Word table conversion, in-text cross-references, native Word/Excel page QA, complete S1/S2 coverage, permitted summary-statistic access instructions and author approval remain open. Use this new numbering consistently; do not append corrected figures beside contradictory historical panels.

The old S2 sparse gene-label figure is superseded by full regional/LD plots. Gene annotation, fine-mapping and colocalization are not implied by these plots. Both Xq12 targets remain one region, and the original omission of the second target from historical sex analyses is still distinguished from later corrected sensitivities.

## Verification limits

Main table cells match their structured receipts and TSVs. The preceding 1,184,954 workbook cells were preserved. Figure source checks and full association receipts were verified before plotting. QQ ranks and lambda use unit weights and full outputs; visual thinning is applied only after ranking. Original REGENIE P values are not replaced by reconstructed Wald P values. The figures and all five table proof pages were visually checked. These are saved plot/table proofs, not native Word/Excel rendering certification.

NLGN4X full-cohort mapping validation, any other retained-claim scientific limitations, fresh-machine/full-pipeline replay, permitted transport and disclosure, and PI/author decisions remain governed by START HERE. Eczema was not modified.
'''
(O/'README.md').write_text(readme)
(O/'S2_table_figure_integration_addendum.private.md').write_text('# S2 integration addendum — secure v8\n\nRetain the scientific methods and evidence in supplement_staging_v7/S2_reconciliation_integration.private.md. Replace its table/figure references using Table_figure_revision_map.private.tsv and the captions below. The new main table sheets supplement the existing S1 contents; the prior 44 worksheets are unchanged. Full private S1 requires coordinated disclosure review before release. This addendum is an integration candidate, not a native Word edit or author approval.\n\n'+(O/'Figure_legends.private.md').read_text()+'\n'+(O/'Table_legends.private.md').read_text())
qa=[]
for r in manifest['figures']:
 name=r['id'];special={'Figure_3':'final2_Figure_3.png','Supplementary_Figure_S7a':'final2_Supplementary_Figure_S7a.png','Supplementary_Figure_S7b':'final2_Supplementary_Figure_S7b.png'}
 proof=O/'review_proofs'/special.get(name,'final_'+name+'.png')
 assert proof.exists();qa.append({'id':name,'source_png_sha256':r['files']['png']['sha256'],'proof':str(proof),'proof_sha256':sha(proof),'review':'full figure visually checked; S7a final PNG identical to reviewed prior rendering' if name=='Supplementary_Figure_S7a' else 'full figure visually checked'})
# Confirm the S7a layout did not change while S7b margin was adjusted.
old=json.loads((O/'LAYOUT_REVISION_COMPLETE.private.json').read_text())['figures'];a=next(x for x in old if x['id']=='Supplementary_Figure_S7a');b=next(x for x in manifest['figures'] if x['id']=='Supplementary_Figure_S7a');assert a['files']['png']['sha256']==b['files']['png']['sha256']
tableproofs=['Table_1_final_01.png','Table_1_02.png','Table_2_01.png','Table_2_02.png','Table_3_01.png']
save(O/'VISUAL_QA.private.json',{'utc':datetime.now(timezone.utc).isoformat(),'figures':qa,'tables':[{'proof':str(O/'review_proofs'/x),'sha256':sha(O/'review_proofs'/x),'review':'visually checked'} for x in tableproofs],'corrections':['chromosome tick crowding','panel-letter margins','ancestry legend completeness','sex model dosage title','log-axis minor-label crowding','Figure 3 title margin','LD color scale clipping','intermediate cohort count disclosure'],'scope':'Saved plot and selected main-table proofs only; no native Word/Excel certification.'})
# Preserve resolved intermediate failure records without presenting them as current failures.
history=O/'build_history';history.mkdir(exist_ok=True)
for p in O.glob('*FAILED.private.json'):shutil.move(str(p),str(history/p.name))
(history/'RESOLVED.md').write_text('Earlier build failures were resolved: QKI source-column fallback, openpyxl read-only style identifier, and final Figure 3 metadata return handling. Final audit is authoritative. Original fully opaque RGBA TIFFs were preserved before lossless RGB normalization; plotted RGB pixels were checked identical.\n')
# Preserve every prior source archive entry under prior_v12 and add the scripts used for this revision.
oldzip=N/'code_release_candidate_20260912_v12/psoriasis_code_source_collection_v12_HELD.zip';dest=N/'code_release_candidate_20260912_v13';dest.mkdir(exist_ok=True);out=dest/'psoriasis_code_source_collection_v13_HELD.zip'
contents={}
with zipfile.ZipFile(oldzip) as z:
 assert z.testzip() is None
 for n in z.namelist():contents['prior_v12/'+n]=z.read(n)
scriptnames=['update_tables_figures_v8.py','update_forests_v8.py','update_supplement_figures_v8.py','update_main_tables_v8.py','finalize_figure_layout_v8.py','prepare_final_proofs_v8.py','audit_tables_figures_v8.py','normalize_tiff_v8.py','package_tables_figures_v8.py']
for n in scriptnames:contents['tables_figures_v8/scripts/'+n]=(I/'scripts'/n).read_bytes()
contents['README.md']=b'Secure source collection v13. All v12 entries preserved under prior_v12; new table/figure scripts added. No scientific tables, figures or participant files added. This is not fresh-machine/full-pipeline replay or release certification. Paths and permissions must be reconciled within the authorized environment.\n'
import importlib.metadata
versions={n:importlib.metadata.version(n) for n in ['numpy','pandas','matplotlib','Pillow','openpyxl','scipy']};versions['python']=platform.python_version()
contents['tables_figures_v8/software_versions.json']=json.dumps(versions,indent=2).encode()
entries=[{'path':n,'sha256':hashlib.sha256(v).hexdigest(),'bytes':len(v)} for n,v in sorted(contents.items())]
contents['MANIFEST.json']=json.dumps({'entries':entries,'prior_archive_sha256':sha(oldzip),'scope':'source collection only'},indent=2).encode()
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for n,v in contents.items():z.writestr(n,v)
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None
 for r in entries:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
 with zipfile.ZipFile(oldzip) as prior:
  for n in prior.namelist():assert z.read('prior_v12/'+n)==prior.read(n)
receipt={'file':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'entries':len(contents),'prior_entries_preserved':289,'CRC_and_manifest_verified':True,'versions':versions};save(dest/'ARCHIVE_RECEIPT.json',receipt)
checkpoint={'utc':datetime.now(timezone.utc).isoformat(),'status':'Corrected table/figure candidate complete in secure Workbench; release and native manuscript integration remain open','tables_figures':str(O),'manifest_sha256':sha(O/'TABLE_FIGURE_MANIFEST.private.json'),'workbook_sha256':sha(O/'S1_reconciled_tables_HELD.xlsx'),'figures':28,'main_tables':3,'sheets':47,'visual_QA_sha256':sha(O/'VISUAL_QA.private.json'),'prior_S2':str(I/'supplement_staging_v7/S2_reconciliation_integration.private.md'),'S2_addendum':str(O/'S2_table_figure_integration_addendum.private.md'),'source_archive':receipt,'no_scientific_export':True,'remaining':['permitted transfer and coordinated disclosure','native Word/S1/S2 integration and final journal placement/formatting','full-cohort mapping and remaining retained-claim scientific evidence or author-approved limitations','full-pipeline/fresh-machine replay and approved access identifiers','PI and author approvals'],'eczema':'untouched'}
save(I/'CONTINUATION_CHECKPOINT_20260913_v8.private.json',checkpoint)
print('V8_PACKAGE_COMPLETE',28,3,47,'ARCHIVE',len(contents),flush=True)

