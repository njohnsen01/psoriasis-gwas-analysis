PSORIASIS GWAS ANALYSIS CODE - 13 NOVEMBER 2026

Version 1.0.0. This archive contains the preserved analysis source code. No license has been assigned.

Contents: 313 preserved source files, including 274 Python/R/SQL sources staged on 3 October and 38 sources for the subsequent expanded-X analysis and the reviewed October 14 matched-refit script. No participant data, genotypes, model objects, raw logs, private result tables, credentials, or input datasets are included. Original protected files remain in All of Us. Three source files build workspace storage URIs from WORKSPACE_BUCKET; SOURCE_MANIFEST.csv records both original and packaged hashes. Historical absolute execution paths are retained for traceability and must be configured in an authorized environment.

Do not batch-run this collection. Historical attempts, validation steps and superseded implementations are preserved because the completed analysis used checkpoints and scientific acceptance gates. Their presence is not a claim that every script produced a reported result. No association analyses were rerun to prepare this archive.

CODE MAP
- Matched 29-variant refits: matched_refit_20261014/fit_joint.py. Reviewed against the recovered Workbench file on 13 November 2026. The copy is byte-identical. Controlled inputs and generated results are external to the archive.
- Primary autosomal results: preserved_sources/historical_sources/.../replacement_primary/scripts/run_autosomes.py, run_autosomes_parallel_resume.py and run_autosomes_stratified.py. The manifest identifies exact paths and distinguishes duplicate historical snapshots. Final source hashes are also recorded in Supplementary Text S3.
- Original two-degree-of-freedom joint X scan and target diagnostics: preserved_sources/later_sources/joint_G_sex_20260914/.
- Targeted ancestry/age models: preserved_sources/later_sources/deeper_analysis_20260919/ and targeted_model_comparison_20260919/. These are not the October 14 matched 29-variant refits.
- Final blood-cell tests: preserved_sources/later_sources/SII_final_cohort_20260914/run_sii_final.py. Earlier scripts in that folder preserve the run history.
- Historical biomarker/HDL work: historical_sources/.../relatedness_sensitivity_20260910/run_biomarkers.py and original_biomarker_cell312.py. This is distinct from the later final blood-cell correction family.
- Technical sensitivities: later_sources/technical_followup_20260919/, mapping_extension_20260919/ and locus_specific_validation_20260919/.
- Expanded full-X follow-up: expanded_full_X/ contains original scan, failure diagnostics, independent validation, ML recovery, Firth sensitivity, lead checks and reviewed-export logic. chain_full_x_rescue.py and review_full_x_rescue.py document the recovery/review sequence. Files containing before_information_refinement preserve earlier versions; they are not final validation. Initial software-valid estimates were not all accepted as final estimates. The final publication counts are 294,445 ordinary and 40,481 Firth sensitivity estimates, with the original 334,926-variant family retained.

REPRODUCTION REQUIREMENTS
Use the matching All of Us Controlled Tier v8 phenotype, QC-selected samples, genotype files, covariate order, sex/ploidy rules, LOCO predictor where specified, ancestry definitions, model family and checkpoint files. A separate authorized researcher must obtain access to controlled inputs. Generated participant records and restricted intermediate outputs must remain within the authorized environment. The scripts' output-writing code does not authorize exporting those outputs. Reference-data paths and software must be configured; no dependency lockfile or complete environment image was recovered. Reported versions include REGENIE 4.1.2, PLINK 2.0.0-a6.9LM (29 January 2025), logistf 1.26.1, minimap2 2.31-r1302 and BWA-MEM2 2.3. Additional imported dependencies are listed in IMPORT_INVENTORY.txt; versions are not guessed.

OPEN REPRODUCIBILITY ITEMS
1. Missing original October calibration and simulation source files: firth.py, nullsim.py, confsim.py, sweep.py, qqrun.py, mkfig.py. All six are absent. See MISSING_SIMULATION_SOURCES.tsv.
2. Original full pooled-X scan minimum P values and original allele identities at chrX:10382681 and chrX:86808645 still require protected source outputs. The known non-PAR primary input subset minimum is not the complete scan minimum.
3. This collection does not establish that every historical analysis command, reference and environment has been preserved. End-to-end reproduction has not been tested from this archive.

LOCAL CONFIGURATION FOR REVIEWED SOURCES
PSO_WORKSPACE_MOUNT: local path to the authorized workspace mount.
WORKSPACE_BUCKET: authorized workspace bucket name or gs:// URI.
GOOGLE_PROJECT: authorized billing project ID.
PSO_PRIMARY_RUN_ID and PSO_X_FOLLOWUP_RUN_ID: matching run IDs without attempt_ or production_ prefixes.
PGENLIB_SITE_DIR: Python site directory containing the validated pgenlib 0.95.1 installation.
Provide the values needed by the script being used. Original and packaged SHA-256 digests are recorded in SOURCE_MANIFEST.csv.
