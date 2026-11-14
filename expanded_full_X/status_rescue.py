"""Aggregate-only status; no raw logs, participant cells or variant-level records."""
from pathlib import Path
import json,datetime,pandas as pd
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'full_x_rescue'
print('UTC',datetime.datetime.now(datetime.timezone.utc).isoformat())
for name in ['diagnose_all_x','refine_benchmark_reference','chain_full_x_rescue','rescue_x_firth_all','run_final_lead_checks']:
 p=R/(name+'.pid')
 if not p.exists():continue
 try:active='\nState:\tZ' not in (Path('/proc')/p.read_text().strip()/'status').read_text()
 except FileNotFoundError:active=False
 log=R/(name+'.console.log');print('JOB',name,'pid',p.read_text().strip(),'active',active,'traceback',('Traceback' in log.read_text()) if log.exists() else False)
for folder in ['benchmark','recovery','firth_recovery','R_ml_benchmark_exact','final_lead_checks']:
 p=O/folder;files=list(p.glob('chunk_*.private.tsv'))
 if files:
  d=pd.concat([pd.read_csv(f,sep='\t',usecols=['status']) for f in files]);print('PROGRESS',folder,'chunks',len(files),'tests',len(d),'status',d.status.value_counts().to_dict())
 else:
  files=[f for f in p.glob('g*.csv') if '.numerics.' not in f.name]
  if files:
   d=pd.concat([pd.read_csv(f,usecols=['status']) for f in files]);print('PROGRESS',folder,'fits',len(d),'status',d.status.value_counts().to_dict())
print('AUDIT_CHUNKS',len(list((R/'full_x_estimability_audit').glob('chunk_*.private.tsv'))))
for name in ['R_ML_BENCHMARK_VALIDATED.json','BENCHMARK_VALIDATED.json','ML_RECOVERY_STARTED.json','RECOVERY_COMPUTED.json','FIRTH_BENCHMARK_VALIDATED.json','FIRTH_RECOVERY_COMPUTED.json','FINAL_LEAD_CHECKS_COMPUTED.json','RESCUE_AND_LEADS_COMPUTED.json']:
 p=O/name
 if p.exists():
  d=json.loads(p.read_text());print('CHECKPOINT',name,json.dumps({k:d[k] for k in ['status','variants','attempted','positive_controls','separation_negative_controls','status_counts','all_status_counts','max_R_Wald_difference','max_separated_R_difference','selected','ML_verified','Firth_verified','ML_threshold_crossers_in_combined_scan','Firth_sensitivity_threshold_crossers'] if k in d}))
p=R/'full_x_estimability_audit/DIAGNOSTICS_COMPLETE.json'
if p.exists():print('FULL_ESTIMABILITY_AUDIT',p.read_text())
