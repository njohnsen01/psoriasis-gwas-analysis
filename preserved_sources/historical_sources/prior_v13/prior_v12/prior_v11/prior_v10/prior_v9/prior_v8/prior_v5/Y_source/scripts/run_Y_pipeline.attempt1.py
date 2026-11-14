from pathlib import Path
import argparse,json,subprocess,sys
Y=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--run-missing',action='store_true');p.add_argument('--download-missing',action='store_true');args=p.parse_args()
prep=[sys.executable,str(Y/'scripts/prepare_Y_inputs.py')]
if args.download_missing:prep.append('--download-missing')
subprocess.run(prep,check=True)
stages=[('site_filters','source_Y_site_filters.py',[]),('counts','source_Y_counts.py',[]),('projection_pilot','project_Y_source.py',['--pilot']),('projected','project_Y_source.py',[]),('results','run_Y_models.py',[])]
for folder,script,extra in stages:
 receipt=Y/folder/'COMPLETE.json'
 if not receipt.exists():
  assert args.run_missing,'Stage '+folder+' missing; use --run-missing for execution'
  cmd=[sys.executable,str(Y/'scripts'/script)]+extra
  with (Y/(folder+'.orchestrator.console.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 r=json.loads(receipt.read_text());assert r['status'] in ['complete','validated']
 assert not (Y/folder/'FAILED.json').exists(),'Inspect failed-stage record before proceeding'
 print('Y_STAGE_VERIFIED',folder,flush=True)
print('Y_PIPELINE_STAGE_RECEIPTS_VERIFIED; quantitative integrity audit is retained in the freeze record')
