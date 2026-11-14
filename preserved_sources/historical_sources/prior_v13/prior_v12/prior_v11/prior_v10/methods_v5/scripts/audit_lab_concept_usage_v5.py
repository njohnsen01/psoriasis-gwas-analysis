from pathlib import Path
import ast,datetime,hashlib,json
import pandas as pd
P=Path('/home/jupyter/psoriasis_reproducibility_20260910/pipeline_v3');R=P.parent;D=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912/supplement_staging_v5/methods_sources')
f=R/'fresh_lab_processing_v3/psoriasis_biomarker_labs_classified.tsv';ids=set()
for n in ast.walk(ast.parse((P/'cells/cell_097.py').read_text())):
 if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id.endswith('_IDS'):ids.update(ast.literal_eval(n.value))
counts={};rows=0
for chunk in pd.read_csv(f,sep=chr(9),usecols=['measurement_concept_id'],dtype=str,chunksize=500000):
 rows+=len(chunk)
 for k,v in chunk.measurement_concept_id.value_counts().items():counts[k]=counts.get(k,0)+int(v)
assert set(counts)<=ids
r={'checked_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':str(f.resolve()),'source_bytes':f.stat().st_size,'declared_concepts':sorted(ids,key=int),'observed_concepts':sorted(counts,key=int),'unobserved_declared_concepts':sorted(ids-set(counts),key=int),'unmatched_vocabulary_concept_rows':counts.get('1497147',0),'classified_rows':rows,'row_counts_by_concept':counts,'release':'HELD'}
(D/'lab_concept_usage.private.json').write_text(json.dumps(r,indent=2)+chr(10));print('CLASSIFIED_CONCEPT_USAGE_VERIFIED; unmatched-ID rows',r['unmatched_vocabulary_concept_rows'])
