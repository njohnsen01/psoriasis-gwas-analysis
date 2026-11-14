from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';R=Path('/home/jupyter/psoriasis_reproducibility_20260910')
G=pd.read_pickle(R/'reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');M=pd.read_csv(R/'reconstructed_inputs_v2/variant_map.tsv',sep='\t');Q=pd.read_csv(I/'estimate_audit/target_QC.private.tsv.gz',sep='\t');E=pd.read_csv(I/'estimate_audit/exact_association_rows.private.tsv.gz',sep='\t')
checks=[]
for m in M.loc[~M.variant_id.str.startswith('chrX:')].itertuples():
 rows=Q[Q.ID.eq(m.variant_id)]
 assert len(rows)>0,m.variant_id
 er=E[E.ID.eq(m.variant_id)]
 if len(er): assert er.ALLELE1.eq(m.effect_allele).all(),(m.variant_id,'effect allele')
 else: assert not rows['PASS'].fillna(False).astype(bool).any(),(m.variant_id,'missing passing association')
 for _,q in rows.iterrows():
  group=q['group'];mask=pd.Series(True,index=G.index) if group=='pooled' else G.sex.eq(0 if group=='female' else 1) if group in ['female','male'] else G.ancestry_pred.str.upper().eq(group)
  assert mask.any(),group
  for kind,val in [('case',1),('control',0)]:
   v=G.loc[mask&G.psoriasis_case.eq(val),m.aligned_column];called=v.notna();ac=float(v.sum());an=int(called.sum()*2);mac=min(ac,an-ac);cr=float(called.mean())
   assert np.isclose(mac,q[kind+'_MAC'],atol=1e-6,rtol=0),(m.variant_id,group,kind,'MAC')
   assert np.isclose(cr,q[kind+'_callrate'],atol=1e-6,rtol=0),(m.variant_id,group,kind,'callrate')
   if pd.notna(q[kind+'_AC']): assert np.isclose(ac,q[kind+'_AC'],atol=1e-6,rtol=0)
   if pd.notna(q[kind+'_AN']): assert an==q[kind+'_AN']
   checks.append({'ID':m.variant_id,'group':group,'phenotype':kind,'allele':m.effect_allele,'MAC':mac,'callrate':cr,'status':'matched','primary_association_available':bool(len(er)),'primary_QC_PASS':bool(q['PASS'])})
out=I/'followup_input_audit';pd.DataFrame(checks).to_csv(out/'autosomal_allele_count_comparison.private.tsv',sep='\t',index=False)
rec={'status':'passed','targets':21,'comparisons':len(checks),'checks':'Case/control MAC and call rate; AC/AN where available. ALLELE1 verified on tested rows. QC-failed targets with no association row remain excluded.','source_files':[]}
for p in [R/'reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl',R/'reconstructed_inputs_v2/variant_map.tsv',I/'estimate_audit/target_QC.private.tsv.gz',I/'estimate_audit/exact_association_rows.private.tsv.gz']:
 rec['source_files'].append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(out/'AUTOSOMAL_ALLELE_CHECK.private.json').write_text(json.dumps(rec,indent=2));print('AUTOSOMAL_ALLELE_AUDIT_PASSED',len(checks))
