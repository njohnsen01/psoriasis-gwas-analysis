from pathlib import Path
import json,subprocess,datetime
import numpy as np,pandas as pd
from scipy.stats import fisher_exact
from mapping_extension import O,REF,sha
from read_qc_core import TARGETS,observe,unique_high_quality
protocol=dict(purpose='Additional descriptive per-allele read-quality diagnostics; no association rerun or QC exclusions',metrics='Median BQ, read-cycle fraction, soft clipping and MQ per REF/ALT; strand Fisher exact only in heterozygotes with both alleles; preserve missing metrics',interpretation='Read-position and strand diagnostics are technical descriptors. They neither prove mapping correctness nor change prior fixed balance flags.',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'ALLELE_QC_PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
ref=''.join(subprocess.check_output(['samtools','faidx',str(REF),'chrX:5996489-5996506'],text=True).splitlines()[1:]);window=(5996489,5996506,ref,ref[:6]+ref[8:])
R=json.loads((O/'results.private.json').read_text());out=[]
for r in R:
 fields=[l.split('	') for l in (O/'reads_private'/r['review_id']/'region.private.sam').read_text().splitlines() if not l.startswith('@')]
 for j,target in enumerate(TARGETS):
  rows=[z for f in fields if (z:=observe(f,target,window if j==1 else None)) is not None];high,_=unique_high_quality(rows);t=r['targets'][j]
  item=dict(IID=r['IID'],target='g'+str(j),sex=r['sex'],ancestry=r['ancestry'],case_status=r['case_status'],expected=t['expected_target_fraction']);strands=[]
  for allele in ['REF','ALT']:
   z=[x for x in high if x['call']==allele];strands.append([sum(not x['reverse'] for x in z),sum(x['reverse'] for x in z)])
   for metric,values in [('BQ',[x['BQ'] for x in z]),('MQ',[x['MQ'] for x in z]),('cycle_fraction',[x['cycle']/x['read_length'] for x in z]),('soft_clip',[x['soft_clipped'] for x in z])]:item[allele+'_'+metric]=float(np.median(values)) if values else np.nan
  item['strand_P']=float(fisher_exact(strands).pvalue) if item['expected']==.5 and all(sum(x)>0 for x in strands) else np.nan
  item['ALT_minus_REF_cycle_fraction']=item['ALT_cycle_fraction']-item['REF_cycle_fraction'];out.append(item)
d=pd.DataFrame(out);d.to_csv(O/'allele_quality.private.tsv',sep='	',index=False);summary=[]
keys=['target','sex','expected'];metrics=['REF_BQ','ALT_BQ','REF_MQ','ALT_MQ','REF_cycle_fraction','ALT_cycle_fraction','ALT_minus_REF_cycle_fraction','REF_soft_clip','ALT_soft_clip']
for key,z in d.groupby(keys):
 if len(z)<=20:continue
 row=dict(zip(keys,key));row['size']='>20; exact denominator withheld'
 for metric in metrics:
  v=z[metric].dropna()
  if len(v)>20:row[metric+'_median']=float(v.median());row[metric+'_q25']=float(v.quantile(.25));row[metric+'_q75']=float(v.quantile(.75))
 summary.append(row)
pd.DataFrame(summary).to_csv(O/'reportable_allele_quality.tsv',sep='	',index=False)
audit=dict(status='allele_quality_complete',all_selected_targets_reviewed=len(d)==3*len(R),heterozygote_strand_flags_after_all_reviewed_target_rows_bonferroni=bool(d.strand_P.lt(.05/len(d)).any()),association_analysis_unchanged=True,limits='Strand tests are post hoc technical diagnostics. Absent flags do not exclude mapping artifacts. No participant counts or individual metrics are released.',source_sha256=sha(Path(__file__)))
(O/'ALLELE_QC_AUDIT.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))

