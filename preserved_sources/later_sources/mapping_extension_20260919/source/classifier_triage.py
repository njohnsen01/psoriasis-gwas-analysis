from pathlib import Path
import json,subprocess
import pandas as pd
from collections import Counter
from mapping_extension import O,REF,sha
from read_qc_core import TARGETS,observe,unique_high_quality
ref=''.join(subprocess.check_output(['samtools','faidx',str(REF),'chrX:5996489-5996506'],text=True).splitlines()[1:]);assert ref[5:8]=='GCT';alt=ref[:6]+ref[8:];window=(5996489,5996506,ref,alt)
R=json.loads((O/'results.private.json').read_text());out=[];ymates=[]
for r in R:
 fields=[l.split('	') for l in (O/'reads_private'/r['review_id']/'region.private.sam').read_text().splitlines() if not l.startswith('@')]
 for j,target in enumerate(TARGETS):
  rows=[z for f in fields if (z:=observe(f,target,window if j==1 else None)) is not None];high,_=unique_high_quality(rows)
  ymates.append(dict(IID=r['IID'],target='g'+str(j),HQ_ALT_Y_mate=any(x['call']=='ALT' and x['mate_chrY'] for x in high),HQ_ALT_Y_alternative=any(x['call']=='ALT' and x['alternative_chrY'] for x in high)))
  if j!=1:continue
  t=r['targets'][j];counts=Counter(x['call'] for x in high);den=counts['ALT']+counts['REF'];canonical=counts['ALT']/den if den else None
  exp=t['expected_target_fraction'];pass_c=canonical is not None and (canonical<=.1 if exp==0 else canonical>=.9 if exp==1 else .25<=canonical<=.75)
  out.append(dict(IID=r['IID'],original_balance_flag=not t['diagnostic_genotype_concordance'],adequate=t['adequate_HQ_depth'],expected=exp,strict_fraction=t['HQ_target_fraction'],OTHER_reads=counts['OTHER'],canonical_fraction=canonical,canonical_denominator=den,flag_persists_using_REF_ALT_only=not pass_c))
pd.DataFrame(out).to_csv(O/'indel_classifier_review.private.tsv',sep='	',index=False);pd.DataFrame(ymates).to_csv(O/'Y_mate_review.private.tsv',sep='	',index=False)
z=pd.DataFrame(out);f=z.loc[z.original_balance_flag];y=pd.DataFrame(ymates)
summary=dict(indel_strict_span_flags_present=not f.empty,all_indel_balance_flags_have_OTHER_haplotypes=bool(f.OTHER_reads.gt(0).all()),any_indel_flag_persists_when_only_REF_ALT_reads_used=bool(f.flag_persists_using_REF_ALT_only.any()),all_flagged_canonical_denominators_at_least_10=bool(f.canonical_denominator.ge(10).all()),HQ_ALT_chrY_mate_present=bool(y.HQ_ALT_Y_mate.any()),HQ_ALT_chrY_alternative_tag_present=bool(y.HQ_ALT_Y_alternative.any()),interpretation='Classifier sensitivity only; removing OTHER from the denominator cannot prove genotype correctness. No GWAS changes. Participant counts withheld.')
(O/'CLASSIFIER_TRIAGE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

