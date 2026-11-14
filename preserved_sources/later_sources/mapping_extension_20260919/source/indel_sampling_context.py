import pandas as pd,json
from scipy.stats import binomtest
from mapping_extension import O
p=O;d=pd.read_csv(p/'indel_classifier_review.private.tsv',sep='	');f=d.loc[d.original_balance_flag].copy();h=f.loc[f.expected.eq(.5)].copy()
h['binomial_P']=[binomtest(round(r.canonical_fraction*r.canonical_denominator),int(r.canonical_denominator),.5).pvalue for r in h.itertuples()];h.to_csv(p/'flagged_heterozygote_binomial.private.tsv',sep='	',index=False)
a=dict(all_deletion_balance_flags_are_heterozygotes=bool(f.expected.eq(.5).all()),any_flagged_heterozygote_nominal_binomial_P_lt_05=bool(h.binomial_P.lt(.05).any()),any_flagged_heterozygote_bonferroni_across_all_reviewed_target_rows=bool(h.binomial_P.lt(.05/(3*len(d))).any()),interpretation='Post hoc descriptive sampling check for fixed diagnostic flags, not removal of flags or modification of GWAS QC. Exact participant counts and individual statistics retained privately.')
(p/'INDEL_FLAG_SAMPLING_CONTEXT.json').write_text(json.dumps(a,indent=2));print(json.dumps(a,indent=2))

