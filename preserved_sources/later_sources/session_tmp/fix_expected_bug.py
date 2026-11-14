import json, collections
from pathlib import Path
O = Path('.')
correct_expected = {
    'reference_AAG': 'REF',
    'SNP_only_GAG': 'ALT',
    'MNV_GAT': 'ALT',
    'neighbor_only_AAT': 'REF',
}
bug_note = "BUG FOUND AND FIXED: the first run's 'expected' truth label incorrectly flagged neighbor_only_AAT as ALT (used 'OR s[1002]!=ref[1002]' instead of matching the actual g2/rs183018930 classification target, which is defined ONLY by position 1000). This inflated apparent error rates uniformly across both aligners and both edge/mid buckets by exactly 25% (1 of 4 haplotypes mislabeled). Recomputed here from haplotype identity alone, matching the ORIGINAL synthetic_haplotype_controls.py convention (expected='ALT' if s[1000]=='G' else 'REF')."
summary_all = []
for model in ['minimap2','bwa_mem2']:
    details = json.loads((O/(model+'_fullref_synthetic_results.json')).read_text())
    for d in details:
        d['expected_FIXED'] = correct_expected[d['haplotype']]
        d['correct_FIXED'] = (d['observed'] == d['expected_FIXED'])
    (O/(model+'_fullref_synthetic_results.json')).write_text(json.dumps(details, indent=2))
    for h in sorted(correct_expected):
        v = [x for x in details if x['haplotype']==h]
        counts = collections.Counter(x['observed'] for x in v)
        summary_all.append(dict(model=model, haplotype=h, n_templates=len(v), n_correct=sum(x['correct_FIXED'] for x in v),
                                 REF=counts['REF'], ALT=counts['ALT'], OTHER=counts['OTHER'], NO_HQ_CALL=counts['NO_HQ_CALL']))
    for bucket_name, pred in [('edge', lambda c: c <= 15 or c >= 135), ('mid', lambda c: 15 < c < 135)]:
        vv = [x for x in details if pred(x['cycle'])]
        summary_all.append(dict(model=model, haplotype='__BUCKET__'+bucket_name, n_templates=len(vv), n_correct=sum(x['correct_FIXED'] for x in vv),
                                 REF=None, ALT=None, OTHER=None, NO_HQ_CALL=sum(1 for x in vv if x['observed']=='NO_HQ_CALL')))
result = dict(bug_note=bug_note, summary=summary_all)
(O/'FULL_REF_EDGE_EXTENSION_RESULTS.json').write_text(json.dumps(result, indent=2))
print(json.dumps([s for s in summary_all], indent=2))
