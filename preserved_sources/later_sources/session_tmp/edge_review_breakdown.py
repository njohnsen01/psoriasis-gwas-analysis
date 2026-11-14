import json
rows = json.loads(open('read_transitions.private.json').read())
INFORMATIVE = ('HQ_REF','HQ_ALT')
conflict_reviews = set()
flagged_reviews = set()
flagged_reviews_with_conflict = set()
all_reviews_with_conflict_nonflagged = set()
for r in rows:
    bwa = r.get('bwa'); mm = r.get('minimap')
    rid = r.get('review_id')
    if r.get('flagged'):
        flagged_reviews.add(rid)
    if not bwa or not mm:
        continue
    bstat = bwa.get('status'); mstat = mm.get('status')
    if bstat not in INFORMATIVE or mstat not in INFORMATIVE:
        continue
    if bstat != mstat:
        conflict_reviews.add(rid)
        if r.get('flagged'):
            flagged_reviews_with_conflict.add(rid)
        else:
            all_reviews_with_conflict_nonflagged.add(rid)
result = {
    'n_distinct_reviews_with_any_status_conflict_read': len(conflict_reviews),
    'n_flagged_reviews_total': len(flagged_reviews),
    'n_flagged_reviews_with_at_least_one_conflict_read': len(flagged_reviews_with_conflict),
    'n_flagged_reviews_with_zero_conflict_reads': len(flagged_reviews) - len(flagged_reviews_with_conflict),
    'n_nonflagged_reviews_with_at_least_one_conflict_read': len(all_reviews_with_conflict_nonflagged),
}
print(json.dumps(result, indent=2))
with open('CLAUDE_EDGE_REVIEW_BREAKDOWN.json','w') as f:
    json.dump(result, f, indent=2)
