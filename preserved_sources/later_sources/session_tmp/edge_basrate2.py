import json, re

TARGET_POS = 5998724
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
INFORMATIVE = ('HQ_REF', 'HQ_ALT')

def query_offset_of_ref_pos(pos, cigar, target):
    ops = CIGAR_RE.findall(cigar)
    ref_pos = pos
    query_pos = 0
    query_len = 0
    found = None
    has_clip = any(op in ('S', 'H') for _, op in ops)
    for length, op in ops:
        length = int(length)
        if op in ('M', '=', 'X'):
            if found is None and ref_pos <= target < ref_pos + length:
                found = query_pos + (target - ref_pos)
            ref_pos += length; query_pos += length; query_len += length
        elif op in ('D', 'N'):
            ref_pos += length
        elif op in ('I', 'S'):
            query_pos += length; query_len += length
    return found, query_len, has_clip

rows = json.loads(open('read_transitions.private.json').read())

def edge_for(side_dict):
    cigar = side_dict.get('CIGAR'); pos = side_dict.get('pos')
    if not cigar or pos is None:
        return None
    off, qlen, has_clip = query_offset_of_ref_pos(pos, cigar, TARGET_POS)
    if off is None or qlen == 0:
        return None
    return (min(off, qlen - 1 - off), has_clip, qlen)

# Set B: literal per-read BWA-vs-minimap status mismatch among informative (HQ_REF/HQ_ALT) calls on both sides
conflict_edge = []
concordant_edge = []
n_conflict_rows = 0
n_concordant_rows = 0
n_review_flagged = 0
n_review_flagged_and_status_conflict = 0

for r in rows:
    bwa = r.get('bwa'); mm = r.get('minimap')
    if r.get('flagged'):
        n_review_flagged += 1
    if not bwa or not mm:
        continue
    bstat = bwa.get('status'); mstat = mm.get('status')
    if bstat not in INFORMATIVE or mstat not in INFORMATIVE:
        continue
    e = edge_for(bwa)
    if e is None:
        continue
    if bstat != mstat:
        n_conflict_rows += 1
        conflict_edge.append(e)
        if r.get('flagged'):
            n_review_flagged_and_status_conflict += 1
    else:
        n_concordant_rows += 1
        concordant_edge.append(e)

def summarize(lst):
    n = len(lst)
    if n == 0:
        return {"n": 0}
    within15 = sum(1 for d, c, q in lst if d <= 15)
    within10 = sum(1 for d, c, q in lst if d <= 10)
    return {"n": n, "n_within_15bp": within15, "frac_within_15bp": round(within15/n, 4), "n_within_10bp": within10}

result = {
    "purpose": "Literal per-read BWA-vs-minimap2 status conflict (bstat != mstat among HQ_REF/HQ_ALT calls) at rs183018930, edge-distance base rate vs concordant reads. Cross-checked against review-level 'flagged' field used by REALIGNMENT_HAPLOTYPE_AGGREGATE.json.",
    "status_conflict_reads": summarize(conflict_edge),
    "status_concordant_reads": summarize(concordant_edge),
    "n_review_level_flagged": n_review_flagged,
    "n_review_flagged_that_are_also_status_conflict_reads": n_review_flagged_and_status_conflict,
}
print(json.dumps(result, indent=2))
with open('CLAUDE_EDGE_BASE_RATE_V2.json', 'w') as f:
    json.dump(result, f, indent=2)
