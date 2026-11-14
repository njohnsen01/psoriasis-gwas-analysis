from pathlib import Path
import sys, json, subprocess, collections, hashlib, datetime, time

B = Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912')
O = B/'locus_specific_validation_20260919'
R = B/'technical_followup_20260919/paired_reads'
M = B/'mapping_extension_20260919'
sys.path.insert(0, str(O/'vendor'))
sys.path.insert(0, str(M/'source'))
import read_qc_core as core

truth = json.loads((O/'fullref_synthetic_truth.json').read_text())
summary_all = []
for model, cmdname in [('minimap2','paired_minimap.command.json'), ('bwa_mem2','paired_bwa.command.json')]:
    cmd = json.loads((R/cmdname).read_text())
    assert str(R/'R1.private.fastq') in cmd and str(R/'R2.private.fastq') in cmd, cmd
    cmd = [str(O/'fullref_synthetic_R1.fastq') if x == str(R/'R1.private.fastq') else str(O/'fullref_synthetic_R2.fastq') if x == str(R/'R2.private.fastq') else x for x in cmd]
    (O/(model+'_fullref_synthetic_command.json')).write_text(json.dumps(cmd, indent=2))
    t0 = time.time()
    print('START', model, flush=True)
    with (O/(model+'_fullref_synthetic.sam')).open('w') as out:
        subprocess.run(cmd, stdout=out, stderr=(O/(model+'_fullref_synthetic.log')).open('w'), check=True, timeout=1800)
    print('ALIGN_DONE', model, round(time.time()-t0,1), 's', flush=True)
    obs = collections.defaultdict(list)
    seen = set()
    for line in (O/(model+'_fullref_synthetic.sam')).open():
        if line.startswith('@'): continue
        a = line.rstrip('\n').split('\t')
        flag = int(a[1]); q = a[0]
        if flag & 2304: continue
        end = '1' if flag & 64 else '2'
        assert (q,end) not in seen, (q,end)
        seen.add((q,end))
        r = core.observe(a, core.TARGETS[2])
        if r: obs[q].append(r)
    missing = {q for q in truth for e in ('1','2')} - {q for (q,e) in seen}
    assert seen == {(q,e) for q in truth for e in ('1','2')}, ('query accounting mismatch', len(seen), len(truth)*2)
    details = []
    for q, t in truth.items():
        v, conf = core.unique_high_quality(obs[q])
        assert len(v) <= 1
        call = v[0]['call'] if v else 'NO_HQ_CALL'
        details.append(dict(q=q, **t, observed=call, correct=call==t['expected']))
    (O/(model+'_fullref_synthetic_results.json')).write_text(json.dumps(details, indent=2))
    for h in sorted(set(x['haplotype'] for x in details)):
        v = [x for x in details if x['haplotype']==h]
        counts = collections.Counter(x['observed'] for x in v)
        summary_all.append(dict(model=model, haplotype=h, n_templates=len(v), n_correct=sum(x['correct'] for x in v),
                                 REF=counts['REF'], ALT=counts['ALT'], OTHER=counts['OTHER'], NO_HQ_CALL=counts['NO_HQ_CALL']))
    # bucket by edge (cycle<=15 or cycle>=135) vs mid, matching prior CLAUDE_EDGE_PLACEMENT_RESULTS.json convention
    for bucket_name, pred in [('edge', lambda c: c <= 15 or c >= 135), ('mid', lambda c: 15 < c < 135)]:
        vv = [x for x in details if pred(x['cycle'])]
        summary_all.append(dict(model=model, haplotype='__BUCKET__'+bucket_name, n_templates=len(vv), n_correct=sum(x['correct'] for x in vv),
                                 REF=None, ALT=None, OTHER=None, NO_HQ_CALL=sum(1 for x in vv if x['observed']=='NO_HQ_CALL')))
    print('DONE', model, flush=True)
(O/'FULL_REF_EDGE_EXTENSION_RESULTS.json').write_text(json.dumps(summary_all, indent=2))
print(json.dumps([s for s in summary_all if s['haplotype'].startswith('__BUCKET__')], indent=2))
