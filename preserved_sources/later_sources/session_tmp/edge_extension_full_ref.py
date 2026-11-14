from pathlib import Path
import sys, json, subprocess, collections, hashlib, datetime

B = Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912')
O = B/'locus_specific_validation_20260919'
R = B/'technical_followup_20260919/paired_reads'
M = B/'mapping_extension_20260919'
sys.path.insert(0, str(O/'vendor'))
sys.path.insert(0, str(M/'source'))
import pysam, read_qc_core as core

fa = pysam.FastaFile(str(M/'reference/Homo_sapiens_assembly38.fasta'))
P = 5998724
lo = P - 1001
ref = fa.fetch('chrX', lo, P+1000).upper()
assert ref[1000:1003] == 'AAG'
haps = {
    'reference_AAG': ref,
    'SNP_only_GAG': ref[:1000] + 'G' + ref[1001:],
    'MNV_GAT': ref[:1000] + 'GAT' + ref[1003:],
    'neighbor_only_AAT': ref[:1002] + 'T' + ref[1003:],
}
for h, s in haps.items():
    assert len(s) == len(ref)
rc = lambda s: s.translate(str.maketrans('ACGT','TGCA'))[::-1]

READ_LEN = 151
FRAG_LEN = 350
GAP = FRAG_LEN - READ_LEN

plan = dict(
    purpose='Extension of synthetic_haplotype_controls.py (Section F negative control) to real-read conditions: (1) read_length corrected to 151bp matching the real conflicting-read length signature in CONFLICT_READ_CONTEXT.json (was 150bp), (2) full cycle sweep 0-150 covering both mid and edge placement in one run, (3) alignment reuses the EXACT real production command lines (paired_minimap.command.json / paired_bwa.command.json) with only the fastq paths substituted, so both aligners run against the SAME full GRCh38 reference (minimap2, with alt contigs) and the same prebuilt full-genome BWA-MEM2 index used for the real paired-read alignments -- NOT the chrX-only reference used in the prior CLAUDE_EDGE_PLACEMENT_RESULTS.json run. This directly tests whether cross-chromosome/genome-wide context (not available to a chrX-only reference) changes the edge-placement misclassification pattern.',
    still_unaddressed=['quality=40 Phred uniform, zero simulated sequencing errors (no realistic per-cycle error/quality decay model)', 'idealized haplotypes only (4 candidates), not a full local-assembly/diplotype model', 'synthetic reads, not real reads -- this remains a controlled simulation, not a replacement for the real-read analysis'],
    read_length=READ_LEN, fragment_length=FRAG_LEN, cycles_swept='0-150 (full read, both edge and mid in one run)',
    orientations=['forward','reverse'], quality=40, sequencing_errors=0,
    haplotypes=list(haps.keys()),
    target_position_grch38='chrX:5998724', reference_note='full GRCh38 (same command line as real paired-read alignment)',
)
(O/'FULL_REF_EDGE_EXTENSION_PLAN.json').write_text(json.dumps(plan, indent=2))

truth = {}
f1 = (O/'fullref_synthetic_R1.fastq').open('w')
f2 = (O/'fullref_synthetic_R2.fastq').open('w')
for h, s in haps.items():
    for orient in ['forward','reverse']:
        for cycle in range(0, READ_LEN):
            target_read_start = 1000 - cycle if orient == 'forward' else 1000 - (READ_LEN - 1 - cycle)
            start = target_read_start if orient == 'forward' else target_read_start - GAP
            if start < 0 or start + FRAG_LEN > len(s):
                continue
            frag = s[start:start+FRAG_LEN]
            assert len(frag) == FRAG_LEN
            q = h + '_' + orient[0] + '_c' + str(cycle)
            truth[q] = dict(haplotype=h, expected='ALT' if s[1000] != ref[1000] or s[1002] != ref[1002] else 'REF', orientation=orient, cycle=cycle)
            r1 = frag[:READ_LEN]
            r2 = rc(frag[-READ_LEN:])
            f1.write('@'+q+'\n'+r1+'\n+\n'+('I'*READ_LEN)+'\n')
            f2.write('@'+q+'\n'+r2+'\n+\n'+('I'*READ_LEN)+'\n')
f1.close(); f2.close()
(O/'fullref_synthetic_truth.json').write_text(json.dumps(truth, indent=2))
print('N_TEMPLATES', len(truth), flush=True)
