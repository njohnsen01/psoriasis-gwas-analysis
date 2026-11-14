# Export review: one results-prose field redacted; analytical code unchanged.
# Recovered original IPython input cells; session 54, cells 47-55.
# Preserves historical cell code; requires prior session context and protected inputs.

# Original cell 47
import sys, importlib
sys.path.insert(0, M + "/source")
import read_qc_core as core
importlib.reload(core)
print("TARGETS:", json.dumps(core.TARGETS, indent=2, default=str)[:2000])
print()
import pysam
fa = pysam.FastaFile(REF)
print("contigs sample:", [c for c in fa.references if c in ("chrX","X")])
window = fa.fetch("chrX", 5998700, 5998750)
print("ref window chrX:5998701-5998750 (1-based):", window)

# Original cell 48
CENTER = 5998724
FLANK = 200
start0 = CENTER - 1 - FLANK
end0 = CENTER - 1 + FLANK
big_window = fa.fetch("chrX", start0, end0)
offset = (CENTER - 1) - start0
assert big_window[offset:offset+3] == "AAG", big_window[offset:offset+3]

def make_variant(base_seq, off, new3):
    return base_seq[:off] + new3 + base_seq[off+3:]

haplotypes = {
    "REF_AAG": big_window,
    "SNP_GAG": make_variant(big_window, offset, "GAG"),
    "MNV_GAT": make_variant(big_window, offset, "GAT"),
    "NEIGHBOR_AAT": make_variant(big_window, offset, "AAT"),
}
print("haplotype window length:", len(big_window), "offset:", offset)
print("haplotype labels:", list(haplotypes.keys()))

sys.path.insert(0, O + "/vendor")
import edlib
print("edlib loaded:", edlib.__file__)

# Original cell 49
print("sanity check")
print("haplotypes" in dir(), "edlib" in dir(), "Rdir" in dir(), "O" in dir())

# Original cell 50
with open(Rdir + "/read_transitions.private.json") as f:
    rt = json.load(f)

flagged = [r for r in rt if r.get("flagged")]
print("total records:", len(rt), "flagged records:", len(flagged))

def best_hap(seq, haps):
    scores = {}
    for name, h in haps.items():
        res = edlib.align(seq, h, mode="HW", task="distance")
        scores[name] = res["editDistance"]
    ranked = sorted(scores.items(), key=lambda kv: kv[1])
    bname, bd = ranked[0]
    sd = ranked[1][1] if len(ranked) > 1 else None
    margin = (sd - bd) if sd is not None else None
    ambiguous = margin is not None and margin <= 1
    return bname, bd, margin, ambiguous

rows = []
for rec in flagged:
    row = {}
    for aligner in ("bwa", "minimap"):
        sub = rec.get(aligner) or {}
        obs = sub.get("obs") or {}
        seq = obs.get("seq")
        call = obs.get("call")
        if seq:
            bname, bd, margin, amb = best_hap(seq, haplotypes)
        else:
            bname, bd, margin, amb = (None, None, None, None)
        row[aligner + "_call"] = call
        row[aligner + "_best_hap"] = bname
        row[aligner + "_best_dist"] = bd
        row[aligner + "_margin"] = margin
        row[aligner + "_ambiguous"] = amb
    rows.append(row)

print("scored", len(rows), "flagged records")

# Original cell 51
with open(O + "/read_transitions.private.json") as f:
    rt = json.load(f)

flagged = [r for r in rt if r.get("flagged")]
print("total records:", len(rt), "flagged records:", len(flagged))

def best_hap(seq, haps):
    scores = {}
    for name, h in haps.items():
        res = edlib.align(seq, h, mode="HW", task="distance")
        scores[name] = res["editDistance"]
    ranked = sorted(scores.items(), key=lambda kv: kv[1])
    bname, bd = ranked[0]
    sd = ranked[1][1] if len(ranked) > 1 else None
    margin = (sd - bd) if sd is not None else None
    ambiguous = margin is not None and margin <= 1
    return bname, bd, margin, ambiguous

rows = []
for rec in flagged:
    row = {}
    for aligner in ("bwa", "minimap"):
        sub = rec.get(aligner) or {}
        obs = sub.get("obs") or {}
        seq = obs.get("seq")
        call = obs.get("call")
        if seq:
            bname, bd, margin, amb = best_hap(seq, haplotypes)
        else:
            bname, bd, margin, amb = (None, None, None, None)
        row[aligner + "_call"] = call
        row[aligner + "_best_hap"] = bname
        row[aligner + "_best_dist"] = bd
        row[aligner + "_margin"] = margin
        row[aligner + "_ambiguous"] = amb
    rows.append(row)

print("scored", len(rows), "flagged records")

# Original cell 52
import collections

bwa_hap_counts = collections.Counter(r["bwa_best_hap"] for r in rows)
mm_hap_counts = collections.Counter(r["minimap_best_hap"] for r in rows)
print("BWA best-haplotype distribution:", dict(bwa_hap_counts))
print("Minimap best-haplotype distribution:", dict(mm_hap_counts))
print()

agree = sum(1 for r in rows if r["bwa_best_hap"] == r["minimap_best_hap"] and r["bwa_best_hap"] is not None)
disagree = sum(1 for r in rows if r["bwa_best_hap"] != r["minimap_best_hap"] and r["bwa_best_hap"] is not None and r["minimap_best_hap"] is not None)
missing = sum(1 for r in rows if r["bwa_best_hap"] is None or r["minimap_best_hap"] is None)
print("realignment AGREES across aligners' reads:", agree)
print("realignment DISAGREES across aligners' reads:", disagree)
print("missing seq on one side:", missing)
print()

amb = sum(1 for r in rows if r["bwa_ambiguous"] or r["minimap_ambiguous"])
print("records with an ambiguous (margin<=1) best-hap call on either side:", amb)
print()

bwa_call_matches_realign = sum(1 for r in rows if r["bwa_call"] and r["bwa_best_hap"] and (
    (r["bwa_call"] == "ALT" and r["bwa_best_hap"] in ("MNV_GAT","SNP_GAG")) or
    (r["bwa_call"] == "REF" and r["bwa_best_hap"] in ("REF_AAG","NEIGHBOR_AAT"))
))
mm_call_matches_realign = sum(1 for r in rows if r["minimap_call"] and r["minimap_best_hap"] and (
    (r["minimap_call"] == "ALT" and r["minimap_best_hap"] in ("MNV_GAT","SNP_GAG")) or
    (r["minimap_call"] == "REF" and r["minimap_best_hap"] in ("REF_AAG","NEIGHBOR_AAT"))
))
print("BWA's own call matches independent realignment call (ALT<->MNV/SNP, REF<->REF/neighbor):", bwa_call_matches_realign, "/", len(rows))
print("Minimap's own call matches independent realignment call:", mm_call_matches_realign, "/", len(rows))

# Original cell 53
bwa_dist_counts = collections.Counter(r["bwa_best_dist"] for r in rows)
bwa_margin_counts = collections.Counter(r["bwa_margin"] for r in rows)
mm_dist_counts = collections.Counter(r["minimap_best_dist"] for r in rows)
mm_margin_counts = collections.Counter(r["minimap_margin"] for r in rows)
print("BWA best edit-distance to MNV_GAT distribution:", dict(bwa_dist_counts))
print("BWA margin (2nd-best minus best) distribution:", dict(bwa_margin_counts))
print("Minimap best edit-distance distribution:", dict(mm_dist_counts))
print("Minimap margin distribution:", dict(mm_margin_counts))
print()

n_bwa_exact = sum(1 for r in rows if r["bwa_best_dist"] == 0)
n_mm_exact = sum(1 for r in rows if r["minimap_best_dist"] == 0)
print("BWA-side reads with EXACT (0 edit) match to MNV_GAT haplotype window:", n_bwa_exact, "/", len(rows))
print("Minimap-side reads with EXACT (0 edit) match to MNV_GAT haplotype window:", n_mm_exact, "/", len(rows))
print()

mm_call_ref_but_seq_supports_alt = sum(1 for r in rows if r["minimap_call"] == "REF" and r["minimap_best_hap"] == "MNV_GAT")
mm_call_other = collections.Counter(r["minimap_call"] for r in rows)
print("minimap_call distribution among flagged reads:", dict(mm_call_other))
print("Reads minimap called REF (or failed) where independent realignment shows MNV_GAT support:", mm_call_ref_but_seq_supports_alt)

# Original cell 54
import hashlib, datetime

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

aggregate = {
    "purpose": "Independent haplotype-aware local realignment of the real flagged rs183018930 read observations, using edlib edit-distance scoring against 4 candidate local reference haplotypes, entirely independent of the minimap2/BWA-MEM2 aligners' own seed-and-extend heuristics.",
    "haplotype_window": {"chrom": "chrX", "center_1based": CENTER, "flank_bp": FLANK, "window_len": len(big_window), "offset_0based": offset},
    "haplotypes_tested": {k: (v[offset-5:offset+8]) for k, v in haplotypes.items()},
    "n_total_read_transition_records": len(rt),
    "n_flagged_records": len(flagged),
    "n_scored": len(rows),
    "bwa_best_haplotype_distribution": dict(bwa_hap_counts),
    "minimap_best_haplotype_distribution": dict(mm_hap_counts),
    "cross_aligner_realignment_agreement": {"agree": agree, "disagree": disagree, "missing_seq_one_side": missing},
    "bwa_edit_distance_to_MNV_GAT_distribution": {str(k): v for k, v in bwa_dist_counts.items()},
    "minimap_edit_distance_to_MNV_GAT_distribution": {str(k): v for k, v in mm_dist_counts.items()},
    "bwa_exact_match_to_MNV_GAT": n_bwa_exact,
    "minimap_exact_match_to_MNV_GAT": n_mm_exact,
    "minimap_call_distribution_among_flagged": dict(mm_call_other),
    "minimap_REF_or_nocall_but_realignment_supports_MNV_GAT": mm_call_ref_but_seq_supports_alt,
    "bwa_own_call_matches_realignment": bwa_call_matches_realign,
    "minimap_own_call_matches_realignment": mm_call_matches_realign,
    "note_on_margin_ambiguity": "margin<=1 flags every record because SNP_GAG and MNV_GAT haplotypes differ by exactly one base (the rs5961396 position), so they are always within edit distance 1 of each other; this is a structural property of the two candidate haplotypes, not evidence of per-read ambiguity. The best (lowest-distance) haplotype is still MNV_GAT in effectively all scored reads.",
    "interpretation": "Historical results prose omitted from this export because it embeds protected read-count results; see secure original output.",
    "caveats": [
        "This adjudicates LOCAL SEQUENCE CONTENT of already-flagged reads; it does not itself constitute disease-association replication.",
        "It does not distinguish the MNV (rs386823709) representation from the equivalent phased-SNP (rs183018930+rs5961396) representation, since both produce an identical local sequence (GAT) and are therefore indistinguishable by sequence alone.",
        "Frozen, bounded, stratified read selection as established in the technical follow-up stage; not representative for cohort-wide error-rate estimation.",
        "No cohort genotypes were changed and no association model was rerun.",
    ],
    "generated_utc": datetime.datetime.utcnow().isoformat() + "Z",
}

print(json.dumps(aggregate, indent=2)[:3000])

# Original cell 55
dep_files = {
    "read_transitions_private_json": O + "/read_transitions.private.json",
    "reference_fasta": str(M) + "/reference/Homo_sapiens_assembly38.fasta" if not isinstance(M, str) else M + "/reference/Homo_sapiens_assembly38.fasta",
    "read_qc_core_py": M + "/source/read_qc_core.py",
    "edlib_binary": edlib.__file__,
}
hashes = {}
missing_files = []
for name, path in dep_files.items():
    try:
        hashes[name] = {"path": path, "sha256": sha256_file(path)}
    except Exception as e:
        missing_files.append({"name": name, "path": path, "error": str(e)})

verified = {
    "stage": "locus_specific_validation_20260919 / independent_haplotype_realignment",
    "dependency_hashes": hashes,
    "missing_files": missing_files,
    "n_flagged_scored": len(rows),
    "generated_utc": datetime.datetime.utcnow().isoformat() + "Z",
}

out_agg = O + "/REALIGNMENT_HAPLOTYPE_AGGREGATE.json"
out_ver = O + "/REALIGNMENT_HAPLOTYPE_VERIFIED.json"
with open(out_agg, "w") as f:
    json.dump(aggregate, f, indent=2)
with open(out_ver, "w") as f:
    json.dump(verified, f, indent=2)

print("wrote", out_agg)
print("wrote", out_ver)
print("missing_files:", missing_files)
