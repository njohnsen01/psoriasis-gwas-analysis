# ============================================================
# Corrected targeted extraction with explicit effect-allele harmonization
# ============================================================
import os
from pathlib import Path
import shutil
import subprocess
import pandas as pd
import numpy as np

ATTEMPT_ID = f"{os.environ['PSO_PRIMARY_RUN_ID']}"
PROJECT = f"{os.environ['GOOGLE_PROJECT']}"
ROOT = Path("/home/jupyter/psoriasis_corrected_downstream_2026-07-31")
PRIVATE = ROOT / "private_participant_level"
AGG = ROOT / "aggregate_outputs"
SCRATCH = PRIVATE / "targeted_chunk_staging"
EXPORTS = PRIVATE / "corrected_missing_lead_exports"
SCRATCH.mkdir(parents=True, exist_ok=True)
EXPORTS.mkdir(parents=True, exist_ok=True)

workspace_bucket = os.environ["WORKSPACE_BUCKET"].removeprefix("gs://").strip("/")
base_uri = (
    f"gs://{workspace_bucket}/"
    "psoriasis_regenie_loco_indexed_range_v1/"
    f"production_{ATTEMPT_ID}"
)
jobs = [
    ("5",49,["chr5:131277907:T:C","chr5:131381397:C:T"]),
    ("5",50,["chr5:131726714:T:C","chr5:131875033:TATC:T","chr5:132004339:C:T"]),
    ("5",59,["chr5:159137047:G:A"]),
    ("X",22,["chrX:67736589:T:C","chrX:67842269:G:T"]),
]
leads = pd.read_csv(
    f"/home/jupyter/psoriasis_regenie_loco_step2/attempt_{ATTEMPT_ID}/final_aggregate/psoriasis_corrected_lead_centered_100kb_groups.tsv",
    sep="\t",
).set_index("ID")

def run(cmd,label):
    cp=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if cp.returncode!=0:
        raise RuntimeError(label+" failed:\n"+cp.stdout[-5000:])
    return cp.stdout

def pvar_records(path, ids):
    text=run(["zstdcat",str(path)],"Read PVAR")
    found={}
    for line in text.splitlines():
        f=line.split("\t")
        if len(f)>=5 and f[2] in ids:
            found[f[2]]={"REF":f[3],"ALT":f[4]}
    if set(found)!=set(ids):
        raise RuntimeError(f"Missing PVAR records: {set(ids)-set(found)}")
    return found

raw_frames=[]
records=[]
for chrom,chunk,ids in jobs:
    tag=f"chr{chrom}_chunk_{chunk:04d}"
    stage=SCRATCH/tag
    stage.mkdir(parents=True,exist_ok=True)
    prefix=stage/f"chr{chrom}_pso_chunk_{chunk:04d}"
    root_uri=f"{base_uri}/chr{chrom}/secure_step2_chunks"
    uri_map={
        Path(str(prefix)+".pgen"):f"{root_uri}/chr{chrom}_pso_chunk_{chunk:04d}.pgen",
        Path(str(prefix)+".pvar.zst"):f"{root_uri}/chr{chrom}_pso_chunk_{chunk:04d}.pvar.zst",
        Path(str(prefix)+".psam"):f"{root_uri}/chr{chrom}_pso_step2.psam",
    }
    for dest,uri in uri_map.items():
        if not dest.exists() or dest.stat().st_size==0:
            run(["gcloud","storage","cp",f"--billing-project={PROJECT}",uri,str(dest)],f"Stage {dest.name}")
    allele_meta=pvar_records(Path(str(prefix)+".pvar.zst"),ids)
    extract_file=stage/f"{tag}_ids.txt"
    extract_file.write_text("\n".join(ids)+"\n")
    out=EXPORTS/f"{tag}_corrected_leads"
    log=run([
        "plink2","--pfile",str(prefix),"vzs",
        "--extract",str(extract_file),"--export","A",
        "--threads","24","--out",str(out)
    ],f"PLINK2 export {tag}")
    raw_path=Path(str(out)+".raw")
    frame=pd.read_csv(raw_path,sep=r"\s+")
    out_frame=frame[["FID","IID"]].copy()
    for vid in ids:
        cols=[c for c in frame.columns if c.startswith(vid+"_")]
        if len(cols)!=1:
            raise RuntimeError(f"{tag}: expected one column for {vid}; found {cols}")
        source_col=cols[0]
        counted=source_col[len(vid)+1:]
        effect=str(leads.loc[vid,"ALT"])
        ref=allele_meta[vid]["REF"]
        alt=allele_meta[vid]["ALT"]
        if effect!=alt:
            raise RuntimeError(f"{vid}: corrected effect allele {effect} != PVAR ALT {alt}")
        source=pd.to_numeric(frame[source_col],errors="coerce")
        male_haploid_scaled=False
        if counted==effect:
            effect_dosage=source.copy()
        elif counted==ref:
            if chrom=="X":
                male=frame["SEX"].eq(1)
                male_max=float(source.loc[male].max())
                if male_max<=1.01:
                    effect_dosage=2.0*(1.0-source)
                    effect_dosage.loc[~male]=2.0-source.loc[~male]
                    male_haploid_scaled=True
                else:
                    effect_dosage=2.0-source
            else:
                effect_dosage=2.0-source
        else:
            raise RuntimeError(f"{vid}: counted allele {counted} is neither REF nor effect ALT")
        target_col=vid+"_"+effect
        out_frame[target_col]=effect_dosage
        observed_af=float(effect_dosage.mean()/2.0)
        expected_af=float(leads.loc[vid,"A1FREQ"])
        af_difference=observed_af-expected_af
        if abs(af_difference)>0.01:
            raise RuntimeError(
                f"{vid}: effect dosage AF {observed_af:.5f} differs from GWAS {expected_af:.5f}"
            )
        records.append({
            "variant_id":vid,"chromosome":chrom,"chunk":chunk,
            "pvar_ref":ref,"pvar_alt":alt,
            "plink_export_column":source_col,
            "plink_counted_allele":counted,
            "corrected_effect_allele":effect,
            "dosage_transformed":bool(counted!=effect),
            "male_haploid_dosage_rescaled_to_0_2":male_haploid_scaled,
            "observed_effect_allele_frequency":observed_af,
            "gwas_effect_allele_frequency":expected_af,
            "frequency_difference":af_difference,
            "n_rows":len(frame),
            "plink_completed":("End time:" in log),
        })
    raw_frames.append(out_frame)

combined=raw_frames[0]
for frame in raw_frames[1:]:
    combined=combined.merge(frame,on=["FID","IID"],how="inner",validate="one_to_one")
if len(combined)!=392155:
    raise RuntimeError(f"Unexpected participant count: {len(combined)}")
dose_cols=[c for c in combined.columns if c not in {"FID","IID"}]
if len(dose_cols)!=8:
    raise RuntimeError(f"Expected 8 harmonized dosages; found {len(dose_cols)}")

combined_path=PRIVATE/"corrected_8_missing_lead_dosages_effect_allele.tsv.gz"
combined.to_csv(combined_path,sep="\t",index=False,compression="gzip")
manifest=pd.DataFrame(records)
manifest_path=AGG/"corrected_8_missing_lead_extraction_manifest.tsv"
manifest.to_csv(manifest_path,sep="\t",index=False)

qc=[]
for col in dose_cols:
    v=pd.to_numeric(combined[col],errors="coerce")
    qc.append({
        "variant_id":col.rsplit("_",1)[0],
        "counted_effect_allele":col.rsplit("_",1)[1],
        "n_total":len(v),"n_nonmissing":int(v.notna().sum()),
        "n_missing":int(v.isna().sum()),"missing_rate":float(v.isna().mean()),
        "dosage_mean":float(v.mean()),"dosage_sd":float(v.std(ddof=1)),
        "dosage_min":float(v.min()),"dosage_max":float(v.max()),
    })
qc=pd.DataFrame(qc)
qc_path=AGG/"corrected_8_missing_lead_dosage_qc.tsv"
qc.to_csv(qc_path,sep="\t",index=False)

shutil.rmtree(SCRATCH)

print("="*72)
print("TARGETED CORRECTED-LEAD DOSAGE EXTRACTION COMPLETE")
print("="*72)
print("Variants extracted and effect-allele harmonized:",len(dose_cols))
print("Participant rows:",len(combined))
print("Allele-frequency validation passed:",True)
print("Temporary genotype staging removed:",not SCRATCH.exists())
print("Private dosage file:",combined_path)
print("Aggregate manifest:",manifest_path)
print("Aggregate QC:",qc_path)
print("\nAllele harmonization summary:")
print(manifest[["variant_id","plink_counted_allele","corrected_effect_allele","dosage_transformed","male_haploid_dosage_rescaled_to_0_2","frequency_difference"]].to_string(index=False))
print("\nDosage QC summary:")
print(qc[["variant_id","n_nonmissing","missing_rate","dosage_min","dosage_max"]].to_string(index=False))
print("Participant identifiers displayed: False")
print("Participant-level files downloaded: False")
