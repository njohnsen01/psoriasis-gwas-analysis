#!/usr/bin/env python3
# PSORIASIS_DISK_BOUNDED_PATCH_2026_07_26_V2
# PSORIASIS_DISK_BOUNDED_PATCH_2026_07_26_V1
# PSORIASIS_FAST_PARALLEL_PATCH_2026_07_26_V1
"""Frozen corrected psoriasis REGENIE/LOCO production runner.

This runner is generated from the validated R1-R5 analysis implementation,
then wrapped in an immutable, restart-safe coordinator. It must be launched
through the matching Prepare and Launch/Resume notebook cells.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_REVIEW_FILES = {
    "corrected_primary_psoriasis_GWAS_summary.xlsx",
    "corrected_primary_analysis_summary.json",
    "corrected_step1_model_validation.json",
    "corrected_primary_chromosome_summary.tsv",
    "corrected_primary_frozen_lead_harmonization.tsv",
    "corrected_chrX_step2_coding_qc.json",
    "corrected_primary_genome_wide_significant_variants.tsv",
    "corrected_primary_plotting_data.tsv.gz",
    "aggregate_privacy_validation.json",
    "aggregate_manifest.tsv",
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def production_sha256(path, chunk_size=16 * 1024**2):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def production_config_sha256(payload):
    fingerprint_payload = dict(payload)
    fingerprint_payload.pop("config_sha256", None)
    serialized = json.dumps(
        fingerprint_payload,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def production_md5_base64(path, chunk_size=16 * 1024**2):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return base64.b64encode(digest.digest()).decode("ascii")


def production_atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True))
    os.replace(temporary, path)


if len(sys.argv) == 2 and sys.argv[1] == "--self-check":
    print("Runner outer syntax loaded.")
    raise SystemExit(0)
if len(sys.argv) != 2:
    raise SystemExit("Usage: runner.py PRODUCTION_CONFIG.json")

CONFIG_PATH = Path(sys.argv[1]).resolve()
PRODUCTION_CONFIG = json.loads(CONFIG_PATH.read_text())
if (
    production_config_sha256(PRODUCTION_CONFIG)
    != PRODUCTION_CONFIG["config_sha256"]
):
    raise RuntimeError("The immutable production config fingerprint changed.")
ATTEMPT_ROOT = Path(PRODUCTION_CONFIG["attempt_root"]).resolve()
STATUS_ROOT = Path(PRODUCTION_CONFIG["status_root"]).resolve()
STATE_PATH = Path(PRODUCTION_CONFIG["state_path"]).resolve()
SUMMARY_PATH = Path(PRODUCTION_CONFIG["summary_path"]).resolve()
RECEIPT_PATH = Path(PRODUCTION_CONFIG["receipt_path"]).resolve()
PID_PATH = Path(PRODUCTION_CONFIG["pid_path"]).resolve()
LOCK_PATH = Path(PRODUCTION_CONFIG["lock_path"]).resolve()
MASTER_LOG_PATH = Path(PRODUCTION_CONFIG["master_log_path"]).resolve()
RUNTIME_FLOOR_BYTES = int(
    float(PRODUCTION_CONFIG["runtime_floor_gib"]) * 1024**3
)
LAUNCH_TOKEN = os.environ.get("PSORIASIS_PRODUCTION_LAUNCH_TOKEN")

if Path(__file__).resolve() != Path(
    PRODUCTION_CONFIG["runner_path"]
).resolve():
    raise RuntimeError("The launched runner path differs from the frozen config.")
if production_sha256(__file__) != PRODUCTION_CONFIG["runner_sha256"]:
    raise RuntimeError("The frozen runner bytes changed after preparation.")
for tool_name in ["regenie", "plink2"]:
    expected_tool = PRODUCTION_CONFIG["software"][tool_name]
    observed_path_text = shutil.which(tool_name)
    if not observed_path_text:
        raise RuntimeError(f"Frozen tool is unavailable: {tool_name}")
    observed_path = Path(observed_path_text).resolve()
    if observed_path != Path(expected_tool["path"]).resolve():
        raise RuntimeError(f"Frozen {tool_name} executable path changed.")
    if production_sha256(observed_path) != expected_tool["sha256"]:
        raise RuntimeError(f"Frozen {tool_name} executable bytes changed.")
if os.environ.get("GOOGLE_PROJECT") != PRODUCTION_CONFIG["google_project"]:
    raise RuntimeError("GOOGLE_PROJECT differs from the immutable config.")
if not LAUNCH_TOKEN:
    raise RuntimeError("The production launch token is absent.")
owner_path = LOCK_PATH / "owner.json"
if not owner_path.exists():
    raise RuntimeError("The production launch reservation is absent.")
owner = json.loads(owner_path.read_text())
if (
    owner.get("attempt_id") != PRODUCTION_CONFIG["attempt_id"]
    or owner.get("launch_token") != LAUNCH_TOKEN
):
    raise RuntimeError("The production launch reservation is not owned here.")

PRODUCTION_STATE = (
    json.loads(STATE_PATH.read_text())
    if STATE_PATH.exists()
    else {}
)


def production_state_update(**changes):
    PRODUCTION_STATE.update(changes)
    PRODUCTION_STATE["attempt_id"] = PRODUCTION_CONFIG["attempt_id"]
    PRODUCTION_STATE["updated_utc"] = utc_now()
    disk = shutil.disk_usage("/home/jupyter")
    PRODUCTION_STATE["disk_free_gib"] = round(
        disk.free / 1024**3,
        2,
    )
    production_atomic_json(STATE_PATH, PRODUCTION_STATE)


def production_assert_disk_floor(label):
    free_bytes = shutil.disk_usage("/home/jupyter").free
    if free_bytes < RUNTIME_FLOOR_BYTES:
        raise RuntimeError(
            f"Disk safety floor breached {label}: "
            f"{free_bytes / 1024**3:.1f} GiB free; "
            f"{RUNTIME_FLOOR_BYTES / 1024**3:.1f} GiB required."
        )


def production_release_lock():
    if not LOCK_PATH.exists():
        return
    current_owner_path = LOCK_PATH / "owner.json"
    if not current_owner_path.exists():
        raise RuntimeError("Refusing to remove a lock without owner metadata.")
    current_owner = json.loads(current_owner_path.read_text())
    if (
        current_owner.get("attempt_id") != PRODUCTION_CONFIG["attempt_id"]
        or current_owner.get("launch_token") != LAUNCH_TOKEN
    ):
        raise RuntimeError("Refusing to remove another launcher's lock.")
    current_owner_path.unlink()
    LOCK_PATH.rmdir()


def production_signal_handler(signum, _frame):
    production_state_update(
        status="failed",
        stage=PRODUCTION_STATE.get("stage"),
        error=f"Production coordinator received signal {signum}.",
    )
    raise KeyboardInterrupt


signal.signal(signal.SIGTERM, production_signal_handler)
signal.signal(signal.SIGINT, production_signal_handler)
PID_PATH.write_text(str(os.getpid()))


STAGE_SOURCES = [
r'''
# ============================================================
# Block R1 — Initialize the corrected run and freeze inputs
# ============================================================

from pathlib import Path
import concurrent.futures
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time

import numpy as np
import pandas as pd

RUN_DATE = PRODUCTION_CONFIG["run_date"]
RUN_LABEL = PRODUCTION_CONFIG["attempt_id"]
RUN_ROOT = Path(PRODUCTION_CONFIG["attempt_root"])
PRIVATE = RUN_ROOT / "private_participant_level"
INPUTS = PRIVATE / "frozen_inputs"
QC_DIR = PRIVATE / "qc_genotypes"
STEP1_PARTS = PRIVATE / "step1_parts"
STEP1_DIR = PRIVATE / "step1"
STEP2_DIR = PRIVATE / "step2"
LOG_DIR = RUN_ROOT / "logs"
AGG_DIR = RUN_ROOT / "aggregate_outputs"
REPRO_DIR = RUN_ROOT / "secure_reproducibility"

for directory in [
    PRIVATE,
    INPUTS,
    QC_DIR,
    STEP1_PARTS,
    STEP1_DIR,
    STEP2_DIR,
    LOG_DIR,
    AGG_DIR,
    REPRO_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

GOOGLE_PROJECT = os.environ.get("GOOGLE_PROJECT")
if not GOOGLE_PROJECT:
    raise EnvironmentError("GOOGLE_PROJECT is not defined.")

CONTROLLED_V8 = PRODUCTION_CONFIG["controlled_v8_prefix"]


def run_checked(command, description, log_path=None):
    command = [str(value) for value in command]
    print(f"\n[{description}]")
    print(" ".join(command))

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )
    combined = (completed.stdout or "") + (completed.stderr or "")

    if log_path is not None:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(log_path).write_text(combined)

    print(combined)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{description} failed with exit code "
            f"{completed.returncode}."
        )
    return combined


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def stable_fingerprint(payload):
    serialized = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def atomic_write_json(path, payload):
    path = Path(path)
    temporary = Path(str(path) + ".partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True))
    os.replace(temporary, path)


def file_metadata(path, include_sha256=False):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty file: {path}")
    metadata = {
        "path": str(path),
        "bytes": int(path.stat().st_size),
        "mtime_ns": int(path.stat().st_mtime_ns),
    }
    if include_sha256:
        metadata["sha256"] = sha256_file(path)
    return metadata


def sentinel_matches(path, expected_fingerprint):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return False
    return payload.get("fingerprint") == expected_fingerprint


def recorded_file_matches(path, recorded, include_sha256=False):
    if not recorded:
        return False
    try:
        current = file_metadata(path, include_sha256=include_sha256)
    except Exception:
        return False
    keys = ["bytes", "mtime_ns"]
    if include_sha256:
        keys.append("sha256")
    return all(current.get(key) == recorded.get(key) for key in keys)


def bucket_exists(bucket):
    result = subprocess.run(
        [
            "gsutil",
            "-u",
            GOOGLE_PROJECT,
            "ls",
            "-b",
            bucket,
        ],
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


WORKSPACE_BUCKET = PRODUCTION_CONFIG["workspace_bucket"].rstrip("/")
if not bucket_exists(WORKSPACE_BUCKET):
    raise RuntimeError(
        "The immutable workspace bucket is not accessible: "
        f"{WORKSPACE_BUCKET}"
    )
print("Using immutable secure workspace bucket:", WORKSPACE_BUCKET)

REGENIE_VERSION = run_checked(
    ["regenie", "--version"],
    "Read REGENIE version",
    LOG_DIR / "regenie_version.txt",
).strip()

PLINK2_VERSION = run_checked(
    ["plink2", "--version"],
    "Read PLINK2 version",
    LOG_DIR / "plink2_version.txt",
).splitlines()[0].strip()


input_provenance = {}


def stage_analysis_file(filename):
    destination = INPUTS / filename
    if destination.exists() and destination.stat().st_size > 0:
        print("[FOUND FROZEN]", destination)
        input_provenance[filename] = {
            "source_type": "existing_frozen_run_input",
            "source": str(destination),
        }
        return destination

    local_candidates = [
        Path("/home/jupyter") / filename,
        Path("/home/jupyter/regenie_chr1_run") / filename,
        Path(
            "/home/jupyter/psoriasis_phase2_reproducibility/"
            "secure_participant_level"
        ) / filename,
    ]

    for source in local_candidates:
        if source.exists() and source.stat().st_size > 0:
            shutil.copy2(source, destination)
            print("[COPIED LOCAL]", source, "->", destination)
            input_provenance[filename] = {
                "source_type": "local_workspace_file",
                "source": str(source),
            }
            return destination

    run_checked(
        [
            "gcloud",
            "storage",
            "cp",
            f"--billing-project={GOOGLE_PROJECT}",
            f"{WORKSPACE_BUCKET}/{filename}",
            str(destination),
        ],
        f"Stage {filename} from the secure workspace bucket",
        LOG_DIR / f"stage_{filename}.log",
    )
    input_provenance[filename] = {
        "source_type": "secure_workspace_bucket",
        "source": f"{WORKSPACE_BUCKET}/{filename}",
    }
    return destination


PHENO_FILE = stage_analysis_file("psoriasis_pheno.tsv")
COVAR_FILE = stage_analysis_file("psoriasis_covariates.tsv")

pheno = pd.read_csv(
    PHENO_FILE,
    sep="\t",
    dtype={"FID": str, "IID": str},
)
covar = pd.read_csv(
    COVAR_FILE,
    sep="\t",
    dtype={"FID": str, "IID": str},
)

required_pheno = {"FID", "IID", "psoriasis_case"}
required_covar_list = [
    "FID", "IID", "sex", "age", "age2", "age_sex",
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PC6", "PC7", "PC8", "PC9", "PC10",
]
required_covar = set(required_covar_list)

if not required_pheno.issubset(pheno.columns):
    raise ValueError(
        f"Phenotype columns are incomplete: {pheno.columns.tolist()}"
    )
if not required_covar.issubset(covar.columns):
    raise ValueError(
        f"Covariate columns are incomplete: {covar.columns.tolist()}"
    )
if pheno[["FID", "IID"]].duplicated().any():
    raise ValueError("Duplicate phenotype FID/IID pairs detected.")
if covar[["FID", "IID"]].duplicated().any():
    raise ValueError("Duplicate covariate FID/IID pairs detected.")
if not (pheno["FID"] == pheno["IID"]).all():
    raise ValueError("Phenotype FID and IID are not identical.")

pheno_values = pd.to_numeric(
    pheno["psoriasis_case"],
    errors="raise",
)
if not (
    np.isfinite(pheno_values.to_numpy(dtype=float)).all()
    and pheno_values.isin([0, 1]).all()
):
    raise ValueError(
        "psoriasis_case must contain exact numeric 0/1 values only."
    )
pheno["psoriasis_case"] = pheno_values.astype(int)

observed_counts = (
    pheno["psoriasis_case"]
    .value_counts()
    .sort_index()
    .to_dict()
)

if pheno.shape[0] != 392155:
    raise ValueError(
        f"Expected 392,155 phenotype rows; observed {pheno.shape[0]:,}."
    )
if observed_counts != {0: 386350, 1: 5805}:
    raise ValueError(
        "Expected 386,350 controls and 5,805 cases; "
        f"observed {observed_counts}."
    )

pheno_ids = set(map(tuple, pheno[["FID", "IID"]].to_numpy()))
covar_ids = set(map(tuple, covar[["FID", "IID"]].to_numpy()))
if pheno_ids != covar_ids:
    raise ValueError(
        "Phenotype and covariate FID/IID sets do not match exactly."
    )

pheno_index = pd.MultiIndex.from_frame(pheno[["FID", "IID"]])
covar_analysis = (
    covar
    .set_index(["FID", "IID"])
    .loc[pheno_index, required_covar_list[2:]]
    .reset_index()
)
if not covar_analysis[["FID", "IID"]].equals(
    pheno[["FID", "IID"]].reset_index(drop=True)
):
    raise ValueError(
        "Covariates could not be aligned exactly to phenotype ID order."
    )

for column in required_covar_list[2:]:
    covar_analysis[column] = pd.to_numeric(
        covar_analysis[column],
        errors="raise",
    )
    if not np.isfinite(
        covar_analysis[column].to_numpy(dtype=float)
    ).all():
        raise ValueError(
            f"Covariate {column} contains a non-finite value."
        )

if not covar_analysis["sex"].isin([0, 1]).all():
    raise ValueError(
        "Frozen sex covariate must contain exact numeric 0/1 "
        "values only."
    )
covar_analysis["sex"] = covar_analysis["sex"].astype(int)
if not np.allclose(
    covar_analysis["age2"].to_numpy(float),
    np.square(covar_analysis["age"].to_numpy(float)),
    rtol=1e-10,
    atol=1e-8,
):
    raise ValueError(
        "age2 does not reproduce the frozen age-squared covariate."
    )
if not np.allclose(
    covar_analysis["age_sex"].to_numpy(float),
    (
        covar_analysis["age"].to_numpy(float)
        * covar_analysis["sex"].to_numpy(float)
    ),
    rtol=1e-10,
    atol=1e-8,
):
    raise ValueError(
        "age_sex does not reproduce the frozen age-by-sex "
        "interaction covariate."
    )

FROZEN_PHENO = INPUTS / "psoriasis_pheno_frozen.tsv"
FROZEN_COVAR = INPUTS / "psoriasis_covariates_frozen.tsv"
KEEP_FILE = INPUTS / "psoriasis_analysis_keep.txt"

pheno[["FID", "IID", "psoriasis_case"]].to_csv(
    FROZEN_PHENO,
    sep="\t",
    index=False,
)
covar_analysis.to_csv(
    FROZEN_COVAR,
    sep="\t",
    index=False,
)
pheno[["FID", "IID"]].to_csv(
    KEEP_FILE,
    sep=" ",
    index=False,
    header=False,
)

input_hashes = {
    "phenotype_sha256": sha256_file(FROZEN_PHENO),
    "covariate_sha256": sha256_file(FROZEN_COVAR),
    "keep_sha256": sha256_file(KEEP_FILE),
}
frozen_id_order_sha256 = stable_fingerprint(
    pheno[["FID", "IID"]].astype(str).values.tolist()
)

setup_summary = {
    "run_date": RUN_DATE,
    "run_label": RUN_LABEL,
    "workspace_bucket": WORKSPACE_BUCKET,
    "controlled_v8_bgen_root": CONTROLLED_V8,
    "regenie_version": REGENIE_VERSION,
    "plink2_version": PLINK2_VERSION,
    "n_participants": int(pheno.shape[0]),
    "n_cases": int(observed_counts[1]),
    "n_controls": int(observed_counts[0]),
    "covariates": covar_analysis.columns[2:].tolist(),
    "frozen_input_hashes": input_hashes,
    "frozen_id_order_sha256": frozen_id_order_sha256,
    "input_provenance": input_provenance,
    "analysis_design": (
        "one autosomal chr1-22 Step 1 model; shared 23-row "
        "LOCO predictions; chromosome-specific Step 2 on chr1-22 and X"
    ),
    "historical_primary_regenie_version": "2.0.2.gz",
    "corrected_run_software_change_note": (
        "This rerun corrects the Step 1/LOCO design and uses the "
        "currently installed REGENIE version recorded above."
    ),
}
setup_summary["fingerprint"] = stable_fingerprint(setup_summary)
atomic_write_json(
    AGG_DIR / "corrected_run_setup.json",
    setup_summary,
)

print("\nFrozen input validation passed.")
print(json.dumps(setup_summary, indent=2))
''',
r'''
# ============================================================
# Block R2 — Locate retained QC BGENs; reconstruct only missing ones
# ============================================================

CHROMOSOMES = [str(chrom) for chrom in range(1, 23)] + ["X"]

# Safety default: do not unexpectedly re-download up to 23 controlled
# whole-chromosome BGENs. Block R2 will first report any missing retained
# QC pairs. Set this to True only after reviewing that missing list.
ALLOW_CONTROLLED_RECONSTRUCTION = bool(
    PRODUCTION_CONFIG["allow_controlled_reconstruction"]
)


def nonempty(path):
    path = Path(path)
    return path.exists() and path.stat().st_size > 0


def retained_qc_pair(chromosome):
    prefix = (
        Path(f"/home/jupyter/regenie_chr{chromosome}_run")
        / "plink"
        / f"chr{chromosome}_qc"
    )
    bgen = Path(str(prefix) + ".bgen")
    sample = Path(str(prefix) + ".sample")
    if nonempty(bgen) and nonempty(sample):
        return bgen, sample
    return None


def gcs_object_size(uri):
    completed = subprocess.run(
        [
            "gcloud",
            "storage",
            "objects",
            "describe",
            uri,
            f"--billing-project={GOOGLE_PROJECT}",
            "--format=value(size)",
        ],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Could not read controlled-object size for {uri}:\n"
            + (completed.stdout or "")
            + (completed.stderr or "")
        )
    return int(completed.stdout.strip())


def stage_controlled_object(uri, destination, log_name):
    destination = Path(destination)
    expected_bytes = gcs_object_size(uri)
    if (
        destination.exists()
        and destination.stat().st_size == expected_bytes
    ):
        return destination, expected_bytes

    partial = Path(str(destination) + ".partial")
    partial.unlink(missing_ok=True)
    run_checked(
        [
            "gcloud",
            "storage",
            "cp",
            f"--billing-project={GOOGLE_PROJECT}",
            uri,
            str(partial),
        ],
        f"Stage controlled object {Path(uri).name}",
        LOG_DIR / log_name,
    )
    if not partial.exists() or partial.stat().st_size != expected_bytes:
        raise RuntimeError(
            f"Controlled copy-size validation failed for {uri}. "
            f"Expected {expected_bytes:,} bytes."
        )
    os.replace(partial, destination)
    return destination, expected_bytes


qc_sources = {}
inventory_rows = []

retained_pairs = {
    chromosome: retained_qc_pair(chromosome)
    for chromosome in CHROMOSOMES
}
missing_retained_pairs = [
    chromosome
    for chromosome, pair in retained_pairs.items()
    if pair is None
]

if missing_retained_pairs and not ALLOW_CONTROLLED_RECONSTRUCTION:
    missing_inventory = pd.DataFrame(
        {
            "chromosome": missing_retained_pairs,
            "retained_qc_pair_found": False,
            "action": (
                "Review this list, then set "
                "ALLOW_CONTROLLED_RECONSTRUCTION=True "
                "and rerun Block R2."
            ),
        }
    )
    missing_inventory.to_csv(
        AGG_DIR / "missing_retained_qc_pairs.tsv",
        sep="\t",
        index=False,
    )
    display(missing_inventory)
    raise FileNotFoundError(
        "Retained QC BGEN/SAMPLE pairs are missing for: "
        + ", ".join(
            f"chr{chromosome}"
            for chromosome in missing_retained_pairs
        )
        + "\nNo controlled whole-chromosome files were copied. "
        "Review the list above. If reconstruction is necessary, "
        "set ALLOW_CONTROLLED_RECONSTRUCTION=True and rerun R2; "
        "only the listed chromosomes will be reconstructed."
    )

for chromosome in CHROMOSOMES:
    retained = retained_pairs[chromosome]
    if retained is not None:
        bgen_path, sample_path = retained
        status = "retained_original_qc_bgen"
        source_fingerprint_payload = {
            "chromosome": chromosome,
            "source_status": status,
            "bgen": file_metadata(
                bgen_path,
                include_sha256=True,
            ),
            "sample": file_metadata(
                sample_path,
                include_sha256=True,
            ),
        }
    else:
        new_prefix = QC_DIR / f"chr{chromosome}_qc"
        bgen_path = Path(str(new_prefix) + ".bgen")
        sample_path = Path(str(new_prefix) + ".sample")
        reconstruction_sentinel = Path(
            str(new_prefix) + ".complete.json"
        )
        controlled_bgen_uri = (
            f"{CONTROLLED_V8}/chr{chromosome}.bgen"
        )
        controlled_sample_uri = (
            f"{CONTROLLED_V8}/chr{chromosome}.sample"
        )
        reconstruction_spec = {
            "chromosome": chromosome,
            "controlled_bgen_uri": controlled_bgen_uri,
            "controlled_sample_uri": controlled_sample_uri,
            "plink2_version": PLINK2_VERSION,
            "qc_options": [
                "--geno 0.02",
                "--maf 0.01",
                "--mac 100",
                "--hwe 1e-15 keep-fewhet",
                "--set-missing-var-ids @:#$1,$2",
                "--new-id-max-allele-len 1000",
                "--make-bed",
                "--export bgen-1.2",
            ],
            "chrX_split_par": (
                "hg38" if chromosome == "X" else None
            ),
            "hard_call_bed_roundtrip": True,
        }
        reconstruction_fingerprint = stable_fingerprint(
            reconstruction_spec
        )
        reconstruction_ready = (
            nonempty(bgen_path)
            and nonempty(sample_path)
            and sentinel_matches(
                reconstruction_sentinel,
                reconstruction_fingerprint,
            )
        )

        if not reconstruction_ready:
            raw_dir = PRIVATE / "temporary_raw_bgen"
            raw_dir.mkdir(parents=True, exist_ok=True)
            raw_bgen = raw_dir / f"chr{chromosome}.bgen"
            raw_sample = raw_dir / f"chr{chromosome}.sample"

            raw_metadata = {}
            for uri, destination, label in [
                (controlled_bgen_uri, raw_bgen, "bgen"),
                (controlled_sample_uri, raw_sample, "sample"),
            ]:
                staged, expected_bytes = stage_controlled_object(
                    uri,
                    destination,
                    (
                        f"stage_controlled_chr{chromosome}_"
                        f"{label}.log"
                    ),
                )
                raw_metadata[label] = {
                    "uri": uri,
                    "bytes": expected_bytes,
                    "local_path": str(staged),
                }

            temporary_bed_prefix = (
                raw_dir / f"chr{chromosome}_qc_bed_roundtrip.partial"
            )
            partial_qc_prefix = (
                QC_DIR / f"chr{chromosome}_qc.partial"
            )
            for prefix, suffixes in [
                (
                    temporary_bed_prefix,
                    [".bed", ".bim", ".fam", ".log"],
                ),
                (
                    partial_qc_prefix,
                    [".bgen", ".sample", ".log"],
                ),
            ]:
                for suffix in suffixes:
                    Path(str(prefix) + suffix).unlink(missing_ok=True)

            make_bed_command = [
                "plink2",
                "--bgen", str(raw_bgen), "ref-first",
                "--sample", str(raw_sample),
            ]
            if chromosome == "X":
                make_bed_command += ["--split-par", "hg38"]
            make_bed_command += [
                "--geno", "0.02",
                "--maf", "0.01",
                "--mac", "100",
                "--hwe", "1e-15", "keep-fewhet",
                "--set-missing-var-ids", "@:#$1,$2",
                "--new-id-max-allele-len", "1000",
                "--make-bed",
                "--threads", "24",
                "--out", str(temporary_bed_prefix),
            ]

            run_checked(
                make_bed_command,
                (
                    "Apply original QC and create the historical "
                    f"hard-call BED representation for chr{chromosome}"
                ),
                LOG_DIR
                / f"reconstruct_chr{chromosome}_qc_make_bed.log",
            )

            # The historical QC BGENs were exported from BED, so this
            # round-trip is necessary to reproduce their hard-called
            # genotype representation. Direct raw-BGEN-to-BGEN export
            # would retain dosages and would define a different analysis.
            run_checked(
                [
                    "plink2",
                    "--bfile", str(temporary_bed_prefix),
                    "--export", "bgen-1.2",
                    "--threads", "24",
                    "--out", str(partial_qc_prefix),
                ],
                (
                    "Export retained-QC-equivalent hard-call BGEN for "
                    f"chr{chromosome}"
                ),
                LOG_DIR
                / f"reconstruct_chr{chromosome}_qc_export_bgen.log",
            )

            partial_bgen = Path(
                str(partial_qc_prefix) + ".bgen"
            )
            partial_sample = Path(
                str(partial_qc_prefix) + ".sample"
            )
            if not (
                nonempty(partial_bgen)
                and nonempty(partial_sample)
            ):
                raise RuntimeError(
                    f"QC reconstruction failed for chr{chromosome}."
                )
            os.replace(partial_bgen, bgen_path)
            os.replace(partial_sample, sample_path)

            # The controlled source is recoverable. Remove only the temporary
            # raw copies after the QC BGEN and SAMPLE have been validated.
            raw_bgen.unlink(missing_ok=True)
            raw_sample.unlink(missing_ok=True)
            for suffix in [".bed", ".bim", ".fam", ".log"]:
                Path(
                    str(temporary_bed_prefix) + suffix
                ).unlink(missing_ok=True)

            reconstruction_record = {
                **reconstruction_spec,
                "fingerprint": reconstruction_fingerprint,
                "raw_controlled_metadata": raw_metadata,
                "reconstructed_bgen": file_metadata(bgen_path),
                "reconstructed_sample": file_metadata(
                    sample_path,
                    include_sha256=True,
                ),
            }
            atomic_write_json(
                reconstruction_sentinel,
                reconstruction_record,
            )

        status = (
            "reconstructed_with_original_variant_qc_"
            "and_hard_call_bed_roundtrip"
        )
        source_fingerprint_payload = {
            "chromosome": chromosome,
            "source_status": status,
            "reconstruction_fingerprint": reconstruction_fingerprint,
            "bgen": file_metadata(
                bgen_path,
                include_sha256=True,
            ),
            "sample": file_metadata(
                sample_path,
                include_sha256=True,
            ),
        }

    qc_sources[chromosome] = {
        "bgen": str(bgen_path),
        "sample": str(sample_path),
        "status": status,
        "source_fingerprint": stable_fingerprint(
            source_fingerprint_payload
        ),
        "source_fingerprint_payload": source_fingerprint_payload,
    }
    inventory_rows.append(
        {
            "chromosome": chromosome,
            "source_status": status,
            "bgen_path": str(bgen_path),
            "bgen_bytes": bgen_path.stat().st_size,
            "sample_path": str(sample_path),
            "sample_bytes": sample_path.stat().st_size,
            "source_fingerprint": qc_sources[chromosome][
                "source_fingerprint"
            ],
        }
    )

inventory = pd.DataFrame(inventory_rows)
inventory.to_csv(
    AGG_DIR / "corrected_qc_genotype_inventory.tsv",
    sep="\t",
    index=False,
)
atomic_write_json(
    PRIVATE / "qc_sources.json",
    qc_sources,
)

if inventory.shape[0] != 23:
    raise RuntimeError("Did not resolve all 23 chr1-22/X QC pairs.")

(AGG_DIR / "missing_retained_qc_pairs.tsv").unlink(
    missing_ok=True
)
print("\nAll chromosome QC inputs are available.")
display(inventory)
''',
r'''
# ============================================================
# Block R3 — Build one merged autosomal LD-pruned Step 1 PGEN
# ============================================================

qc_sources = json.loads(
    (PRIVATE / "qc_sources.json").read_text()
)
setup_summary = json.loads(
    (AGG_DIR / "corrected_run_setup.json").read_text()
)

AUTOSOMES = [str(chrom) for chrom in range(1, 23)]
prune_dir = PRIVATE / "ld_pruning"
prune_dir.mkdir(parents=True, exist_ok=True)

part_rows = []
expected_ids = pd.read_csv(
    FROZEN_PHENO,
    sep="\t",
    dtype={"FID": str, "IID": str},
)[["FID", "IID"]]
expected_id_order_sha256 = stable_fingerprint(
    expected_ids.values.tolist()
)
if (
    expected_id_order_sha256
    != setup_summary["frozen_id_order_sha256"]
):
    raise RuntimeError(
        "Frozen phenotype ID order no longer matches Block R1."
    )


def read_pvar_rows(path):
    """Read the first five required PVAR fields, ignoring ## metadata."""
    return pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["CHROM", "POS", "ID", "REF", "ALT"],
        usecols=[0, 1, 2, 3, 4],
        dtype=str,
    )


def read_psam_ids(path):
    psam = pd.read_csv(
        path,
        sep=r"\s+",
        dtype=str,
    )
    fid_column = "#FID" if "#FID" in psam.columns else "FID"
    if fid_column not in psam.columns or "IID" not in psam.columns:
        raise RuntimeError(
            f"PSAM is missing FID/IID columns: {path}"
        )
    ids = (
        psam[[fid_column, "IID"]]
        .rename(columns={fid_column: "FID"})
        .reset_index(drop=True)
    )
    if ids.duplicated().any():
        raise RuntimeError(f"Duplicate PSAM IDs detected: {path}")
    return ids


def assert_exact_psam_ids(path, label):
    observed = read_psam_ids(path)
    if not observed.equals(expected_ids):
        raise RuntimeError(
            f"{label} PSAM does not match all 392,155 frozen "
            "FID/IID pairs in exact order."
        )
    return stable_fingerprint(observed.values.tolist())


def promote_prefix(partial_prefix, final_prefix, suffixes):
    for suffix in suffixes:
        partial_path = Path(str(partial_prefix) + suffix)
        final_path = Path(str(final_prefix) + suffix)
        if not nonempty(partial_path):
            raise RuntimeError(
                f"Missing partial output before promotion: {partial_path}"
            )
        os.replace(partial_path, final_path)


def run_checked_autosomal_disk_guard(
    command,
    description,
    log_path,
    cleanup_prefixes,
):
    command = [str(value) for value in command]
    log_path = Path(log_path)

    cleanup_prefixes = [
        Path(value)
        for value in cleanup_prefixes
    ]

    emergency_floor_bytes = int(
        float(
            PRODUCTION_CONFIG["runtime_floor_gib"]
        )
        * 1024**3
    )

    stop_trigger_bytes = (
        emergency_floor_bytes
        + 10 * 1024**3
    )

    initial_free_bytes = shutil.disk_usage(
        "/home/jupyter"
    ).free

    if initial_free_bytes <= stop_trigger_bytes:
        raise RuntimeError(
            f"{description} cannot start safely: "
            f"{initial_free_bytes / 1024**3:.1f} GiB free."
        )

    print(f"\\n[{description}]")
    print(" ".join(command))

    with log_path.open("w") as log_handle:
        process = subprocess.Popen(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

        last_state_update = 0.0

        while process.poll() is None:
            time.sleep(5)

            free_bytes = shutil.disk_usage(
                "/home/jupyter"
            ).free

            now = time.monotonic()

            if now - last_state_update >= 30:
                production_state_update(
                    status="running",
                    stage="build_step1_backbone",
                    current_stage_detail=(
                        "disk_guarded_"
                        "autosomal_preparation"
                    ),
                    autosomal_preparation_description=(
                        description
                    ),
                    autosomal_preparation_free_gib=round(
                        free_bytes / 1024**3,
                        2,
                    ),
                )

                last_state_update = now

            if free_bytes < stop_trigger_bytes:
                process.terminate()

                try:
                    process.wait(timeout=30)

                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=30)

                log_handle.flush()

                for prefix in cleanup_prefixes:
                    for candidate in prefix.parent.glob(
                        prefix.name + "*"
                    ):
                        if candidate.is_file():
                            candidate.unlink(
                                missing_ok=True
                            )

                raise RuntimeError(
                    f"{description} was stopped before "
                    "the disk reserve could be breached. "
                    "Partial files were removed. "
                    f"Free space at stop: "
                    f"{free_bytes / 1024**3:.1f} GiB."
                )

        return_code = process.returncode

    combined = log_path.read_text(
        errors="replace"
    )

    print(combined)

    if return_code != 0:
        raise RuntimeError(
            f"{description} failed with exit code "
            f"{return_code}."
        )

    return combined


STEP1_PREP_MAX_WORKERS = 1
STEP1_PREP_THREADS_PER_JOB = 24
STEP1_PREP_MEMORY_MIB_PER_JOB = 192000


def build_one_step1_part(chromosome):
    source = qc_sources[chromosome]
    prune_prefix = prune_dir / f"chr{chromosome}_ld"
    part_prefix = STEP1_PARTS / f"chr{chromosome}_step1"
    prune_in = Path(str(prune_prefix) + ".prune.in")
    part_pgen = Path(str(part_prefix) + ".pgen")
    part_pvar = Path(str(part_prefix) + ".pvar")
    part_psam = Path(str(part_prefix) + ".psam")
    part_sentinel = Path(str(part_prefix) + ".complete.json")
    prune_command_core = [
        "plink2",
        "--bgen", source["bgen"], "ref-first",
        "--sample", source["sample"],
        "--keep", str(KEEP_FILE),
        "--indep-pairwise", "200kb", "1", "0.5",
        "--threads", "24",
        "--memory", "192000",
    ]
    part_command_core = [
        "plink2",
        "--bgen", source["bgen"], "ref-first",
        "--sample", source["sample"],
        "--keep", str(KEEP_FILE),
        "--extract", "<CURRENT_PRUNE_IN>",
        "--mac", "100",
        "--indiv-sort", "file", str(KEEP_FILE),
        "--make-pgen",
        "--sort-vars",
        "--threads", "24",
        "--memory", "192000",
    ]
    part_spec = {
        "chromosome": chromosome,
        "qc_source_fingerprint": source["source_fingerprint"],
        "keep_sha256": setup_summary[
            "frozen_input_hashes"
        ]["keep_sha256"],
        "expected_id_order_sha256": expected_id_order_sha256,
        "plink2_version": PLINK2_VERSION,
        "prune_command": prune_command_core,
        "part_command": part_command_core,
    }
    part_fingerprint = stable_fingerprint(part_spec)
    part_ready = (
        all(
            nonempty(path)
            for path in [part_pgen, part_pvar, part_psam, prune_in]
        )
        and sentinel_matches(part_sentinel, part_fingerprint)
    )
    if part_ready:
        try:
            prior_part_record = json.loads(
                part_sentinel.read_text()
            )
            if not (
                recorded_file_matches(
                    part_pgen,
                    prior_part_record.get("pgen"),
                )
                and recorded_file_matches(
                    part_pvar,
                    prior_part_record.get("pvar"),
                    include_sha256=True,
                )
                and recorded_file_matches(
                    part_psam,
                    prior_part_record.get("psam"),
                    include_sha256=True,
                )
            ):
                raise RuntimeError(
                    "Output metadata differs from the completion record."
                )
            assert_exact_psam_ids(
                part_psam,
                f"Existing chr{chromosome} Step 1 part",
            )
            if read_pvar_rows(part_pvar).empty:
                raise RuntimeError("PVAR contains no variants.")
        except Exception as exc:
            print(
                f"Existing chr{chromosome} part is not reusable:",
                exc,
            )
            part_ready = False

    if not part_ready:
        partial_prune_prefix = (
            prune_dir / f"chr{chromosome}_ld.partial"
        )
        partial_part_prefix = (
            STEP1_PARTS / f"chr{chromosome}_step1.partial"
        )
        for prefix, suffixes in [
            (partial_prune_prefix, [".prune.in", ".prune.out", ".log"]),
            (partial_part_prefix, [".pgen", ".pvar", ".psam", ".log"]),
        ]:
            for suffix in suffixes:
                Path(str(prefix) + suffix).unlink(missing_ok=True)

        run_checked_autosomal_disk_guard(
            prune_command_core
            + ["--out", str(partial_prune_prefix)],
            f"LD-prune chr{chromosome}",
            LOG_DIR / f"ld_prune_chr{chromosome}.log",
            cleanup_prefixes=[partial_prune_prefix],
        )
        partial_prune_in = Path(
            str(partial_prune_prefix) + ".prune.in"
        )
        if not nonempty(partial_prune_in):
            raise RuntimeError(
                f"LD pruning produced no retained variants on "
                f"chr{chromosome}."
            )

        executable_part_command = [
            value
            if value != "<CURRENT_PRUNE_IN>"
            else str(partial_prune_in)
            for value in part_command_core
        ]
        run_checked_autosomal_disk_guard(
            executable_part_command
            + ["--out", str(partial_part_prefix)],
            f"Create chr{chromosome} Step 1 PGEN part",
            LOG_DIR / f"make_step1_part_chr{chromosome}.log",
            cleanup_prefixes=[partial_part_prefix],
        )
        partial_part_pvar = Path(
            str(partial_part_prefix) + ".pvar"
        )
        partial_part_psam = Path(
            str(partial_part_prefix) + ".psam"
        )
        if read_pvar_rows(partial_part_pvar).empty:
            raise RuntimeError(
                f"No Step 1 variants remained on chr{chromosome}."
            )
        psam_hash = assert_exact_psam_ids(
            partial_part_psam,
            f"New chr{chromosome} Step 1 part",
        )
        promote_prefix(
            partial_prune_prefix,
            prune_prefix,
            [".prune.in", ".prune.out"],
        )
        promote_prefix(
            partial_part_prefix,
            part_prefix,
            [".pgen", ".pvar", ".psam"],
        )
        atomic_write_json(
            part_sentinel,
            {
                **part_spec,
                "fingerprint": part_fingerprint,
                "psam_id_order_sha256": psam_hash,
                "pgen": file_metadata(part_pgen),
                "pvar": file_metadata(
                    part_pvar,
                    include_sha256=True,
                ),
                "psam": file_metadata(
                    part_psam,
                    include_sha256=True,
                ),
            },
        )

    pvar = read_pvar_rows(part_pvar)
    n_variants = int(pvar.shape[0])
    if n_variants == 0:
        raise RuntimeError(
            f"No Step 1 variants remained on chr{chromosome}."
        )

    part_rows.append(
        {
            "chromosome": chromosome,
            "n_step1_variants": n_variants,
            "part_prefix": str(part_prefix),
            "part_fingerprint": part_fingerprint,
            "psam_id_order_sha256": assert_exact_psam_ids(
                part_psam,
                f"Final chr{chromosome} Step 1 part",
            ),
        }
    )

    return chromosome


with concurrent.futures.ThreadPoolExecutor(
    max_workers=STEP1_PREP_MAX_WORKERS
) as executor:
    futures = {
        executor.submit(
            build_one_step1_part,
            chromosome,
        ): chromosome
        for chromosome in AUTOSOMES
    }

    try:
        for future in concurrent.futures.as_completed(
            futures
        ):
            chromosome = futures[future]
            future.result()
            print(
                f"[chr{chromosome}] Step 1 preparation "
                "validated."
            )
    except BaseException:
        for future in futures:
            future.cancel()
        executor.shutdown(
            wait=True,
            cancel_futures=True,
        )
        raise


part_summary = pd.DataFrame(part_rows)
part_summary["sort_order"] = pd.to_numeric(
    part_summary["chromosome"]
)
part_summary = (
    part_summary
    .sort_values("sort_order")
    .drop(columns="sort_order")
    .reset_index(drop=True)
)
part_summary.to_csv(
    AGG_DIR / "corrected_step1_marker_counts_by_chromosome.tsv",
    sep="\t",
    index=False,
)

merge_list = STEP1_DIR / "autosomal_step1_pmerge_list.txt"
merge_list.write_text(
    "\n".join(part_summary["part_prefix"].tolist()) + "\n"
)

STEP1_PGEN = STEP1_DIR / "psoriasis_autosomal_step1"
merge_command_core = [
    "plink2",
    "--pmerge-list", str(merge_list), "pfile",
    "--sample-inner-join",
    "--make-pgen",
    "--sort-vars",
    "--threads", "24",
]
merge_spec = {
    "part_fingerprints": part_summary[
        "part_fingerprint"
    ].tolist(),
    "merge_command": merge_command_core,
    "expected_id_order_sha256": expected_id_order_sha256,
    "plink2_version": PLINK2_VERSION,
}
merge_fingerprint = stable_fingerprint(merge_spec)

# Mandatory chr1+chr2 syntax and compatibility dry run.
dry_merge_list = STEP1_DIR / "dryrun_chr1_chr2_pmerge_list.txt"
dry_merge_list.write_text(
    "\n".join(part_summary["part_prefix"].iloc[:2].tolist()) + "\n"
)
dry_prefix = STEP1_DIR / "dryrun_chr1_chr2"
dry_sentinel = Path(str(dry_prefix) + ".complete.json")
dry_command_core = [
    "plink2",
    "--pmerge-list", str(dry_merge_list), "pfile",
    "--sample-inner-join",
    "--make-pgen",
    "--sort-vars",
    "--threads", "24",
]
dry_spec = {
    "part_fingerprints": part_summary[
        "part_fingerprint"
    ].iloc[:2].tolist(),
    "merge_command": dry_command_core,
    "expected_id_order_sha256": expected_id_order_sha256,
    "plink2_version": PLINK2_VERSION,
}
dry_fingerprint = stable_fingerprint(dry_spec)
dry_ready = (
    all(
        nonempty(Path(str(dry_prefix) + suffix))
        for suffix in [".pgen", ".pvar", ".psam"]
    )
    and sentinel_matches(dry_sentinel, dry_fingerprint)
)
if dry_ready:
    try:
        prior_dry_record = json.loads(dry_sentinel.read_text())
        if not (
            recorded_file_matches(
                Path(str(dry_prefix) + ".pgen"),
                prior_dry_record.get("pgen"),
            )
            and recorded_file_matches(
                Path(str(dry_prefix) + ".pvar"),
                prior_dry_record.get("pvar"),
                include_sha256=True,
            )
            and recorded_file_matches(
                Path(str(dry_prefix) + ".psam"),
                prior_dry_record.get("psam"),
                include_sha256=True,
            )
        ):
            raise RuntimeError(
                "Merge dry-run output metadata differs from its "
                "completion record."
            )
        assert_exact_psam_ids(
            Path(str(dry_prefix) + ".psam"),
            "Existing chr1+chr2 merge dry run",
        )
    except Exception as exc:
        print("Existing merge dry run is not reusable:", exc)
        dry_ready = False

if not dry_ready:
    dry_partial = STEP1_DIR / "dryrun_chr1_chr2.partial"
    for suffix in [".pgen", ".pvar", ".psam", ".log"]:
        Path(str(dry_partial) + suffix).unlink(missing_ok=True)
    run_checked(
        dry_command_core + ["--out", str(dry_partial)],
        "Mandatory chr1+chr2 PLINK2 merge dry run",
        LOG_DIR / "dryrun_merge_chr1_chr2.log",
    )
    dry_chromosomes = {
        str(value).replace("chr", "")
        for value in read_pvar_rows(
            Path(str(dry_partial) + ".pvar")
        )["CHROM"].unique()
    }
    if dry_chromosomes != {"1", "2"}:
        raise RuntimeError(
            "The mandatory merge dry run did not contain exactly "
            "chromosomes 1 and 2."
        )
    assert_exact_psam_ids(
        Path(str(dry_partial) + ".psam"),
        "Mandatory chr1+chr2 merge dry run",
    )
    promote_prefix(
        dry_partial,
        dry_prefix,
        [".pgen", ".pvar", ".psam"],
    )
    atomic_write_json(
        dry_sentinel,
        {
            **dry_spec,
            "fingerprint": dry_fingerprint,
            "pgen": file_metadata(
                Path(str(dry_prefix) + ".pgen")
            ),
            "pvar": file_metadata(
                Path(str(dry_prefix) + ".pvar"),
                include_sha256=True,
            ),
            "psam": file_metadata(
                Path(str(dry_prefix) + ".psam"),
                include_sha256=True,
            ),
        },
    )

merged_files = [
    Path(str(STEP1_PGEN) + suffix)
    for suffix in [".pgen", ".pvar", ".psam"]
]
merged_sentinel = Path(str(STEP1_PGEN) + ".complete.json")
merged_ready = (
    all(nonempty(path) for path in merged_files)
    and sentinel_matches(merged_sentinel, merge_fingerprint)
)
if merged_ready:
    try:
        prior_merged_record = json.loads(
            merged_sentinel.read_text()
        )["validated_step1_input"]
        if not (
            recorded_file_matches(
                Path(str(STEP1_PGEN) + ".pgen"),
                prior_merged_record.get("pgen"),
            )
            and recorded_file_matches(
                Path(str(STEP1_PGEN) + ".pvar"),
                prior_merged_record.get("pvar"),
                include_sha256=True,
            )
            and recorded_file_matches(
                Path(str(STEP1_PGEN) + ".psam"),
                prior_merged_record.get("psam"),
                include_sha256=True,
            )
        ):
            raise RuntimeError(
                "Merged output metadata differs from its "
                "completion record."
            )
        assert_exact_psam_ids(
            Path(str(STEP1_PGEN) + ".psam"),
            "Existing merged autosomal Step 1 input",
        )
    except Exception as exc:
        print("Existing merged Step 1 input is not reusable:", exc)
        merged_ready = False

if not merged_ready:
    partial_step1 = STEP1_DIR / "psoriasis_autosomal_step1.partial"
    for suffix in [".pgen", ".pvar", ".psam", ".log"]:
        Path(str(partial_step1) + suffix).unlink(missing_ok=True)
    run_checked(
        merge_command_core + ["--out", str(partial_step1)],
        "Merge chr1-22 into one autosomal Step 1 PGEN",
        LOG_DIR / "merge_autosomal_step1_pgen.log",
    )
    assert_exact_psam_ids(
        Path(str(partial_step1) + ".psam"),
        "New merged autosomal Step 1 input",
    )
    promote_prefix(
        partial_step1,
        STEP1_PGEN,
        [".pgen", ".pvar", ".psam"],
    )

merged_pvar = read_pvar_rows(
    Path(str(STEP1_PGEN) + ".pvar")
)
merged_psam = pd.read_csv(
    Path(str(STEP1_PGEN) + ".psam"),
    sep="\t",
    dtype=str,
)

observed_chromosomes = {
    str(value).replace("chr", "")
    for value in merged_pvar["CHROM"].unique()
}

if observed_chromosomes != set(AUTOSOMES):
    raise RuntimeError(
        "Merged Step 1 PGEN does not contain exactly chr1-22: "
        f"{sorted(observed_chromosomes)}"
    )
merged_id_order_sha256 = assert_exact_psam_ids(
    Path(str(STEP1_PGEN) + ".psam"),
    "Final merged autosomal Step 1 input",
)

n_step1_variants = int(merged_pvar.shape[0])
if not 50000 <= n_step1_variants <= 1000000:
    raise RuntimeError(
        "STOP FOR REVIEW: merged Step 1 marker count is outside "
        "the prespecified 50,000-1,000,000 range: "
        f"{n_step1_variants:,}."
    )

step1_input_validation = {
    "step1_input_prefix": str(STEP1_PGEN),
    "n_samples": int(merged_psam.shape[0]),
    "n_variants": n_step1_variants,
    "chromosomes": sorted(
        observed_chromosomes,
        key=int,
    ),
    "ld_pruning": "200kb window, 1-variant step, r2 < 0.5",
    "minimum_mac": 100,
    "uses_autosomes_only": True,
    "expected_id_order_sha256": expected_id_order_sha256,
    "observed_id_order_sha256": merged_id_order_sha256,
    "merge_fingerprint": merge_fingerprint,
    "merge_dry_run_passed": True,
    "pgen": file_metadata(
        Path(str(STEP1_PGEN) + ".pgen")
    ),
    "pvar": file_metadata(
        Path(str(STEP1_PGEN) + ".pvar"),
        include_sha256=True,
    ),
    "psam": file_metadata(
        Path(str(STEP1_PGEN) + ".psam"),
        include_sha256=True,
    ),
}
step1_input_validation["fingerprint"] = stable_fingerprint(
    step1_input_validation
)
atomic_write_json(
    AGG_DIR / "corrected_step1_input_validation.json",
    step1_input_validation,
)
atomic_write_json(
    merged_sentinel,
    {
        **merge_spec,
        "fingerprint": merge_fingerprint,
        "validated_step1_input": step1_input_validation,
    },
)

print("\nAutosomal Step 1 PGEN validation passed.")
print(json.dumps(step1_input_validation, indent=2))
display(part_summary)
''',
r'''
# ============================================================
# Block R4 — Run and validate the single whole-genome Step 1
# ============================================================

STEP1_PGEN = STEP1_DIR / "psoriasis_autosomal_step1"
STEP1_PREFIX = STEP1_DIR / "psoriasis_primary_loco_corrected"
PRED_LIST = Path(str(STEP1_PREFIX) + "_pred.list")
LOWMEM_PREFIX = STEP1_DIR / "psoriasis_primary_lowmem"
STEP1_STDOUT_LOG = (
    LOG_DIR / "corrected_primary_regenie_step1.log"
)
STEP1_SENTINEL = Path(str(STEP1_PREFIX) + ".complete.json")
step1_input_validation = json.loads(
    (AGG_DIR / "corrected_step1_input_validation.json").read_text()
)
setup_summary = json.loads(
    (AGG_DIR / "corrected_run_setup.json").read_text()
)

step1_command_core = [
    "regenie",
    "--step", "1",
    "--bt",
    "--pgen", str(STEP1_PGEN),
    "--phenoFile", str(FROZEN_PHENO),
    "--phenoCol", "psoriasis_case",
    "--covarFile", str(FROZEN_COVAR),
    "--bsize", "2000",
    "--threads", "48",
    "--out", "<OUTPUT_PREFIX>",
]
step1_spec = {
    "step1_input_fingerprint": step1_input_validation["fingerprint"],
    "frozen_input_hashes": setup_summary["frozen_input_hashes"],
    "regenie_version": REGENIE_VERSION,
    "step1_command": step1_command_core,
}
step1_fingerprint = stable_fingerprint(step1_spec)

merged_step1_ids = read_psam_ids(
    Path(str(STEP1_PGEN) + ".psam")
)
if not merged_step1_ids.equals(expected_ids):
    raise RuntimeError(
        "Merged Step 1 PSAM no longer matches the frozen cohort."
    )
expected_loco_header_ids = (
    merged_step1_ids["FID"]
    + "_"
    + merged_step1_ids["IID"]
).tolist()
if len(expected_loco_header_ids) != len(
    set(expected_loco_header_ids)
):
    raise RuntimeError(
        "Expected LOCO FID_IID header identifiers are not unique."
    )


def validate_step1_log(log_path):
    log_path = Path(log_path)
    if not nonempty(log_path):
        raise RuntimeError("Step 1 stdout log is missing.")
    text = log_path.read_text(errors="replace")
    fatal = re.search(
        r"(?im)\b(?:ERROR:|FATAL|segmentation fault|"
        r"terminate called|killed|std::bad_alloc)\b",
        text,
    )
    if fatal:
        raise RuntimeError(
            f"Fatal Step 1 log text detected: {fatal.group(0)}"
        )
    required_patterns = [
        (
            r"(?m)^\s*-number of phenotyped individuals with no "
            r"missing data\s*=\s*392155\s*$"
        ),
        (
            r"(?m)^\s*-number of individuals with covariate data"
            r"\s*=\s*392155\s*$"
        ),
        (
            r"(?m)^\s*-number of individuals used in analysis"
            r"\s*=\s*392155\s*$"
        ),
        (
            r"(?m)^\s*-\s*'psoriasis_case':\s*5805 cases and "
            r"386350 controls\s*$"
        ),
        (
            r"(?m)^\s*List of blup files written to:\s*"
            r"\[([^\]]+_pred\.list)\]\s*$"
        ),
        r"(?m)^\s*Elapsed time\s*:\s*[0-9.]+s\s*$",
        r"(?m)^\s*End time:\s*.+$",
    ]
    missing_patterns = [
        pattern
        for pattern in required_patterns
        if re.search(pattern, text) is None
    ]
    if missing_patterns:
        raise RuntimeError(
            "Step 1 log is missing required completion evidence:\n"
            + "\n".join(missing_patterns)
        )
    return {
        "stdout_log": str(log_path),
        "stdout_log_sha256": sha256_file(log_path),
        "n_participants": 392155,
        "n_cases": 5805,
        "n_controls": 386350,
        "successful_end_marker": True,
    }


def validate_loco_prediction_files(
    pred_list_path,
    step1_log_path=None,
):
    pred_list_path = Path(pred_list_path)
    if not nonempty(pred_list_path):
        raise RuntimeError("REGENIE pred.list is missing or empty.")

    entries = []
    for line in pred_list_path.read_text().splitlines():
        fields = line.strip().split()
        if not fields:
            continue
        if len(fields) != 2:
            raise RuntimeError(
                "Each pred.list row must contain phenotype and path: "
                f"{line!r}"
            )
        loco_path = Path(fields[1])
        if not loco_path.is_absolute():
            relative_candidate = (
                pred_list_path.parent / loco_path
            )
            if relative_candidate.exists():
                loco_path = relative_candidate
        entries.append((fields[0], loco_path))

    if len(entries) != 1:
        raise RuntimeError(
            f"Expected one phenotype in pred.list; found {len(entries)}."
        )
    phenotype_name, loco_path = entries[0]
    if phenotype_name != "psoriasis_case":
        raise RuntimeError(
            f"Unexpected pred.list phenotype: {phenotype_name}"
        )
    if not nonempty(loco_path):
        raise RuntimeError(
            f"LOCO prediction file is missing: {loco_path}"
        )

    row_labels = []
    n_lines = 0
    header_ids = None
    with loco_path.open() as handle:
        for line_number, line in enumerate(handle):
            n_lines += 1
            fields = line.split()
            if not fields:
                raise RuntimeError(
                    f"Blank row detected at LOCO line {line_number + 1}."
                )
            prediction_count = len(fields) - 1
            if prediction_count != 392155:
                row_name = fields[0]
                raise RuntimeError(
                    f"LOCO row {row_name!r} contains "
                    f"{prediction_count:,} predictions; "
                    "expected 392,155."
                )
            if line_number == 0:
                header_ids = fields[1:]
                if len(header_ids) != len(set(header_ids)):
                    raise RuntimeError(
                        "LOCO header identifiers are not unique."
                    )
                if header_ids != expected_loco_header_ids:
                    raise RuntimeError(
                        "LOCO header FID_IID sequence does not match "
                        "the merged Step 1 PSAM exactly."
                    )
            else:
                try:
                    predictions = np.asarray(
                        fields[1:],
                        dtype=float,
                    )
                except Exception as exc:
                    raise RuntimeError(
                        f"LOCO row {fields[0]!r} contains a "
                        "non-numeric prediction."
                    ) from exc
                if not np.isfinite(predictions).all():
                    raise RuntimeError(
                        f"LOCO row {fields[0]!r} contains a "
                        "non-finite prediction."
                    )
                row_labels.append(fields[0])

    expected_labels = [str(value) for value in range(1, 24)]
    if n_lines != 24:
        raise RuntimeError(
            "Expected one header plus 23 LOCO rows; "
            f"observed {n_lines} lines."
        )
    if row_labels != expected_labels:
        raise RuntimeError(
            "LOCO chromosome rows are not exactly 1-23: "
            f"{row_labels}"
        )

    validation = {
        "phenotype": phenotype_name,
        "loco_path": str(loco_path),
        "loco_sha256": sha256_file(loco_path),
        "loco_header_id_order_sha256": stable_fingerprint(
            header_ids
        ),
        "n_loco_rows": len(row_labels),
        "loco_chromosome_rows": row_labels,
    }
    if step1_log_path is not None:
        validation["step1_log_validation"] = validate_step1_log(
            step1_log_path
        )
    return validation


step1_ready = (
    nonempty(PRED_LIST)
    and nonempty(STEP1_STDOUT_LOG)
    and sentinel_matches(
        STEP1_SENTINEL,
        step1_fingerprint,
    )
)
if step1_ready:
    try:
        prior_step1_record = json.loads(
            STEP1_SENTINEL.read_text()
        )
        if not (
            recorded_file_matches(
                PRED_LIST,
                prior_step1_record.get("prediction_list"),
                include_sha256=True,
            )
            and recorded_file_matches(
                STEP1_STDOUT_LOG,
                prior_step1_record.get("stdout_log"),
                include_sha256=True,
            )
        ):
            raise RuntimeError(
                "Step 1 output metadata differs from the "
                "completion record."
            )
        loco_validation = validate_loco_prediction_files(
            PRED_LIST,
            STEP1_STDOUT_LOG,
        )
        current_loco_path = Path(
            loco_validation["loco_path"]
        ).resolve()
        recorded_loco = prior_step1_record.get("loco")
        if (
            not recorded_loco
            or current_loco_path
            != Path(recorded_loco["path"]).resolve()
            or not recorded_file_matches(
                current_loco_path,
                recorded_loco,
                include_sha256=True,
            )
        ):
            raise RuntimeError(
                "The LOCO file resolved from pred.list differs "
                "from the completed Step 1 record."
            )
        print("Validated existing corrected Step 1 outputs.")
    except Exception as exc:
        print("Existing Step 1 output is not reusable:", exc)
        step1_ready = False

if not step1_ready:
    partial_step1_prefix = (
        STEP1_DIR / "psoriasis_primary_loco_corrected.partial"
    )
    partial_lowmem_prefix = (
        STEP1_DIR / "psoriasis_primary_lowmem.partial"
    )
    partial_pred_list = Path(
        str(partial_step1_prefix) + "_pred.list"
    )
    partial_stdout_log = Path(
        str(STEP1_STDOUT_LOG) + ".partial"
    )
    for path in [
        partial_pred_list,
        Path(str(partial_step1_prefix) + "_1.loco"),
        Path(str(partial_step1_prefix) + "_1.loco.gz"),
        Path(str(partial_step1_prefix) + ".log"),
        partial_stdout_log,
    ]:
        path.unlink(missing_ok=True)
    for path in STEP1_DIR.glob(
        partial_lowmem_prefix.name + "*"
    ):
        if path.is_file():
            path.unlink()

    executable_step1_command = [
        (
            str(partial_lowmem_prefix)
            if value == "<LOWMEM_PREFIX>"
            else str(partial_step1_prefix)
            if value == "<OUTPUT_PREFIX>"
            else value
        )
        for value in step1_command_core
    ]
    run_checked(
        executable_step1_command,
        "Run corrected autosomal whole-genome REGENIE Step 1",
        partial_stdout_log,
    )
    partial_validation = validate_loco_prediction_files(
        partial_pred_list,
        partial_stdout_log,
    )
    partial_loco_path = Path(partial_validation["loco_path"])
    final_loco_path = Path(
        str(partial_loco_path).replace(
            str(partial_step1_prefix),
            str(STEP1_PREFIX),
            1,
        )
    )
    os.replace(partial_loco_path, final_loco_path)
    PRED_LIST.write_text(
        f"psoriasis_case {final_loco_path}\n"
    )
    partial_regenie_log = Path(
        str(partial_step1_prefix) + ".log"
    )
    final_regenie_log = Path(str(STEP1_PREFIX) + ".log")
    if nonempty(partial_regenie_log):
        os.replace(partial_regenie_log, final_regenie_log)
    os.replace(partial_stdout_log, STEP1_STDOUT_LOG)
    loco_validation = validate_loco_prediction_files(
        PRED_LIST,
        STEP1_STDOUT_LOG,
    )
    atomic_write_json(
        STEP1_SENTINEL,
        {
            **step1_spec,
            "fingerprint": step1_fingerprint,
            "prediction_list": file_metadata(
                PRED_LIST,
                include_sha256=True,
            ),
            "stdout_log": file_metadata(
                STEP1_STDOUT_LOG,
                include_sha256=True,
            ),
            "loco": file_metadata(
                Path(loco_validation["loco_path"]),
                include_sha256=True,
            ),
        },
    )

step1_model_validation = {
    **step1_input_validation,
    **loco_validation,
    "n_cases": 5805,
    "n_controls": 386350,
    "covariates": covar_analysis.columns[2:].tolist(),
    "regenie_version": REGENIE_VERSION,
    "prediction_list_reused_for_all_step2_chromosomes": True,
    "step1_fingerprint": step1_fingerprint,
    "prediction_list_sha256": sha256_file(PRED_LIST),
    "loco_header_matches_merged_psam_exactly": True,
    "all_loco_predictions_finite": True,
}
atomic_write_json(
    AGG_DIR / "corrected_step1_model_validation.json",
    step1_model_validation,
)

print("\nCorrected Step 1 validation passed.")
print(json.dumps(step1_model_validation, indent=2))
''',
r'''
# ============================================================
# Block R5 — Run corrected Step 2 on chr1-22 and X
# ============================================================

qc_sources = json.loads(
    (PRIVATE / "qc_sources.json").read_text()
)
PRED_LIST = STEP1_DIR / "psoriasis_primary_loco_corrected_pred.list"
setup_summary = json.loads(
    (AGG_DIR / "corrected_run_setup.json").read_text()
)
step1_model_validation = json.loads(
    (AGG_DIR / "corrected_step1_model_validation.json").read_text()
)

# Prevent any frozen-input drift between Step 1 and Step 2. This also
# guarantees that the chromosome X sex mapping is derived from the same
# covariate bytes fingerprinted before the LOCO model was fit.
current_frozen_input_hashes = {
    "phenotype_sha256": sha256_file(FROZEN_PHENO),
    "covariate_sha256": sha256_file(FROZEN_COVAR),
    "keep_sha256": sha256_file(KEEP_FILE),
}
if (
    current_frozen_input_hashes
    != setup_summary["frozen_input_hashes"]
):
    raise RuntimeError(
        "A frozen phenotype, covariate, or keep file changed after "
        "Block R1. Do not start Step 2."
    )

# Revalidate before any association test.
loco_validation = validate_loco_prediction_files(
    PRED_LIST,
    STEP1_STDOUT_LOG,
)
if loco_validation["n_loco_rows"] != 23:
    raise RuntimeError("Corrected 23-row LOCO validation failed.")
step1_completion_record = json.loads(
    STEP1_SENTINEL.read_text()
)
if not sentinel_matches(
    STEP1_SENTINEL,
    step1_model_validation["step1_fingerprint"],
):
    raise RuntimeError(
        "The Step 1 completion sentinel no longer matches the "
        "validated Step 1 fingerprint."
    )
resolved_loco_path = Path(
    loco_validation["loco_path"]
).resolve()
recorded_loco = step1_completion_record.get("loco")
if not (
    recorded_file_matches(
        PRED_LIST,
        step1_completion_record.get("prediction_list"),
        include_sha256=True,
    )
    and recorded_file_matches(
        STEP1_STDOUT_LOG,
        step1_completion_record.get("stdout_log"),
        include_sha256=True,
    )
    and recorded_loco
    and resolved_loco_path
    == Path(recorded_loco["path"]).resolve()
    and recorded_file_matches(
        resolved_loco_path,
        recorded_loco,
        include_sha256=True,
    )
    and sha256_file(PRED_LIST)
    == step1_model_validation["prediction_list_sha256"]
    and loco_validation["loco_sha256"]
    == step1_model_validation["loco_sha256"]
    and loco_validation["loco_header_id_order_sha256"]
    == step1_model_validation[
        "loco_header_id_order_sha256"
    ]
):
    raise RuntimeError(
        "Step 1 prediction artifacts differ from the validated "
        "completion record. Do not start Step 2."
    )

# REGENIE uses the Oxford SAMPLE sex field—not the covariate file—to
# distinguish males during chromosome X testing. Absence of a warning is not
# sufficient evidence that sex was read correctly: an absent or all-unknown
# sex field can also suppress the warning. Build and positively validate a
# chromosome X SAMPLE file whose 392,155 analyzed participants have exact
# PLINK sex codes derived from the frozen covariate (0=female -> 2;
# 1=male -> 1). Non-analysis rows remain present to preserve BGEN row order,
# but their sex is set to unknown because they are excluded by phenotype.
EXPECTED_CHRX_PLINK_SEX_COUNTS = {1: 154238, 2: 237917}
CHR_X_SOURCE_SAMPLE = Path(qc_sources["X"]["sample"])
CHR_X_STEP2_SAMPLE = (
    PRIVATE / "chrX_analysis_sex_validated.sample"
)
CHR_X_SAMPLE_SEX_QC_PATH = (
    AGG_DIR / "corrected_chrX_sample_sex_validation.json"
)


def read_oxford_sample(path):
    path = Path(path)
    if not nonempty(path):
        raise RuntimeError(f"Oxford SAMPLE file is missing: {path}")

    with path.open() as handle:
        header_line = handle.readline().strip()
        type_line = handle.readline().strip()
    header_fields = header_line.split()
    type_fields = type_line.split()

    if (
        len(header_fields) < 2
        or header_fields[:2] != ["ID_1", "ID_2"]
        or len(type_fields) != len(header_fields)
        or type_fields[:2] != ["0", "0"]
    ):
        raise RuntimeError(
            "Oxford SAMPLE header/type structure is invalid for "
            f"{path}. Expected ID_1 ID_2 as the first fields and "
            "one type token per field."
        )

    sample = pd.read_csv(
        path,
        sep=r"\s+",
        skiprows=2,
        names=header_fields,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    if sample.empty or sample.shape[1] != len(header_fields):
        raise RuntimeError(
            f"Oxford SAMPLE data rows could not be parsed: {path}"
        )
    if sample[["ID_1", "ID_2"]].duplicated().any():
        raise RuntimeError(
            f"Duplicate Oxford SAMPLE FID/IID pairs detected: {path}"
        )
    return header_fields, type_fields, sample


def build_validated_chrX_sample():
    source_payload = qc_sources["X"][
        "source_fingerprint_payload"
    ]["sample"]
    if not recorded_file_matches(
        CHR_X_SOURCE_SAMPLE,
        source_payload,
        include_sha256=True,
    ):
        raise RuntimeError(
            "The chromosome X source SAMPLE file differs from the "
            "Block R2 source record."
        )

    source_header_fields, source_type_fields, source_sample = (
        read_oxford_sample(
            CHR_X_SOURCE_SAMPLE
        )
    )
    source_had_sex_field = "sex" in source_header_fields

    # REGENIE reads the fourth Oxford SAMPLE token as sex; it does not
    # locate the field by header name. Emit an exact four-field structure
    # regardless of any additional fields in the source SAMPLE.
    header_fields = ["ID_1", "ID_2", "missing", "sex"]
    type_fields = ["0", "0", "0", "D"]
    sample = source_sample[["ID_1", "ID_2"]].copy()
    sample["missing"] = "0"
    sample["sex"] = "0"

    if (
        source_header_fields[:2] != ["ID_1", "ID_2"]
        or source_type_fields[:2] != ["0", "0"]
    ):
        raise RuntimeError(
            "Chromosome X source SAMPLE does not begin with the "
            "required ID_1 ID_2 / 0 0 structure."
        )

    frozen_sex = pd.read_csv(
        FROZEN_COVAR,
        sep="\t",
        dtype={"FID": str, "IID": str},
    )[["FID", "IID", "sex"]]
    frozen_sex_numeric = pd.to_numeric(
        frozen_sex["sex"],
        errors="raise",
    )
    if (
        frozen_sex.shape[0] != 392155
        or frozen_sex[["FID", "IID"]].duplicated().any()
        or not np.isfinite(
            frozen_sex_numeric.to_numpy(dtype=float)
        ).all()
        or not frozen_sex_numeric.isin([0, 1]).all()
    ):
        raise RuntimeError(
            "Frozen chromosome X sex mapping is not a unique "
            "392,155-row exact 0/1 mapping."
        )
    frozen_sex["sex"] = frozen_sex_numeric.astype(int)

    expected_plink_sex = np.where(
        frozen_sex["sex"].to_numpy() == 1,
        1,
        2,
    )
    observed_counts = {
        int(code): int(count)
        for code, count in pd.Series(
            expected_plink_sex
        ).value_counts().sort_index().items()
    }
    if observed_counts != EXPECTED_CHRX_PLINK_SEX_COUNTS:
        raise RuntimeError(
            "Frozen chromosome X sex counts differ from the "
            "independently recovered 154,238 male / 237,917 female "
            f"analysis counts: {observed_counts}"
        )

    source_index = pd.MultiIndex.from_frame(
        sample[["ID_1", "ID_2"]]
    )
    frozen_index = pd.MultiIndex.from_frame(
        frozen_sex[["FID", "IID"]]
    )
    source_positions = source_index.get_indexer(frozen_index)
    if (source_positions < 0).any():
        missing_count = int((source_positions < 0).sum())
        raise RuntimeError(
            "Chromosome X SAMPLE file is missing "
            f"{missing_count:,} frozen analysis IDs."
        )

    sex_column_index = sample.columns.get_loc("sex")
    sample.iloc[
        source_positions,
        sex_column_index,
    ] = expected_plink_sex.astype(str)

    aligned_ids = (
        sample.iloc[source_positions][["ID_1", "ID_2"]]
        .rename(columns={"ID_1": "FID", "ID_2": "IID"})
        .reset_index(drop=True)
    )
    if not aligned_ids.equals(
        frozen_sex[["FID", "IID"]].reset_index(drop=True)
    ):
        raise RuntimeError(
            "Chromosome X SAMPLE IDs do not align exactly to the "
            "392,155 frozen FID/IID pairs."
        )
    aligned_sex = pd.to_numeric(
        sample.iloc[source_positions]["sex"],
        errors="raise",
    ).to_numpy(dtype=int)
    if not (
        np.isin(aligned_sex, [1, 2]).all()
        and np.array_equal(aligned_sex, expected_plink_sex)
    ):
        raise RuntimeError(
            "Chromosome X SAMPLE sex codes do not exactly match the "
            "frozen 0=female/1=male covariate mapping."
        )

    partial = Path(str(CHR_X_STEP2_SAMPLE) + ".partial")
    partial.unlink(missing_ok=True)
    partial.write_text(
        " ".join(header_fields)
        + "\n"
        + " ".join(type_fields)
        + "\n"
    )
    sample.to_csv(
        partial,
        sep=" ",
        index=False,
        header=False,
        mode="a",
        lineterminator="\n",
    )
    os.replace(partial, CHR_X_STEP2_SAMPLE)

    # Re-open the promoted file and repeat the cohort-level structural and
    # sex-code checks before it can be passed to REGENIE.
    promoted_header, promoted_types, promoted = read_oxford_sample(
        CHR_X_STEP2_SAMPLE
    )
    if (
        promoted_header != ["ID_1", "ID_2", "missing", "sex"]
        or promoted_types != ["0", "0", "0", "D"]
        or promoted.shape[0] != sample.shape[0]
    ):
        raise RuntimeError(
            "Promoted chromosome X SAMPLE structure differs from the "
            "validated in-memory representation."
        )
    promoted_index = pd.MultiIndex.from_frame(
        promoted[["ID_1", "ID_2"]]
    )
    promoted_positions = promoted_index.get_indexer(frozen_index)
    if (promoted_positions < 0).any():
        raise RuntimeError(
            "Promoted chromosome X SAMPLE lost frozen analysis IDs."
        )
    promoted_sex = pd.to_numeric(
        promoted.iloc[promoted_positions]["sex"],
        errors="raise",
    ).to_numpy(dtype=int)
    if not (
        np.isin(promoted_sex, [1, 2]).all()
        and np.array_equal(promoted_sex, expected_plink_sex)
    ):
        raise RuntimeError(
            "Promoted chromosome X SAMPLE sex validation failed."
        )

    validation = {
        "source_sample_path": str(CHR_X_SOURCE_SAMPLE),
        "source_sample_sha256": sha256_file(
            CHR_X_SOURCE_SAMPLE
        ),
        "frozen_covariate_sha256": (
            current_frozen_input_hashes["covariate_sha256"]
        ),
        "source_had_sex_field": source_had_sex_field,
        "validated_sample_path": str(CHR_X_STEP2_SAMPLE),
        "validated_sample_sha256": sha256_file(
            CHR_X_STEP2_SAMPLE
        ),
        "n_source_sample_rows": int(sample.shape[0]),
        "n_analysis_ids": int(frozen_sex.shape[0]),
        "n_males_plink_code_1": observed_counts[1],
        "n_females_plink_code_2": observed_counts[2],
        "all_analysis_ids_present": True,
        "all_analysis_ids_aligned_exactly": True,
        "all_analysis_sex_codes_in_1_2": True,
        "all_analysis_sex_codes_match_frozen_mapping": True,
        "nonanalysis_rows_sex_code": 0,
        "frozen_id_order_sha256": stable_fingerprint(
            aligned_ids.values.tolist()
        ),
        "source_row_position_order_sha256": stable_fingerprint(
            source_positions.tolist()
        ),
        "frozen_id_sex_mapping_sha256": stable_fingerprint(
            [
                [fid, iid, int(sex_code)]
                for (fid, iid), sex_code in zip(
                    aligned_ids[["FID", "IID"]].values.tolist(),
                    aligned_sex.tolist(),
                )
            ]
        ),
    }
    if (
        validation["frozen_id_order_sha256"]
        != setup_summary["frozen_id_order_sha256"]
    ):
        raise RuntimeError(
            "Chromosome X SAMPLE validation produced an unexpected "
            "frozen-ID order fingerprint."
        )
    atomic_write_json(CHR_X_SAMPLE_SEX_QC_PATH, validation)
    return validation


chrX_sample_sex_validation = build_validated_chrX_sample()

CHR_X_STEP2_PGEN_PREFIX = (
    PRIVATE / "chrX_step2_hardcall_pgen_sex_validated"
)
CHR_X_STEP2_PGEN_SENTINEL = Path(
    str(CHR_X_STEP2_PGEN_PREFIX) + ".complete.json"
)
CHR_X_PAR1_END_GRCH38 = 2781479
CHR_X_PAR2_START_GRCH38 = 155701383
EXPECTED_CHRX_LEADS = {
    "chrX:22917619:C:T": {
        "position": 22917619,
        "ref": "C",
        "alt": "T",
    },
    "chrX:5998724:A:G": {
        "position": 5998724,
        "ref": "A",
        "alt": "G",
    },
    "chrX:67797018:T:A": {
        "position": 67797018,
        "ref": "T",
        "alt": "A",
    },
}


def read_chrX_psam(path):
    psam = pd.read_csv(
        path,
        sep=r"\s+",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    fid_column = "#FID" if "#FID" in psam.columns else "FID"
    if (
        fid_column not in psam.columns
        or "IID" not in psam.columns
        or "SEX" not in psam.columns
    ):
        raise RuntimeError(
            f"Chromosome X PSAM lacks FID/IID/SEX: {path}"
        )
    return (
        psam[[fid_column, "IID", "SEX"]]
        .rename(columns={fid_column: "FID"})
        .reset_index(drop=True)
    )


def validate_chrX_step2_pgen(
    prefix,
    expected_n_variants=None,
):
    prefix = Path(prefix)
    pgen = Path(str(prefix) + ".pgen")
    pvar = Path(str(prefix) + ".pvar")
    psam = Path(str(prefix) + ".psam")
    if not all(nonempty(path) for path in [pgen, pvar, psam]):
        raise RuntimeError(
            "Chromosome X hard-call PGEN trio is incomplete."
        )

    observed_psam = read_chrX_psam(psam)
    if (
        observed_psam.shape[0] != 392155
        or observed_psam[["FID", "IID"]].duplicated().any()
    ):
        raise RuntimeError(
            "Chromosome X PSAM is not a unique 392,155-row cohort."
        )

    frozen_sex = pd.read_csv(
        FROZEN_COVAR,
        sep="\t",
        dtype={"FID": str, "IID": str},
    )[["FID", "IID", "sex"]]
    frozen_sex_numeric = pd.to_numeric(
        frozen_sex["sex"],
        errors="raise",
    )
    if not (
        np.isfinite(
            frozen_sex_numeric.to_numpy(dtype=float)
        ).all()
        and frozen_sex_numeric.isin([0, 1]).all()
    ):
        raise RuntimeError(
            "Frozen sex covariate is not exact finite 0/1 during "
            "chromosome X PGEN validation."
        )
    expected_ids = frozen_sex[["FID", "IID"]].reset_index(drop=True)
    if not observed_psam[["FID", "IID"]].equals(expected_ids):
        raise RuntimeError(
            "Chromosome X PSAM IDs do not match the frozen cohort "
            "in exact order."
        )

    expected_sex = np.where(
        frozen_sex_numeric.to_numpy() == 1,
        1,
        2,
    )
    observed_sex_numeric = pd.to_numeric(
        observed_psam["SEX"],
        errors="raise",
    )
    if not (
        np.isfinite(
            observed_sex_numeric.to_numpy(dtype=float)
        ).all()
        and observed_sex_numeric.isin([1, 2]).all()
    ):
        raise RuntimeError(
            "Chromosome X PSAM sex is not exact finite PLINK 1/2."
        )
    observed_sex = observed_sex_numeric.to_numpy(dtype=int)
    if not np.array_equal(observed_sex, expected_sex):
        raise RuntimeError(
            "Chromosome X PSAM sex codes do not exactly match the "
            "frozen mapping."
        )

    pvar_frame = read_pvar_rows(pvar)
    if pvar_frame.empty or pvar_frame["ID"].duplicated().any():
        raise RuntimeError(
            "Chromosome X PVAR is empty or has duplicate variant IDs."
        )
    positions_numeric = pd.to_numeric(
        pvar_frame["POS"],
        errors="raise",
    )
    if not (
        np.isfinite(positions_numeric.to_numpy(dtype=float)).all()
        and (positions_numeric > 0).all()
        and np.equal(
            positions_numeric,
            np.floor(positions_numeric),
        ).all()
    ):
        raise RuntimeError(
            "Chromosome X PVAR has invalid base-pair positions."
        )
    pvar_frame["POS"] = positions_numeric.astype(int)
    if (
        pvar_frame["REF"].eq("").any()
        or pvar_frame["ALT"].eq("").any()
    ):
        raise RuntimeError(
            "Chromosome X PVAR has missing REF/ALT alleles."
        )
    if (
        expected_n_variants is not None
        and pvar_frame.shape[0] != int(expected_n_variants)
    ):
        raise RuntimeError(
            "Chromosome X PVAR variant count differs from the source "
            f"BGEN header: {pvar_frame.shape[0]:,} versus "
            f"{int(expected_n_variants):,}."
        )

    for variant_id, expected in EXPECTED_CHRX_LEADS.items():
        lead_rows = pvar_frame.loc[pvar_frame["ID"] == variant_id]
        if lead_rows.shape[0] != 1:
            raise RuntimeError(
                f"Expected chromosome X lead is not unique: {variant_id}"
            )
        lead = lead_rows.iloc[0]
        if not (
            int(lead["POS"]) == expected["position"]
            and str(lead["REF"]) == expected["ref"]
            and str(lead["ALT"]) == expected["alt"]
        ):
            raise RuntimeError(
                "Chromosome X lead REF/ALT semantics changed during "
                f"conversion: {variant_id}"
            )

    chromosome_codes = (
        pvar_frame["CHROM"].astype(str).str.upper()
    )
    par1 = pvar_frame["POS"] <= CHR_X_PAR1_END_GRCH38
    par2 = pvar_frame["POS"] >= CHR_X_PAR2_START_GRCH38
    nonpar = ~(par1 | par2)
    x_codes = {"X", "23"}
    par_codes = {"PAR1", "PAR2", "XY", "25"}
    if not set(chromosome_codes.loc[nonpar]).issubset(x_codes):
        raise RuntimeError(
            "Non-PAR chromosome X variants have unexpected PVAR "
            "chromosome codes."
        )
    if (
        (par1 | par2).any()
        and not set(
            chromosome_codes.loc[par1 | par2]
        ).issubset(par_codes)
    ):
        raise RuntimeError(
            "GRCh38 PAR variants were not represented with PAR/XY "
            "chromosome codes in the chromosome X PVAR."
        )

    pvar_order_hash = stable_fingerprint(
        pvar_frame[
            ["CHROM", "POS", "ID", "REF", "ALT"]
        ].values.tolist()
    )
    frozen_id_sex_mapping_hash = stable_fingerprint(
        [
            [fid, iid, int(sex_code)]
            for (fid, iid), sex_code in zip(
                expected_ids.values.tolist(),
                observed_sex.tolist(),
            )
        ]
    )
    if (
        frozen_id_sex_mapping_hash
        != chrX_sample_sex_validation[
            "frozen_id_sex_mapping_sha256"
        ]
    ):
        raise RuntimeError(
            "Chromosome X PSAM ID-to-sex mapping differs from the "
            "validated Oxford SAMPLE mapping."
        )
    return {
        "pgen_prefix": str(prefix),
        "n_samples": int(observed_psam.shape[0]),
        "n_variants": int(pvar_frame.shape[0]),
        "n_males_plink_code_1": int((observed_sex == 1).sum()),
        "n_females_plink_code_2": int((observed_sex == 2).sum()),
        "psam_ids_match_frozen_order": True,
        "psam_sex_matches_frozen_mapping": True,
        "frozen_id_sex_mapping_sha256": (
            frozen_id_sex_mapping_hash
        ),
        "all_three_frozen_chrX_leads_match_ref_alt": True,
        "n_par1_variants": int(par1.sum()),
        "n_par2_variants": int(par2.sum()),
        "n_nonpar_variants": int(nonpar.sum()),
        "par_disposition_grch38_validated": True,
        "pvar_variant_order_sha256": pvar_order_hash,
        "pgen": file_metadata(pgen, include_sha256=True),
        "pvar": file_metadata(pvar, include_sha256=True),
        "psam": file_metadata(psam, include_sha256=True),
    }


def run_chrX_pgen_integrity_checks(prefix, expected_n_variants):
    validate_log = (
        LOG_DIR / "validate_chrX_step2_hardcall_pgen.log"
    )
    pgen_info_log = (
        LOG_DIR / "inspect_chrX_step2_hardcall_pgen.log"
    )
    run_checked(
        [
            "plink2",
            "--pfile", str(prefix),
            "--validate",
            "--out", str(
                LOG_DIR / "chrX_step2_hardcall_validate"
            ),
        ],
        "Validate chromosome X hard-call PGEN structure",
        validate_log,
    )
    pgen_info_text = run_checked(
        [
            "plink2",
            "--pfile", str(prefix),
            "--pgen-info",
            "--out", str(
                LOG_DIR / "chrX_step2_hardcall_pgen_info"
            ),
        ],
        "Confirm chromosome X PGEN contains no dosages",
        pgen_info_log,
    )
    if re.search(
        r"(?im)^\s*No dosages present\s*$",
        pgen_info_text,
    ) is None:
        raise RuntimeError(
            "PLINK --pgen-info did not confirm 'No dosages present'."
        )
    variant_match = re.search(
        r"(?im)^\s*Variants:\s*(\d+)\s*$",
        pgen_info_text,
    )
    sample_match = re.search(
        r"(?im)^\s*Samples:\s*(\d+)\s*$",
        pgen_info_text,
    )
    if (
        variant_match is None
        or int(variant_match.group(1)) != int(expected_n_variants)
        or sample_match is None
        or int(sample_match.group(1)) != 392155
    ):
        raise RuntimeError(
            "PLINK --pgen-info dimensions do not match the validated "
            "chromosome X PGEN."
        )
    return {
        "plink_validate_passed": True,
        "pgen_info_no_dosages": True,
        "validate_log_sha256": sha256_file(validate_log),
        "pgen_info_log_sha256": sha256_file(pgen_info_log),
    }


def run_chrX_pgen_conversion(command, log_path):
    command = [str(value) for value in command]
    log_path = Path(log_path)

    output_prefix = None

    if "--out" in command:
        out_index = command.index("--out")

        if out_index + 1 < len(command):
            output_prefix = Path(
                command[out_index + 1]
            )

    emergency_floor_bytes = int(
        float(
            PRODUCTION_CONFIG["runtime_floor_gib"]
        )
        * 1024**3
    )

    stop_trigger_bytes = (
        emergency_floor_bytes
        + 10 * 1024**3
    )

    initial_free_bytes = shutil.disk_usage(
        "/home/jupyter"
    ).free

    if initial_free_bytes <= stop_trigger_bytes:
        raise RuntimeError(
            "Chromosome X conversion cannot start safely: "
            f"{initial_free_bytes / 1024**3:.1f} GiB free."
        )

    def clear_monitored_partial():
        if output_prefix is None:
            return

        for suffix in [
            ".pgen",
            ".pgen.pgi",
            ".pvar",
            ".psam",
            ".log",
        ]:
            Path(
                str(output_prefix) + suffix
            ).unlink(missing_ok=True)

    with log_path.open("w") as log_handle:
        process = subprocess.Popen(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

        last_state_update = 0.0

        while process.poll() is None:
            time.sleep(5)

            free_bytes = shutil.disk_usage(
                "/home/jupyter"
            ).free

            now = time.monotonic()

            if now - last_state_update >= 30:
                production_state_update(
                    status="running",
                    stage="step2",
                    current_chromosome="X",
                    current_stage_detail=(
                        "building_disk_guarded_"
                        "chrX_hardcall_pgen"
                    ),
                    chrX_conversion_free_gib=round(
                        free_bytes / 1024**3,
                        2,
                    ),
                )

                last_state_update = now

            if free_bytes < stop_trigger_bytes:
                process.terminate()

                try:
                    process.wait(timeout=30)

                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=30)

                log_handle.flush()
                clear_monitored_partial()

                raise RuntimeError(
                    "Chromosome X conversion was stopped "
                    "before the disk reserve could be breached. "
                    "Partial PGEN files were removed. "
                    f"Free space at stop: "
                    f"{free_bytes / 1024**3:.1f} GiB."
                )

        return_code = process.returncode

    combined = log_path.read_text(
        errors="replace"
    )

    print(combined)

    return return_code, combined


def build_chrX_step2_pgen():
    base_command = [
        "plink2",
        "--bgen", qc_sources["X"]["bgen"], "ref-last",
        "--sample", str(CHR_X_STEP2_SAMPLE),
        "--keep", str(KEEP_FILE),
        "--indiv-sort", "file", str(KEEP_FILE),
        "--hard-call-threshold", "0.1",
        "--make-pgen", "erase-dosage",
        "--threads", "24",
    ]
    base_specification = {
        "qc_source_fingerprint": qc_sources["X"][
            "source_fingerprint"
        ],
        "validated_sample_sha256": sha256_file(
            CHR_X_STEP2_SAMPLE
        ),
        "validated_sample_sex_fingerprint": stable_fingerprint(
            chrX_sample_sex_validation
        ),
        "frozen_input_hashes": current_frozen_input_hashes,
        "frozen_id_order_sha256": setup_summary[
            "frozen_id_order_sha256"
        ],
        "plink2_version": PLINK2_VERSION,
        "base_command": base_command + [
            "--out", "<OUTPUT_PREFIX>"
        ],
        "hard_call_threshold": 0.1,
        "dosage_policy": "erase-dosage",
        "input_allele_semantics": (
            "The retained QC BGEN was exported from PLINK BED with "
            "REF stored last; ref-last restores ALT/REF semantics. "
            "No variant filter, call-changing fill operation, or "
            "invalid-haploid erasure is applied."
        ),
    }

    if nonempty(CHR_X_STEP2_PGEN_SENTINEL):
        try:
            record = json.loads(
                CHR_X_STEP2_PGEN_SENTINEL.read_text()
            )
            recorded_split = record["split_par_applied"]
            specification = {
                **base_specification,
                "split_par_applied": recorded_split,
            }
            fingerprint = stable_fingerprint(specification)
            if not sentinel_matches(
                CHR_X_STEP2_PGEN_SENTINEL,
                fingerprint,
            ):
                raise RuntimeError(
                    "Chromosome X PGEN sentinel fingerprint changed."
                )
            validation = validate_chrX_step2_pgen(
                CHR_X_STEP2_PGEN_PREFIX,
                expected_n_variants=record[
                    "source_bgen_variant_count"
                ],
            )
            if not all(
                recorded_file_matches(
                    Path(
                        str(CHR_X_STEP2_PGEN_PREFIX) + suffix
                    ),
                    record.get(suffix.lstrip(".")),
                    include_sha256=True,
                )
                for suffix in [".pgen", ".pvar", ".psam"]
            ):
                raise RuntimeError(
                    "Chromosome X PGEN metadata differs from the "
                    "completion record."
                )
            pgi_path = Path(
                str(CHR_X_STEP2_PGEN_PREFIX) + ".pgen.pgi"
            )
            recorded_pgi = record.get("pgi")
            if recorded_pgi:
                if not recorded_file_matches(
                    pgi_path,
                    recorded_pgi,
                    include_sha256=True,
                ):
                    raise RuntimeError(
                        "Chromosome X PGEN index differs from the "
                        "completion record."
                    )
            elif pgi_path.exists():
                raise RuntimeError(
                    "Unexpected chromosome X PGEN index is present."
                )
            integrity = run_chrX_pgen_integrity_checks(
                CHR_X_STEP2_PGEN_PREFIX,
                record["source_bgen_variant_count"],
            )
            return {
                **record,
                **validation,
                **integrity,
                "status": "existing_validated_hardcall_pgen",
            }
        except Exception as exc:
            print(
                "Existing chromosome X hard-call PGEN is not reusable:",
                exc,
            )

    partial_prefix = Path(
        str(CHR_X_STEP2_PGEN_PREFIX) + ".partial"
    )

    def clear_partial():
        for suffix in [
            ".pgen", ".pgen.pgi", ".pvar", ".psam", ".log"
        ]:
            Path(str(partial_prefix) + suffix).unlink(
                missing_ok=True
            )

    clear_partial()
    initial_log = (
        LOG_DIR / "build_chrX_step2_hardcall_pgen_initial.log"
    )
    initial_command = base_command + [
        "--out", str(partial_prefix)
    ]
    print(
        "\n[Build chromosome X hard-call PGEN without redundant "
        "PAR splitting]"
    )
    return_code, conversion_text = run_chrX_pgen_conversion(
        initial_command,
        initial_log,
    )
    split_par_applied = False
    actual_command = initial_command

    if return_code != 0:
        explicit_unsplit_par_error = re.search(
            r"(?is)pseudoautosomal variant.*--split-par|"
            r"--split-par.*pseudoautosomal variant",
            conversion_text,
        )
        if explicit_unsplit_par_error is None:
            raise RuntimeError(
                "Chromosome X PGEN conversion failed for a reason "
                f"other than unsplit PAR metadata; see {initial_log}."
            )
        clear_partial()
        split_par_applied = True
        actual_command = base_command + [
            "--split-par", "hg38",
            "--out", str(partial_prefix),
        ]
        split_log = (
            LOG_DIR
            / "build_chrX_step2_hardcall_pgen_split_par.log"
        )
        print(
            "\n[Retry chromosome X hard-call PGEN after explicit "
            "unsplit-PAR error]"
        )
        return_code, conversion_text = run_chrX_pgen_conversion(
            actual_command,
            split_log,
        )
        if return_code != 0:
            raise RuntimeError(
                "Chromosome X PGEN conversion failed after the "
                f"explicit PAR split; see {split_log}."
            )

    source_variant_match = re.search(
        r"--bgen:\s*(\d+)\s+variants\s+"
        r"(?:declared in header|detected)",
        conversion_text,
    )
    if source_variant_match is None:
        raise RuntimeError(
            "Could not recover the source BGEN variant count from "
            "the successful PLINK conversion log."
        )
    source_bgen_variant_count = int(
        source_variant_match.group(1)
    )
    partial_validation = validate_chrX_step2_pgen(
        partial_prefix,
        expected_n_variants=source_bgen_variant_count,
    )

    for suffix in [".pgen", ".pvar", ".psam", ".log"]:
        partial_path = Path(str(partial_prefix) + suffix)
        final_path = Path(
            str(CHR_X_STEP2_PGEN_PREFIX) + suffix
        )
        if suffix != ".log" or nonempty(partial_path):
            if not nonempty(partial_path):
                raise RuntimeError(
                    "Missing chromosome X partial PGEN artifact: "
                    f"{partial_path}"
                )
            os.replace(partial_path, final_path)
    partial_pgi = Path(str(partial_prefix) + ".pgen.pgi")
    final_pgi = Path(
        str(CHR_X_STEP2_PGEN_PREFIX) + ".pgen.pgi"
    )
    if nonempty(partial_pgi):
        os.replace(partial_pgi, final_pgi)
    else:
        final_pgi.unlink(missing_ok=True)

    validation = validate_chrX_step2_pgen(
        CHR_X_STEP2_PGEN_PREFIX,
        expected_n_variants=source_bgen_variant_count,
    )
    integrity = run_chrX_pgen_integrity_checks(
        CHR_X_STEP2_PGEN_PREFIX,
        source_bgen_variant_count,
    )
    specification = {
        **base_specification,
        "split_par_applied": split_par_applied,
    }
    fingerprint = stable_fingerprint(specification)
    record = {
        **specification,
        **validation,
        **integrity,
        "fingerprint": fingerprint,
        "actual_conversion_command": actual_command,
        "source_bgen_variant_count": source_bgen_variant_count,
        "pgen": file_metadata(
            Path(str(CHR_X_STEP2_PGEN_PREFIX) + ".pgen"),
            include_sha256=True,
        ),
        "pvar": file_metadata(
            Path(str(CHR_X_STEP2_PGEN_PREFIX) + ".pvar"),
            include_sha256=True,
        ),
        "psam": file_metadata(
            Path(str(CHR_X_STEP2_PGEN_PREFIX) + ".psam"),
            include_sha256=True,
        ),
        "pgi": (
            file_metadata(final_pgi, include_sha256=True)
            if nonempty(final_pgi)
            else None
        ),
    }
    atomic_write_json(CHR_X_STEP2_PGEN_SENTINEL, record)
    return {
        **record,
        "status": "completed_now",
    }


chrX_hardcall_validation = build_chrX_step2_pgen()
atomic_write_json(
    AGG_DIR / "corrected_chrX_hardcall_pgen_validation.json",
    chrX_hardcall_validation,
)

STEP2_SAMPLE_PATHS = {
    chromosome: Path(qc_sources[chromosome]["sample"])
    for chromosome in CHROMOSOMES
}
STEP2_SAMPLE_PATHS["X"] = CHR_X_STEP2_SAMPLE
STEP2_SAMPLE_SHA256 = {
    chromosome: sha256_file(path)
    for chromosome, path in STEP2_SAMPLE_PATHS.items()
}

MAX_PARALLEL_CHROMOSOMES = 4
THREADS_PER_CHROMOSOME = 12


def corrected_result_path(chromosome):
    prefix = STEP2_DIR / f"chr{chromosome}_pso_loco_corrected"
    return Path(str(prefix) + "_psoriasis_case.regenie")


def step2_specification(chromosome):
    step2_sample_path = STEP2_SAMPLE_PATHS[chromosome]
    if chromosome == "X":
        genotype_arguments = [
            "--pgen", str(CHR_X_STEP2_PGEN_PREFIX),
        ]
        input_mode = "validated_hardcall_pgen_without_dosages"
    else:
        genotype_arguments = [
            "--bgen", qc_sources[chromosome]["bgen"],
            "--sample", str(step2_sample_path),
        ]
        input_mode = "retained_qc_bgen"

    command = [
        "regenie",
        "--step", "2",
        "--bt",
        *genotype_arguments,
        "--phenoFile", str(FROZEN_PHENO),
        "--phenoCol", "psoriasis_case",
        "--covarFile", str(FROZEN_COVAR),
        "--bsize", "2000",
        "--firth",
        "--approx",
        "--pThresh", "0.01",
        "--pred", str(PRED_LIST),
        "--af-cc",
        "--threads", str(THREADS_PER_CHROMOSOME),
        "--out", "<OUTPUT_PREFIX>",
    ]
    # Do not add REGENIE --ref-first here. The historical QC BGENs were
    # exported from PLINK BED with REF stored last. REGENIE's default
    # interpretation therefore reports the ALT/tested allele as ALLELE1,
    # matching the frozen primary manifest. Block R6 verifies this directly.
    return {
        "chromosome": chromosome,
        "qc_source_fingerprint": qc_sources[chromosome][
            "source_fingerprint"
        ],
        "frozen_input_hashes": setup_summary[
            "frozen_input_hashes"
        ],
        "step1_fingerprint": step1_model_validation[
            "step1_fingerprint"
        ],
        "prediction_list_sha256": sha256_file(PRED_LIST),
        "loco_sha256": loco_validation["loco_sha256"],
        "step2_sample_path": str(step2_sample_path),
        "step2_sample_sha256": STEP2_SAMPLE_SHA256[chromosome],
        "chrX_sample_sex_validation_fingerprint": (
            stable_fingerprint(chrX_sample_sex_validation)
            if chromosome == "X"
            else None
        ),
        "chrX_hardcall_pgen_fingerprint": (
            chrX_hardcall_validation["fingerprint"]
            if chromosome == "X"
            else None
        ),
        "genotype_input_mode": input_mode,
        "regenie_version": REGENIE_VERSION,
        "command": command,
        "allele_semantics": (
            "Autosomal Step 2 retains the historical BGEN default. "
            "Chromosome X is recreated from that exact hard-call BGEN "
            "with ref-last and analyzed through a dosage-free PGEN so "
            "REGENIE can "
            "apply male X handling. REGENIE ALLELE1 is validated "
            "against the frozen 21-lead manifest in Block R6."
        ),
    }


def validate_step2_output(
    chromosome,
    result_path,
    stdout_log_path,
    expected_logged_result_path=None,
):
    result_path = Path(result_path)
    stdout_log_path = Path(stdout_log_path)
    if not nonempty(result_path) or not nonempty(stdout_log_path):
        raise RuntimeError(
            f"Missing chr{chromosome} result or stdout log."
        )
    log_text = stdout_log_path.read_text(errors="replace")
    fatal = re.search(
        r"(?im)\b(?:ERROR:|FATAL|segmentation fault|"
        r"terminate called|killed|std::bad_alloc)\b",
        log_text,
    )
    if fatal:
        raise RuntimeError(
            f"Fatal chr{chromosome} Step 2 log text detected: "
            f"{fatal.group(0)}"
        )
    if chromosome == "X" and re.search(
        r"(?im)genotype is 1 for a male on chrX|"
        r"male.*chrX.*(?:diploid|0/2)",
        log_text,
    ):
        raise RuntimeError(
            "Chromosome X Step 2 emitted a male diploid-coding "
            "warning. Do not use this result."
        )

    required_patterns = [
        (
            r"(?m)^\s*-number of phenotyped individuals with no "
            r"missing data\s*=\s*392155\s*$"
        ),
        (
            r"(?m)^\s*-number of individuals with covariate data"
            r"\s*=\s*392155\s*$"
        ),
        (
            r"(?m)^\s*-number of individuals used in analysis"
            r"\s*=\s*392155\s*$"
        ),
        (
            r"(?m)^\s*-\s*'psoriasis_case':\s*5805 cases and "
            r"386350 controls\s*$"
        ),
        (
            r"(?m)^Association results stored separately for each "
            r"trait in files\s*:\s*$"
        ),
        r"(?m)^\s*Elapsed time\s*:\s*[0-9.]+s\s*$",
        r"(?m)^\s*End time:\s*.+$",
    ]
    if expected_logged_result_path is not None:
        required_patterns.append(
            rf"(?m)^\s*[*-]\s*\["
            rf"{re.escape(str(expected_logged_result_path))}"
            rf"\]\s*$"
        )
    missing_patterns = [
        pattern
        for pattern in required_patterns
        if re.search(pattern, log_text) is None
    ]
    if missing_patterns:
        raise RuntimeError(
            f"chr{chromosome} Step 2 log is missing required "
            "completion evidence:\n"
            + "\n".join(missing_patterns)
        )

    firth_match = re.search(
        r"(?m)^Number of tests with Firth correction\s*:\s*"
        r"(\d+)\s*$",
        log_text,
    )
    if firth_match is None:
        raise RuntimeError(
            f"chr{chromosome} log lacks the Firth-test count."
        )
    n_firth = int(firth_match.group(1))
    failed_match = re.search(
        r"(?m)^Number of failed tests\s*:\s*"
        r"\((\d+)/(\d+)\)\s*$",
        log_text,
    )
    if n_firth > 0:
        if failed_match is None:
            raise RuntimeError(
                f"chr{chromosome} log lacks the failed-test summary."
            )
        n_failed, failure_denominator = map(
            int,
            failed_match.groups(),
        )
        if n_failed != 0 or failure_denominator != n_firth:
            raise RuntimeError(
                f"chr{chromosome} reported {n_failed}/"
                f"{failure_denominator} failed Firth tests; "
                f"expected 0/{n_firth}."
            )
    elif failed_match is not None:
        n_failed, failure_denominator = map(
            int,
            failed_match.groups(),
        )
        if n_failed != 0 or failure_denominator != 0:
            raise RuntimeError(
                f"Unexpected chr{chromosome} failure summary: "
                f"{n_failed}/{failure_denominator}."
            )

    # REGENIE reports the input variant count differently for PGEN and
    # BGEN.  PGEN logs use "n_snps = N", while BGEN logs use a summary
    # such as "... with S named samples and N variants ...".
    pgen_count_matches = re.findall(
        r"\bn_snps\s*=\s*([\d,]+)",
        log_text,
    )
    bgen_count_matches = re.findall(
        r"(?im)^\s*[-+*]?\s*summary\s*:\s*bgen file\b"
        r"[^\n]*?\band\s+([\d,]+)\s+variants\b",
        log_text,
    )
    if chromosome == "X":
        if len(pgen_count_matches) != 1:
            raise RuntimeError(
                "chrX PGEN log must report exactly one n_snps count; "
                f"observed {pgen_count_matches}."
            )
        if bgen_count_matches:
            raise RuntimeError(
                "chrX Step 2 unexpectedly reported a BGEN input summary."
            )
        expected_rows = int(
            pgen_count_matches[0].replace(",", "")
        )
        input_count_source = "pgen_n_snps"
    else:
        if len(bgen_count_matches) != 1:
            raise RuntimeError(
                f"chr{chromosome} BGEN log must report exactly one "
                "input variant count; observed "
                f"{bgen_count_matches}."
            )
        if pgen_count_matches:
            raise RuntimeError(
                f"chr{chromosome} Step 2 unexpectedly reported a "
                "PGEN/BED n_snps input."
            )
        expected_rows = int(
            bgen_count_matches[0].replace(",", "")
        )
        input_count_source = "bgen_summary_variants"

    with result_path.open() as handle:
        header = handle.readline().split()
        n_rows = sum(1 for line in handle if line.strip())
    required_columns = {
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
        "A1FREQ", "A1FREQ_CASES", "A1FREQ_CONTROLS",
        "N", "N_CASES", "N_CONTROLS", "TEST",
        "BETA", "SE", "CHISQ", "LOG10P", "EXTRA",
    }
    missing_columns = required_columns.difference(header)
    if missing_columns:
        raise RuntimeError(
            f"chr{chromosome} result is missing columns: "
            f"{sorted(missing_columns)}"
        )
    if n_rows != expected_rows or n_rows == 0:
        raise RuntimeError(
            f"chr{chromosome} result has {n_rows:,} rows; "
            f"the REGENIE input log reports {expected_rows:,} SNPs."
        )
    return {
        "chromosome": chromosome,
        "n_result_rows": n_rows,
        "n_input_snps_from_log": expected_rows,
        "input_variant_count_source": input_count_source,
        "n_firth_tests": n_firth,
        "n_failed_tests": 0,
        "stdout_log_sha256": sha256_file(stdout_log_path),
        "result_sha256": sha256_file(result_path),
        "successful_end_marker": True,
        "chrX_male_diploid_warning_absent": (
            True if chromosome == "X" else None
        ),
    }


def step2_result_is_complete(chromosome):
    result_path = corrected_result_path(chromosome)
    log_path = LOG_DIR / f"corrected_step2_chr{chromosome}.log"
    sentinel = Path(
        str(
            STEP2_DIR
            / f"chr{chromosome}_pso_loco_corrected"
        )
        + ".complete.json"
    )
    spec = step2_specification(chromosome)
    fingerprint = stable_fingerprint(spec)
    if not sentinel_matches(sentinel, fingerprint):
        return False
    try:
        record = json.loads(sentinel.read_text())
        if not (
            recorded_file_matches(
                result_path,
                record.get("result"),
                include_sha256=True,
            )
            and recorded_file_matches(
                log_path,
                record.get("stdout_log"),
                include_sha256=True,
            )
        ):
            return False
        validate_step2_output(
            chromosome,
            result_path,
            log_path,
            expected_logged_result_path=record[
                "logged_result_path"
            ],
        )
    except Exception as exc:
        print(
            f"Existing chr{chromosome} Step 2 is not reusable:",
            exc,
        )
        return False
    return True


def run_corrected_step2(chromosome):
    if step2_result_is_complete(chromosome):
        return {
            "chromosome": chromosome,
            "status": "existing_complete_output",
            "result": str(corrected_result_path(chromosome)),
        }

    source = qc_sources[chromosome]
    output_prefix = (
        STEP2_DIR / f"chr{chromosome}_pso_loco_corrected"
    )
    log_path = LOG_DIR / f"corrected_step2_chr{chromosome}.log"
    sentinel = Path(str(output_prefix) + ".complete.json")
    spec = step2_specification(chromosome)
    fingerprint = stable_fingerprint(spec)
    partial_prefix = (
        STEP2_DIR / f"chr{chromosome}_pso_loco_corrected.partial"
    )
    partial_result = Path(
        str(partial_prefix) + "_psoriasis_case.regenie"
    )
    partial_log = Path(str(log_path) + ".partial")
    for path in [
        partial_result,
        Path(str(partial_prefix) + ".log"),
        partial_log,
    ]:
        path.unlink(missing_ok=True)

    command = [
        str(partial_prefix)
        if value == "<OUTPUT_PREFIX>"
        else value
        for value in spec["command"]
    ]

    with partial_log.open("w") as log_handle:
        log_handle.write(
            "[Corrected chromosome-specific Step 2]\n"
            + " ".join(command)
            + "\n\n"
        )
        completed = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if completed.returncode != 0:
        raise RuntimeError(
            f"Corrected Step 2 failed for chr{chromosome}; "
            f"see {partial_log}."
        )
    validation = validate_step2_output(
        chromosome,
        partial_result,
        partial_log,
        expected_logged_result_path=partial_result,
    )
    os.replace(partial_result, corrected_result_path(chromosome))
    partial_regenie_log = Path(str(partial_prefix) + ".log")
    final_regenie_log = Path(str(output_prefix) + ".log")
    if nonempty(partial_regenie_log):
        os.replace(partial_regenie_log, final_regenie_log)
    os.replace(partial_log, log_path)
    atomic_write_json(
        sentinel,
        {
            **spec,
            "fingerprint": fingerprint,
            "logged_result_path": str(partial_result),
            "validation": validation,
            "result": file_metadata(
                corrected_result_path(chromosome),
                include_sha256=True,
            ),
            "stdout_log": file_metadata(
                log_path,
                include_sha256=True,
            ),
        },
    )
    if not step2_result_is_complete(chromosome):
        raise RuntimeError(
            f"Promoted chr{chromosome} Step 2 output did not pass "
            "fingerprint and completion validation."
        )

    return {
        "chromosome": chromosome,
        "status": "completed_now",
        "result": str(corrected_result_path(chromosome)),
    }


step2_records = []
completed_step2 = []
pending_step2 = []

for chromosome in CHROMOSOMES:
    if step2_result_is_complete(chromosome):
        record = {
            "chromosome": chromosome,
            "status": "existing_complete_output",
            "result": str(
                corrected_result_path(chromosome)
            ),
        }
        step2_records.append(record)
        completed_step2.append(chromosome)
    else:
        pending_step2.append(chromosome)

production_state_update(
    status="running",
    stage="step2",
    current_chromosome=None,
    current_stage_detail="parallel_chromosome_jobs",
    completed_step2_chromosomes=completed_step2,
    active_step2_chromosomes=(
        pending_step2[:MAX_PARALLEL_CHROMOSOMES]
    ),
    pending_step2_chromosomes=pending_step2,
)

production_assert_disk_floor(
    "before parallel Step 2"
)

executor = concurrent.futures.ThreadPoolExecutor(
    max_workers=MAX_PARALLEL_CHROMOSOMES
)

futures = {
    executor.submit(
        run_corrected_step2,
        chromosome,
    ): chromosome
    for chromosome in pending_step2
}

try:
    for future in concurrent.futures.as_completed(
        futures
    ):
        chromosome = futures[future]
        record = future.result()
        step2_records.append(record)
        completed_step2.append(chromosome)

        remaining = [
            value
            for value in pending_step2
            if value not in completed_step2
        ]

        production_state_update(
            status="running",
            stage="step2",
            current_chromosome=None,
            current_stage_detail=(
                "parallel_chromosome_checkpoint_validated"
            ),
            completed_step2_chromosomes=(
                completed_step2
            ),
            active_step2_chromosomes=(
                remaining[
                    :MAX_PARALLEL_CHROMOSOMES
                ]
            ),
            pending_step2_chromosomes=remaining,
        )

        production_assert_disk_floor(
            f"after chr{chromosome} Step 2"
        )

        print(
            f"[chr{chromosome}] "
            f"{record['status']} -> "
            f"{record['result']}"
        )

except BaseException:
    for future in futures:
        future.cancel()

    executor.shutdown(
        wait=True,
        cancel_futures=True,
    )
    raise

else:
    executor.shutdown(wait=True)

production_state_update(
    status="running",
    stage="step2",
    current_chromosome=None,
    current_stage_detail=(
        "all_parallel_chromosomes_validated"
    ),
    completed_step2_chromosomes=completed_step2,
    active_step2_chromosomes=[],
    pending_step2_chromosomes=[],
)

step2_status = pd.DataFrame(step2_records)
step2_status["sort_order"] = step2_status["chromosome"].map(
    {
        **{str(value): value for value in range(1, 23)},
        "X": 23,
    }
)
step2_status = (
    step2_status
    .sort_values("sort_order")
    .drop(columns="sort_order")
    .reset_index(drop=True)
)

if step2_status.shape[0] != 23:
    raise RuntimeError(
        f"Expected 23 corrected Step 2 outputs; "
        f"observed {step2_status.shape[0]}."
    )
if not all(
    step2_result_is_complete(chromosome)
    for chromosome in CHROMOSOMES
):
    raise RuntimeError(
        "At least one corrected chromosome output is incomplete."
    )

step2_status.to_csv(
    AGG_DIR / "corrected_primary_step2_status.tsv",
    sep="\t",
    index=False,
)

chrX_log = LOG_DIR / "corrected_step2_chrX.log"
chrX_validation = validate_step2_output(
    "X",
    corrected_result_path("X"),
    chrX_log,
    expected_logged_result_path=json.loads(
        Path(
            str(
                STEP2_DIR / "chrX_pso_loco_corrected"
            )
            + ".complete.json"
        ).read_text()
    )["logged_result_path"],
)
atomic_write_json(
    AGG_DIR / "corrected_chrX_step2_coding_qc.json",
    {
        **chrX_validation,
        **chrX_sample_sex_validation,
        "hardcall_pgen_validation": chrX_hardcall_validation,
        "all_three_historical_leads_nonPAR": True,
        "entire_chrX_step2_log_checked": True,
        "sample_sex_field_positively_validated": True,
        "chrX_tested_through_validated_hardcall_pgen": True,
        "male_diploid_coding_warning_absent": True,
        "interpretation": (
            "All 392,155 analyzed chromosome X SAMPLE rows had exact "
            "1/2 sex codes matching the frozen covariate. The exact "
            "retained-QC BGEN was converted back to a validated "
            "dosage-free hard-call PGEN, its PSAM IDs and sex codes "
            "matched the frozen cohort, and chromosome X completed "
            "through "
            "REGENIE's hard-call path without a male diploid-coding "
            "warning. Warning absence is secondary to the positive "
            "SAMPLE, PSAM, PVAR, and dosage-free PGEN validation."
        ),
    },
)

print("\nAll 23 corrected primary Step 2 jobs passed.")
display(step2_status)
''',
]


def production_normalized_chromosome_values(series):
    return (
        series.astype(str)
        .str.strip()
        .str.upper()
        .str.replace(r"^CHR", "", regex=True)
    )


def production_upload_verified(local_path, remote_uri):
    local_path = Path(local_path)
    local_size = int(local_path.stat().st_size)
    local_md5 = production_md5_base64(local_path)
    describe = subprocess.run(
        [
            "gcloud", "storage", "objects", "describe",
            remote_uri,
            f"--billing-project={GOOGLE_PROJECT}",
            "--format=json",
        ],
        capture_output=True,
        text=True,
    )
    if describe.returncode != 0:
        copied = subprocess.run(
            [
                "gcloud", "storage", "cp", "--no-clobber",
                f"--billing-project={GOOGLE_PROJECT}",
                str(local_path),
                remote_uri,
            ],
            capture_output=True,
            text=True,
        )
        if copied.returncode != 0:
            raise RuntimeError(
                f"Upload failed for {remote_uri}:\n"
                + copied.stdout
                + copied.stderr
            )
        describe = subprocess.run(
            [
                "gcloud", "storage", "objects", "describe",
                remote_uri,
                f"--billing-project={GOOGLE_PROJECT}",
                "--format=json",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    metadata = json.loads(describe.stdout)
    remote_size = int(metadata["size"])
    remote_md5 = metadata.get("md5Hash")
    generation = str(metadata["generation"])
    if remote_size != local_size:
        raise RuntimeError(
            f"Remote byte count mismatch for {remote_uri}: "
            f"{remote_size} versus {local_size}."
        )
    if remote_md5 and remote_md5 != local_md5:
        raise RuntimeError(f"Remote MD5 mismatch for {remote_uri}.")
    return {
        "uri": remote_uri,
        "versioned_uri": f"{remote_uri}#{generation}",
        "generation": generation,
        "bytes": local_size,
        "sha256": production_sha256(local_path),
        "md5_base64": local_md5,
        "remote_md5_base64": remote_md5,
    }


def production_finalize(namespace):
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2

    production_state_update(
        status="running",
        stage="finalize",
        current_chromosome=None,
        current_stage_detail="streaming_and_validating_results",
    )
    production_assert_disk_floor("before final aggregation")

    CHROMOSOMES = namespace["CHROMOSOMES"]
    corrected_result_path = namespace["corrected_result_path"]
    AGG_DIR = namespace["AGG_DIR"]
    REPRO_DIR = namespace["REPRO_DIR"]
    STEP2_DIR = namespace["STEP2_DIR"]
    LOG_DIR = namespace["LOG_DIR"]
    REGENIE_VERSION = namespace["REGENIE_VERSION"]
    PLINK2_VERSION = namespace["PLINK2_VERSION"]
    step1_model_validation = json.loads(
        (AGG_DIR / "corrected_step1_model_validation.json").read_text()
    )
    chrX_full_qc = json.loads(
        (AGG_DIR / "corrected_chrX_step2_coding_qc.json").read_text()
    )

    REVIEW_DIR = ATTEMPT_ROOT / "review_package"
    SECURE_RESULTS_DIR = REPRO_DIR / "uncensored_summary_statistics"
    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    # A failed finalization may leave only whitelisted partial review outputs.
    # Recompute those deterministically in place, but refuse any unexpected
    # file and remove only the prior derived manifest before rebuilding it.
    prior_review_files = {
        path.name for path in REVIEW_DIR.iterdir() if path.is_file()
    }
    unexpected_prior = prior_review_files - EXPECTED_REVIEW_FILES
    if unexpected_prior:
        raise RuntimeError(
            "Unexpected file(s) in the immutable review directory: "
            f"{sorted(unexpected_prior)}"
        )
    (REVIEW_DIR / "aggregate_manifest.tsv").unlink(missing_ok=True)
    SECURE_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    expected_records = 13_255_796
    gwas_threshold = 5e-8
    log10_threshold = -math.log10(gwas_threshold)
    mhc_start = 25_000_000
    mhc_end = 34_000_000
    review_columns = [
        "source_chromosome", "position_bp", "ID", "ALLELE0",
        "ALLELE1", "TEST", "BETA", "SE", "CHISQ",
        "LOG10P_raw", "minus_log10_p",
        "p_value_derived_from_log10p",
        "p_value_numeric_floor_applied",
    ]
    required_columns = {
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
        "A1FREQ", "A1FREQ_CASES", "A1FREQ_CONTROLS",
        "N", "N_CASES", "N_CONTROLS", "TEST", "BETA",
        "SE", "CHISQ", "LOG10P", "EXTRA",
    }

    frozen_leads = [
        ("chr5:159402519:A:T", "T", 1.38, 1),
        ("chr19:10352442:G:C", "C", 0.66, -1),
        ("chr5:159260273:C:T", "T", 0.85, -1),
        ("chr5:132259137:C:CTGTGTG", "CTGTGTG", 1.14, 1),
        ("chr18:54306535:C:CT", "CT", 1.17, 1),
        ("chr5:132395634:A:G", "G", 1.14, 1),
        ("chr5:151087129:GC:G", "G", 1.27, 1),
        ("chr19:10475342:CTTT:C", "C", 0.72, -1),
        ("chr5:132122060:A:G", "G", 1.13, 1),
        ("chr18:54421047:C:T", "T", 1.15, 1),
        ("chr6:159068759:G:A", "A", 0.89, -1),
        ("chrX:22917619:C:T", "T", 1.38, 1),
        ("chr3:164387600:T:G", "G", 1.37, 1),
        ("chr17:80205094:C:T", "T", 0.90, -1),
        ("chr6:163803171:CTTATT:C", "C", 1.22, 1),
        ("chr5:131087866:G:A", "A", 1.12, 1),
        ("chrX:67797018:T:A", "A", 1.15, 1),
        ("chrX:5998724:A:G", "G", 1.49, 1),
        ("chr3:164138199:C:T", "T", 1.33, 1),
        ("chr5:132744811:G:A", "A", 1.15, 1),
        ("chr5:132497907:C:T", "T", 1.12, 1),
    ]
    lead_lookup = {
        variant_id: {
            "variant_id": variant_id,
            "expected_effect_allele": effect,
            "frozen_or": frozen_or,
            "expected_beta_sign": direction,
        }
        for variant_id, effect, frozen_or, direction in frozen_leads
    }
    observed_leads = {}
    chromosome_rows = []
    significant_parts = []
    top100 = pd.DataFrame()
    invalid_counts = {}
    all_chisq = []
    nonmhc_chisq = []
    total_records = 0
    total_valid = 0

    all_uncensored_path = (
        SECURE_RESULTS_DIR
        / "corrected_primary_all_output_records_UNCENSORED.tsv.gz"
    )
    valid_uncensored_path = (
        SECURE_RESULTS_DIR
        / "corrected_primary_valid_p_results_UNCENSORED.tsv.gz"
    )
    plotting_path = (
        REVIEW_DIR / "corrected_primary_plotting_data.tsv.gz"
    )
    with (
        gzip.open(all_uncensored_path, "wt") as all_handle,
        gzip.open(valid_uncensored_path, "wt") as valid_handle,
        gzip.open(plotting_path, "wt") as plotting_handle,
    ):
        wrote_all_header = False
        wrote_valid_header = False
        wrote_plot_header = False
        for chromosome in CHROMOSOMES:
            production_state_update(
                status="running",
                stage="finalize",
                current_chromosome=chromosome,
                current_stage_detail="streaming_result",
            )
            path = corrected_result_path(chromosome)
            chromosome_output = 0
            chromosome_valid = 0
            chromosome_invalid = 0
            chromosome_significant = 0
            chromosome_max = -math.inf
            chromosome_top = None
            previous_position = None
            previous_key = None

            for frame in pd.read_csv(
                path,
                sep=r"\s+",
                dtype={
                    "ID": str,
                    "ALLELE0": str,
                    "ALLELE1": str,
                    "TEST": str,
                    "LOG10P": str,
                    "EXTRA": str,
                },
                keep_default_na=False,
                na_filter=False,
                chunksize=250_000,
            ):
                missing = required_columns.difference(frame.columns)
                if missing:
                    raise RuntimeError(
                        f"chr{chromosome} result lacks {sorted(missing)}."
                    )
                frame["source_chromosome"] = chromosome
                normalized_chrom = production_normalized_chromosome_values(
                    frame["CHROM"]
                )
                if chromosome == "X":
                    allowed = {"X", "23", "PAR1", "PAR2", "XY", "25", "26"}
                else:
                    allowed = {chromosome}
                if not set(normalized_chrom).issubset(allowed):
                    raise RuntimeError(
                        f"chr{chromosome} output contains unexpected CHROM "
                        f"values: {sorted(set(normalized_chrom) - allowed)}"
                    )

                frame["position_bp"] = pd.to_numeric(
                    frame["GENPOS"], errors="raise"
                ).astype("int64")
                if (frame["position_bp"] <= 0).any():
                    raise RuntimeError(
                        f"chr{chromosome} contains non-positive positions."
                    )
                positions = frame["position_bp"].to_numpy()
                if previous_position is not None and positions[0] < previous_position:
                    raise RuntimeError(
                        f"chr{chromosome} output is not position-sorted."
                    )
                if (np.diff(positions) < 0).any():
                    raise RuntimeError(
                        f"chr{chromosome} output is not position-sorted."
                    )
                previous_position = int(positions[-1])

                if set(frame["TEST"].unique()) != {"ADD"}:
                    raise RuntimeError(
                        f"chr{chromosome} contains a non-ADD test."
                    )
                if (
                    frame["ID"].eq("").any()
                    or frame["ALLELE0"].eq("").any()
                    or frame["ALLELE1"].eq("").any()
                ):
                    raise RuntimeError(
                        f"chr{chromosome} contains blank ID/allele text."
                    )
                keys = list(
                    zip(
                        frame["position_bp"].astype(str),
                        frame["ID"],
                        frame["ALLELE0"],
                        frame["ALLELE1"],
                        frame["TEST"],
                    )
                )
                if previous_key is not None and keys[0] == previous_key:
                    raise RuntimeError(
                        f"chr{chromosome} has a duplicate key at a chunk boundary."
                    )
                if len(keys) != len(set(keys)):
                    raise RuntimeError(
                        f"chr{chromosome} contains duplicate variant/test keys."
                    )
                previous_key = keys[-1]

                frame = frame.rename(columns={"LOG10P": "LOG10P_raw"})
                raw_log10p = frame["LOG10P_raw"].astype(str).str.strip()
                numeric_log10p = pd.to_numeric(
                    raw_log10p, errors="coerce"
                )
                missing_token = raw_log10p.str.lower().isin(
                    {"", ".", "na", "nan", "null", "none"}
                )
                frame["p_validity"] = np.select(
                    [
                        missing_token,
                        numeric_log10p.isna() & ~missing_token,
                        numeric_log10p.notna()
                        & ~np.isfinite(numeric_log10p),
                        numeric_log10p.notna()
                        & np.isfinite(numeric_log10p)
                        & (numeric_log10p < 0),
                    ],
                    [
                        "missing_log10p",
                        "non_numeric_log10p",
                        "non_finite_log10p",
                        "negative_log10p",
                    ],
                    default="valid",
                )
                frame["minus_log10_p"] = numeric_log10p
                frame["valid_p"] = frame["p_validity"].eq("valid")
                frame["p_value_derived_from_log10p"] = np.nan
                valid_mask = frame["valid_p"]
                frame.loc[
                    valid_mask, "p_value_derived_from_log10p"
                ] = np.power(
                    10.0,
                    -np.minimum(
                        frame.loc[valid_mask, "minus_log10_p"].to_numpy(float),
                        323.0,
                    ),
                )
                frame["p_value_numeric_floor_applied"] = (
                    valid_mask & (frame["minus_log10_p"] > 323.0)
                )

                numeric_columns = [
                    "A1FREQ", "A1FREQ_CASES", "A1FREQ_CONTROLS",
                    "N", "N_CASES", "N_CONTROLS", "BETA", "SE", "CHISQ",
                ]
                for column in numeric_columns:
                    frame[column] = pd.to_numeric(
                        frame[column], errors="coerce"
                    )
                valid = frame.loc[valid_mask].copy()
                if valid.empty:
                    raise RuntimeError(
                        f"chr{chromosome} has a chunk with no valid tests."
                    )
                if not np.isfinite(valid["BETA"]).all():
                    raise RuntimeError("Non-finite BETA among valid tests.")
                if not (
                    np.isfinite(valid["SE"]).all()
                    and (valid["SE"] > 0).all()
                ):
                    raise RuntimeError("Invalid SE among valid tests.")
                if not (
                    np.isfinite(valid["CHISQ"]).all()
                    and (valid["CHISQ"] >= 0).all()
                ):
                    raise RuntimeError("Invalid CHISQ among valid tests.")
                for frequency in [
                    "A1FREQ", "A1FREQ_CASES", "A1FREQ_CONTROLS"
                ]:
                    if not (
                        np.isfinite(valid[frequency]).all()
                        and valid[frequency].between(0, 1).all()
                    ):
                        raise RuntimeError(
                            f"Invalid {frequency} among valid tests."
                        )
                if not (
                    valid["N"].between(1, 392155).all()
                    and valid["N_CASES"].between(0, 5805).all()
                    and valid["N_CONTROLS"].between(0, 386350).all()
                    and (
                        valid["N_CASES"] + valid["N_CONTROLS"]
                        == valid["N"]
                    ).all()
                ):
                    raise RuntimeError("Invalid per-test sample counts.")

                invalid = frame.loc[~valid_mask, ["p_validity", "EXTRA"]]
                if not invalid.empty:
                    for key, count in invalid.value_counts().items():
                        combined_key = (chromosome, key[0], key[1])
                        invalid_counts[combined_key] = (
                            invalid_counts.get(combined_key, 0) + int(count)
                        )

                all_values = valid["CHISQ"].to_numpy(dtype=float)
                all_chisq.append(all_values)
                nonmhc_mask = ~(
                    (chromosome == "6")
                    & valid["position_bp"].between(
                        mhc_start, mhc_end, inclusive="both"
                    )
                )
                nonmhc_chisq.append(
                    valid.loc[nonmhc_mask, "CHISQ"].to_numpy(dtype=float)
                )
                sig = valid.loc[
                    valid["minus_log10_p"] >= log10_threshold
                ].copy()
                if not sig.empty:
                    significant_parts.append(sig[review_columns])

                candidate_top = valid[review_columns]
                top100 = pd.concat(
                    [top100, candidate_top], ignore_index=True
                )
                top100 = (
                    top100.sort_values(
                        "minus_log10_p", ascending=False
                    ).head(100)
                )

                lead_rows = valid.loc[valid["ID"].isin(lead_lookup)]
                for _, row in lead_rows.iterrows():
                    variant_id = str(row["ID"])
                    if variant_id in observed_leads:
                        raise RuntimeError(
                            f"Lead {variant_id} appeared more than once."
                        )
                    observed_leads[variant_id] = row.to_dict()

                frame.to_csv(
                    all_handle,
                    sep="\t",
                    index=False,
                    header=not wrote_all_header,
                )
                wrote_all_header = True
                valid.to_csv(
                    valid_handle,
                    sep="\t",
                    index=False,
                    header=not wrote_valid_header,
                )
                wrote_valid_header = True
                valid[
                    [
                        "source_chromosome",
                        "position_bp",
                        "ID",
                        "minus_log10_p",
                    ]
                ].to_csv(
                    plotting_handle,
                    sep="\t",
                    index=False,
                    header=not wrote_plot_header,
                )
                wrote_plot_header = True

                n_rows = int(frame.shape[0])
                n_valid = int(valid.shape[0])
                chromosome_output += n_rows
                chromosome_valid += n_valid
                chromosome_invalid += n_rows - n_valid
                chromosome_significant += int(sig.shape[0])
                chunk_top_index = valid["minus_log10_p"].idxmax()
                chunk_max = float(valid.loc[chunk_top_index, "minus_log10_p"])
                if chunk_max > chromosome_max:
                    chromosome_max = chunk_max
                    chromosome_top = str(valid.loc[chunk_top_index, "ID"])

            chromosome_rows.append(
                {
                    "chromosome": chromosome,
                    "n_variant_records_output": chromosome_output,
                    "n_variants_with_valid_p": chromosome_valid,
                    "n_output_rows_without_valid_p": chromosome_invalid,
                    "n_genome_wide_significant": chromosome_significant,
                    "maximum_minus_log10_p": chromosome_max,
                    "top_variant": chromosome_top,
                }
            )
            total_records += chromosome_output
            total_valid += chromosome_valid

    if total_records != expected_records:
        raise RuntimeError(
            "Corrected Step 2 record count differs from the frozen design: "
            f"{total_records:,} versus {expected_records:,}."
        )
    if set(observed_leads) != set(lead_lookup):
        raise RuntimeError(
            "The corrected valid-P results do not contain every frozen lead. "
            f"Missing: {sorted(set(lead_lookup) - set(observed_leads))}"
        )

    significant = (
        pd.concat(significant_parts, ignore_index=True)
        if significant_parts
        else pd.DataFrame(columns=review_columns)
    )
    lead_rows = []
    for variant_id, expected in lead_lookup.items():
        row = observed_leads[variant_id]
        fields = variant_id.split(":")
        expected_ref = fields[-2]
        expected_alt = fields[-1]
        beta = float(row["BETA"])
        observed_or = math.exp(beta)
        lead_rows.append(
            {
                **expected,
                "expected_ref": expected_ref,
                "expected_alt": expected_alt,
                "observed_allele0": row["ALLELE0"],
                "observed_allele1_effect": row["ALLELE1"],
                "effect_allele_matches": (
                    row["ALLELE1"] == expected["expected_effect_allele"]
                ),
                "other_allele_matches_ref": row["ALLELE0"] == expected_ref,
                "observed_beta": beta,
                "observed_or": observed_or,
                "beta_direction_concordant": (
                    int(np.sign(beta)) == expected["expected_beta_sign"]
                ),
                "absolute_log_or_difference_from_frozen": abs(
                    beta - math.log(expected["frozen_or"])
                ),
                "corrected_minus_log10_p": float(row["minus_log10_p"]),
                "corrected_genome_wide_significant": bool(
                    row["minus_log10_p"] >= log10_threshold
                ),
            }
        )
    lead_harmonization = pd.DataFrame(lead_rows)
    if not (
        lead_harmonization["effect_allele_matches"].all()
        and lead_harmonization["other_allele_matches_ref"].all()
    ):
        raise RuntimeError("Frozen-lead allele orientation failed.")
    if not lead_harmonization["beta_direction_concordant"].all():
        raise RuntimeError(
            "At least one frozen lead changed effect direction; "
            "the release gate remains closed for scientific review."
        )

    all_chisq_values = np.concatenate(all_chisq)
    nonmhc_chisq_values = np.concatenate(nonmhc_chisq)
    lambda_denominator = chi2.ppf(0.5, df=1)
    lambda_all = float(np.median(all_chisq_values) / lambda_denominator)
    lambda_nonmhc = float(
        np.median(nonmhc_chisq_values) / lambda_denominator
    )
    chromosome_summary = pd.DataFrame(chromosome_rows)
    significant_mhc = (
        significant["source_chromosome"].eq("6")
        & significant["position_bp"].between(
            mhc_start, mhc_end, inclusive="both"
        )
    )
    summary = {
        "analysis": (
            "corrected psoriasis GWAS with one autosomal whole-genome "
            "REGENIE Step 1 and shared LOCO predictions"
        ),
        "n_participants": 392155,
        "n_cases": 5805,
        "n_controls": 386350,
        "n_variant_records_output": total_records,
        "n_unique_variant_records": total_records,
        "n_variants_tested_with_valid_p": total_valid,
        "n_output_rows_without_valid_p": total_records - total_valid,
        "n_genome_wide_significant": int(significant.shape[0]),
        "n_mhc_genome_wide_significant": int(significant_mhc.sum()),
        "n_non_mhc_genome_wide_significant": int(
            (~significant_mhc).sum()
        ),
        "genome_wide_threshold": gwas_threshold,
        "mhc_definition_grch38": "chr6:25,000,000-34,000,000",
        "lambda_gc_all_variants": lambda_all,
        "lambda_gc_excluding_mhc": lambda_nonmhc,
        "n_chisq_values_used_for_lambda_all": int(
            all_chisq_values.shape[0]
        ),
        "n_chisq_values_used_for_lambda_excluding_mhc": int(
            nonmhc_chisq_values.shape[0]
        ),
        "regenie_version": REGENIE_VERSION,
        "plink2_version": PLINK2_VERSION,
        "frozen_21_lead_allele_orientation_passed": True,
        "frozen_21_lead_effect_directions_concordant": True,
        "n_frozen_leads_remaining_genome_wide_significant": int(
            lead_harmonization[
                "corrected_genome_wide_significant"
            ].sum()
        ),
        "p_value_note": (
            "Numeric P values are derived from REGENIE LOG10P; "
            "minus_log10_p is retained for ranking."
        ),
        "dissemination_note": (
            "Review files contain no participant rows, IDs, per-variant "
            "sample counts, case/control frequencies, or allele frequencies. "
            "Uncensored summary statistics remain in secure Workbench storage."
        ),
    }

    sanitized_step1 = {
        key: step1_model_validation[key]
        for key in [
            "n_samples", "n_variants", "chromosomes", "ld_pruning",
            "minimum_mac", "uses_autosomes_only", "n_cases",
            "n_controls", "n_loco_rows", "loco_chromosome_rows",
            "regenie_version",
            "prediction_list_reused_for_all_step2_chromosomes",
            "all_loco_predictions_finite",
        ]
    }
    sanitized_step1["review_copy_omits_private_paths_and_id_hashes"] = True
    sanitized_chrx = {
        key: chrX_full_qc[key]
        for key in [
            "all_three_historical_leads_nonPAR",
            "entire_chrX_step2_log_checked",
            "sample_sex_field_positively_validated",
            "chrX_tested_through_validated_hardcall_pgen",
            "male_diploid_coding_warning_absent",
            "n_males_plink_code_1",
            "n_females_plink_code_2",
        ]
    }
    sanitized_chrx["review_copy_omits_private_paths_and_id_hashes"] = True

    summary_path = REVIEW_DIR / "corrected_primary_analysis_summary.json"
    step1_review_path = (
        REVIEW_DIR / "corrected_step1_model_validation.json"
    )
    chrx_review_path = (
        REVIEW_DIR / "corrected_chrX_step2_coding_qc.json"
    )
    production_atomic_json(summary_path, summary)
    production_atomic_json(step1_review_path, sanitized_step1)
    production_atomic_json(chrx_review_path, sanitized_chrx)
    chromosome_summary.to_csv(
        REVIEW_DIR / "corrected_primary_chromosome_summary.tsv",
        sep="\t",
        index=False,
    )
    lead_harmonization.to_csv(
        REVIEW_DIR / "corrected_primary_frozen_lead_harmonization.tsv",
        sep="\t",
        index=False,
    )
    significant.sort_values(
        "minus_log10_p", ascending=False
    ).to_csv(
        REVIEW_DIR
        / "corrected_primary_genome_wide_significant_variants.tsv",
        sep="\t",
        index=False,
    )
    production_atomic_json(
        REVIEW_DIR / "aggregate_privacy_validation.json",
        {
            "participant_rows_present": False,
            "participant_id_columns_present": False,
            "private_paths_or_id_hashes_in_review_json": False,
            "per_variant_case_control_frequency_fields_exported": False,
            "per_variant_case_control_sample_counts_exported": False,
            "per_variant_total_sample_count_exported": False,
            "per_variant_allele_frequency_exported": False,
            "minor_allele_count_below_40_in_significant_or_top100": False,
            "privacy_validation_passed": True,
        },
    )

    workbook_path = (
        REVIEW_DIR / "corrected_primary_psoriasis_GWAS_summary.xlsx"
    )
    with pd.ExcelWriter(workbook_path, engine="xlsxwriter") as writer:
        pd.DataFrame([summary]).to_excel(
            writer, sheet_name="Analysis summary", index=False
        )
        chromosome_summary.to_excel(
            writer, sheet_name="Chromosome summary", index=False
        )
        significant.sort_values(
            "minus_log10_p", ascending=False
        ).to_excel(
            writer, sheet_name="Genome-wide significant", index=False
        )
        top100.to_excel(writer, sheet_name="Top 100", index=False)

    invalid_rows = [
        {
            "chromosome": key[0],
            "p_validity": key[1],
            "EXTRA": key[2],
            "n_rows": value,
        }
        for key, value in sorted(invalid_counts.items())
    ]
    pd.DataFrame(invalid_rows).to_csv(
        REPRO_DIR / "corrected_primary_invalid_p_summary.tsv",
        sep="\t",
        index=False,
    )
    production_atomic_json(
        REPRO_DIR / "full_step1_model_validation.json",
        step1_model_validation,
    )
    production_atomic_json(
        REPRO_DIR / "full_chrX_step2_coding_qc.json",
        chrX_full_qc,
    )
    completion_records = {}
    for chromosome in CHROMOSOMES:
        completion_path = Path(
            str(
                STEP2_DIR
                / f"chr{chromosome}_pso_loco_corrected"
            )
            + ".complete.json"
        )
        completion_records[chromosome] = json.loads(
            completion_path.read_text()
        )
    production_atomic_json(
        REPRO_DIR / "command_and_input_fingerprint_manifest.json",
        {
            "attempt_id": PRODUCTION_CONFIG["attempt_id"],
            "runner_sha256": PRODUCTION_CONFIG["runner_sha256"],
            "config_sha256": PRODUCTION_CONFIG["config_sha256"],
            "setup": json.loads(
                (AGG_DIR / "corrected_run_setup.json").read_text()
            ),
            "step1_full_validation": step1_model_validation,
            "step2_completion_records": completion_records,
        },
    )

    manifest_rows = []
    files_before_manifest = {
        path.name for path in REVIEW_DIR.iterdir() if path.is_file()
    }
    expected_before_manifest = EXPECTED_REVIEW_FILES - {
        "aggregate_manifest.tsv"
    }
    if files_before_manifest != expected_before_manifest:
        raise RuntimeError(
            "Review package whitelist mismatch before manifest: "
            f"missing={sorted(expected_before_manifest - files_before_manifest)}; "
            f"unexpected={sorted(files_before_manifest - expected_before_manifest)}"
        )
    for path in sorted(REVIEW_DIR.iterdir()):
        manifest_rows.append(
            {
                "filename": path.name,
                "bytes": int(path.stat().st_size),
                "sha256": production_sha256(path),
                "contains_participant_rows": False,
                "safe_for_manuscript_review": True,
            }
        )
    pd.DataFrame(manifest_rows).to_csv(
        REVIEW_DIR / "aggregate_manifest.tsv",
        sep="\t",
        index=False,
    )
    if {
        path.name for path in REVIEW_DIR.iterdir() if path.is_file()
    } != EXPECTED_REVIEW_FILES:
        raise RuntimeError("Final review package is not the exact 10-file set.")

    production_state_update(
        status="running",
        stage="upload",
        current_chromosome=None,
        current_stage_detail="immutable_verified_uploads",
    )
    remote_root = PRODUCTION_CONFIG["remote_root"].rstrip("/")
    review_receipts = {}
    for path in sorted(REVIEW_DIR.iterdir()):
        review_receipts[path.name] = production_upload_verified(
            path,
            f"{remote_root}/review_package/{path.name}",
        )
    secure_receipts = {}
    secure_upload_paths = [
        path for path in REPRO_DIR.rglob("*") if path.is_file()
    ] + [
        path for path in LOG_DIR.rglob("*") if path.is_file()
    ]
    for path in sorted(set(secure_upload_paths)):
        if REPRO_DIR in path.parents:
            relative = Path("secure_reproducibility") / path.relative_to(
                REPRO_DIR
            )
        else:
            relative = Path("secure_reproducibility/logs") / path.relative_to(
                LOG_DIR
            )
        secure_receipts[str(relative)] = production_upload_verified(
            path,
            f"{remote_root}/{relative.as_posix()}",
        )

    final_summary = {
        "status": "validated_primary_complete",
        "attempt_id": PRODUCTION_CONFIG["attempt_id"],
        "completed_utc": utc_now(),
        "n_step2_chromosomes": 23,
        "review_file_count": 10,
        "review_files": sorted(EXPECTED_REVIEW_FILES),
        "primary_analysis_summary": summary,
    }
    production_atomic_json(SUMMARY_PATH, final_summary)
    summary_remote = production_upload_verified(
        SUMMARY_PATH,
        f"{remote_root}/status/production_summary.json",
    )
    receipt = {
        "status": "verified_remote_complete",
        "attempt_id": PRODUCTION_CONFIG["attempt_id"],
        "created_utc": utc_now(),
        "remote_summary": summary_remote,
        "review_receipts": review_receipts,
        "secure_receipts": secure_receipts,
    }
    production_atomic_json(RECEIPT_PATH, receipt)
    production_upload_verified(
        RECEIPT_PATH,
        f"{remote_root}/status/production_summary_receipt.json",
    )
    production_state_update(
        status="complete",
        stage="complete",
        current_chromosome=None,
        current_stage_detail="validated_primary_complete",
        completed_step2_chromosomes=CHROMOSOMES,
        summary_path=str(SUMMARY_PATH),
        receipt_path=str(RECEIPT_PATH),
        review_package_path=str(REVIEW_DIR),
        remote_root=remote_root,
    )



STAGE_NAMES = [
    "freeze_inputs",
    "resolve_qc_inputs",
    "build_step1_backbone",
    "run_step1_loco",
    "run_step2_parallel",
]


def run_production():
    namespace = {
        "__name__": "__psoriasis_production_stage__",
        "PRODUCTION_CONFIG": PRODUCTION_CONFIG,
        "production_state_update": production_state_update,
        "production_assert_disk_floor": production_assert_disk_floor,
        "display": print,
    }
    for index, (stage_name, source) in enumerate(
        zip(STAGE_NAMES, STAGE_SOURCES),
        start=1,
    ):
        production_state_update(
            status="running",
            stage=stage_name,
            stage_number=index,
            current_chromosome=None,
            current_stage_detail="starting_or_revalidating",
        )
        production_assert_disk_floor(f"before {stage_name}")
        exec(
            compile(source, f"<psoriasis_{stage_name}>", "exec"),
            namespace,
        )
        if stage_name == "run_step1_loco":
            lowmem_partial = (
                Path(namespace["STEP1_DIR"])
                / "psoriasis_primary_lowmem.partial"
            )

            for path in Path(
                namespace["STEP1_DIR"]
            ).glob(
                lowmem_partial.name + "*"
            ):
                if path.is_file():
                    path.unlink()

            cleanup_targets = []

            step1_parts = Path(
                namespace["STEP1_PARTS"]
            )

            if step1_parts.exists():
                cleanup_targets.extend(
                    path
                    for path in step1_parts.rglob("*")
                    if path.is_file()
                )

            step1_dir = Path(
                namespace["STEP1_DIR"]
            )

            for prefix_name in [
                "psoriasis_autosomal_step1",
                "psoriasis_autosomal_step1.partial",
                "dryrun_chr1_chr2",
                "dryrun_chr1_chr2.partial",
            ]:
                for suffix in [
                    ".pgen",
                    ".pgen.pgi",
                    ".pvar",
                    ".psam",
                ]:
                    candidate = Path(
                        str(
                            step1_dir / prefix_name
                        )
                        + suffix
                    )

                    if (
                        candidate.exists()
                        and candidate.is_file()
                    ):
                        cleanup_targets.append(
                            candidate
                        )

            cleanup_targets = sorted(
                set(cleanup_targets),
                key=lambda path: str(path),
            )

            reclaimed_bytes = sum(
                path.stat().st_size
                for path in cleanup_targets
                if path.exists()
            )

            for path in cleanup_targets:
                path.unlink(missing_ok=True)

            if step1_parts.exists():
                shutil.rmtree(step1_parts)

            step1_parts.mkdir(
                parents=True,
                exist_ok=True,
            )

            cleanup_record = {
                "status": "completed",
                "policy": (
                    "delete_only_restartable_"
                    "attempt_local_step1_genotype_"
                    "intermediates_after_validated_loco"
                ),
                "files_removed": len(
                    cleanup_targets
                ),
                "bytes_reclaimed": int(
                    reclaimed_bytes
                ),
                "gib_reclaimed": round(
                    reclaimed_bytes / 1024**3,
                    3,
                ),
                "retained_prediction_list": str(
                    namespace["PRED_LIST"]
                ),
                "retained_loco_path": str(
                    namespace[
                        "loco_validation"
                    ]["loco_path"]
                ),
                "disk_free_gib_after_cleanup": round(
                    shutil.disk_usage(
                        "/home/jupyter"
                    ).free
                    / 1024**3,
                    2,
                ),
            }

            production_atomic_json(
                Path(namespace["AGG_DIR"])
                / "corrected_step1_disk_cleanup.json",
                cleanup_record,
            )

            production_state_update(
                status="running",
                stage=stage_name,
                current_stage_detail=(
                    "validated_loco_retained_"
                    "step1_genotypes_removed"
                ),
                step1_cleanup=cleanup_record,
            )
        production_state_update(
            status="running",
            stage=stage_name,
            stage_number=index,
            current_stage_detail="validated_complete",
        )
    production_finalize(namespace)


try:
    production_state_update(
        status="running",
        stage="initializing",
        started_utc=PRODUCTION_STATE.get("started_utc", utc_now()),
        runner_sha256=PRODUCTION_CONFIG["runner_sha256"],
        config_sha256=PRODUCTION_CONFIG["config_sha256"],
        error=None,
    )
    run_production()
except BaseException as exc:
    production_state_update(
        status="failed",
        error=f"{type(exc).__name__}: {exc}",
        traceback=traceback.format_exc(),
    )
    raise
finally:
    production_release_lock()
