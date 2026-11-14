"""Resume-safe disk-bounded assembly of autosomal REGENIE Step 1 BED."""

import base64
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


if len(sys.argv) != 2:
    raise SystemExit("Use: assembly_runner.py PLAN.json")

PLAN_PATH = Path(sys.argv[1])
PLAN = json.loads(PLAN_PATH.read_text())
PREFIX = Path(PLAN["assembled_prefix"])
BED = Path(str(PREFIX) + ".bed")
BIM = Path(str(PREFIX) + ".bim")
FAM = Path(str(PREFIX) + ".fam")
STATE_PATH = Path(PLAN["state_path"])
PID_PATH = Path(PLAN["pid_path"])
TEMPORARY_ROOT = Path(PLAN["temporary_root"])
LOCK_PATH = STATE_PATH.parent / "step1_assembly.lock"
RECEIPT_PATH = (
    STATE_PATH.parent / "step1_assembly_validation.json"
)
BED_MAGIC = b"\x6c\x1b\x01"


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2))
    temporary.replace(path)


def md5_base64(path):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(
            lambda: handle.read(32 * 1024**2),
            b"",
        ):
            digest.update(block)
    return base64.b64encode(digest.digest()).decode()


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(
            lambda: handle.read(32 * 1024**2),
            b"",
        ):
            digest.update(block)
    return digest.hexdigest()


def current_free_bytes():
    return shutil.disk_usage("/home/jupyter").free


def update_state(state, **updates):
    state.update(updates)
    state["updated_utc"] = now()
    state["free_disk_gib"] = round(
        current_free_bytes() / 1024**3,
        2,
    )
    write_json(STATE_PATH, state)


def run_checked(command, label):
    print(f"\n[{label}]", flush=True)
    print(" ".join(str(value) for value in command), flush=True)
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.stdout:
        print(completed.stdout, flush=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{label} failed with exit code "
            f"{completed.returncode}."
        )


def verify_artifact(path, artifact, label):
    path = Path(path)
    expected_size = int(artifact["size_bytes"])
    if not path.exists() or path.stat().st_size != expected_size:
        raise RuntimeError(
            f"{label} has the wrong local byte size."
        )
    observed_md5 = md5_base64(path)
    if observed_md5 != artifact["md5_base64"]:
        raise RuntimeError(
            f"{label} failed its frozen MD5 verification."
        )


def download_artifact(artifact, destination, label):
    destination = Path(destination)
    partial = destination.with_suffix(
        destination.suffix + ".partial"
    )
    destination.unlink(missing_ok=True)
    partial.unlink(missing_ok=True)
    run_checked(
        [
            "gsutil",
            "-u",
            PLAN["project"],
            "-o",
            "GSUtil:parallel_composite_upload_threshold=0",
            "cp",
            artifact["versioned_uri"],
            str(partial),
        ],
        label,
    )
    verify_artifact(partial, artifact, label)
    partial.replace(destination)


def count_nonempty_lines(path):
    count = 0
    with Path(path).open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        for line in handle:
            if line.strip():
                count += 1
    return count


def truncate_to_checkpoint(state):
    completed = state.get("completed_chromosomes", [])
    expected_bed = int(state["assembled_bed_bytes"])
    expected_bim = int(state["assembled_bim_bytes"])

    if not completed:
        PREFIX.parent.mkdir(parents=True, exist_ok=True)
        with BED.open("wb") as handle:
            handle.write(BED_MAGIC)
            handle.flush()
            os.fsync(handle.fileno())
        with BIM.open("wb"):
            pass
        state["assembled_bed_bytes"] = 3
        state["assembled_bim_bytes"] = 0
        state["assembled_variants"] = 0
        write_json(STATE_PATH, state)
        return

    for path, expected, label in [
        (BED, expected_bed, "assembled BED"),
        (BIM, expected_bim, "assembled BIM"),
    ]:
        if not path.exists():
            raise RuntimeError(
                f"Checkpointed {label} is missing."
            )
        observed = path.stat().st_size
        if observed < expected:
            raise RuntimeError(
                f"Checkpointed {label} is shorter than its "
                "validated checkpoint."
            )
        if observed > expected:
            with path.open("r+b") as handle:
                handle.truncate(expected)
                handle.flush()
                os.fsync(handle.fileno())

    with BED.open("rb") as handle:
        if handle.read(3) != BED_MAGIC:
            raise RuntimeError(
                "Checkpointed assembled BED has invalid magic bytes."
            )


def prepare_fam():
    first = PLAN["artifact_plan"][0]
    artifact = first["fam"]
    if FAM.exists():
        verify_artifact(FAM, artifact, "assembled FAM")
    else:
        download_artifact(
            artifact,
            FAM,
            "Download and verify canonical autosomal FAM",
        )
    observed_rows = count_nonempty_lines(FAM)
    if observed_rows != int(
        PLAN["expected_analysis_samples"]
    ):
        raise RuntimeError(
            "Canonical FAM has the wrong sample count."
        )


def append_chromosome(row, state):
    chromosome = str(row["chromosome"])
    required_free = (
        int(row["bed"]["size_bytes"])
        + int(
            float(PLAN["minimum_protected_free_gib"])
            * 1024**3
        )
    )
    if current_free_bytes() < required_free:
        raise RuntimeError(
            f"Insufficient protected space before chromosome "
            f"{chromosome}."
        )

    TEMPORARY_ROOT.mkdir(parents=True, exist_ok=True)
    temporary_bed = (
        TEMPORARY_ROOT / f"chr{chromosome}_step1.bed"
    )
    temporary_bim = (
        TEMPORARY_ROOT / f"chr{chromosome}_step1.bim"
    )

    update_state(
        state,
        status="running",
        current_chromosome=chromosome,
        current_stage="download_bed",
        error=None,
    )
    download_artifact(
        row["bed"],
        temporary_bed,
        f"Download and verify chr{chromosome} Step 1 BED",
    )

    update_state(
        state,
        current_stage="download_bim",
    )
    download_artifact(
        row["bim"],
        temporary_bim,
        f"Download and verify chr{chromosome} Step 1 BIM",
    )

    expected_variants = int(row["step1_variant_count"])
    observed_bim_rows = count_nonempty_lines(temporary_bim)
    if observed_bim_rows != expected_variants:
        raise RuntimeError(
            f"chr{chromosome} BIM has "
            f"{observed_bim_rows:,} rows; expected "
            f"{expected_variants:,}."
        )

    with temporary_bed.open("rb") as source:
        if source.read(3) != BED_MAGIC:
            raise RuntimeError(
                f"chr{chromosome} BED is not variant-major BED."
            )
        update_state(
            state,
            current_stage="append_bed",
        )
        with BED.open("ab") as destination:
            shutil.copyfileobj(
                source,
                destination,
                32 * 1024**2,
            )
            destination.flush()
            os.fsync(destination.fileno())

    update_state(
        state,
        current_stage="append_bim",
    )
    with temporary_bim.open("rb") as source:
        content = source.read()
    if content and not content.endswith(b"\n"):
        content += b"\n"
    with BIM.open("ab") as destination:
        destination.write(content)
        destination.flush()
        os.fsync(destination.fileno())

    new_bed_bytes = BED.stat().st_size
    new_bim_bytes = BIM.stat().st_size
    expected_new_bed_bytes = (
        int(state["assembled_bed_bytes"])
        + int(row["bed"]["size_bytes"])
        - 3
    )
    if new_bed_bytes != expected_new_bed_bytes:
        raise RuntimeError(
            f"chr{chromosome} BED append size mismatch."
        )

    temporary_bed.unlink()
    temporary_bim.unlink()

    completed = list(state["completed_chromosomes"])
    completed.append(chromosome)
    update_state(
        state,
        completed_chromosomes=completed,
        assembled_bed_bytes=new_bed_bytes,
        assembled_bim_bytes=new_bim_bytes,
        assembled_variants=(
            int(state["assembled_variants"])
            + expected_variants
        ),
        current_chromosome=None,
        current_stage="chromosome_checkpointed",
    )


def validate_final(state):
    update_state(
        state,
        current_stage="validate_final_sizes",
    )
    if BED.stat().st_size != int(
        PLAN["expected_merged_bed_bytes"]
    ):
        raise RuntimeError(
            "Final assembled BED has the wrong byte size."
        )
    if int(state["assembled_variants"]) != int(
        PLAN["expected_step1_variants"]
    ):
        raise RuntimeError(
            "Final assembled variant count is incorrect."
        )
    if count_nonempty_lines(BIM) != int(
        PLAN["expected_step1_variants"]
    ):
        raise RuntimeError(
            "Final assembled BIM row count is incorrect."
        )
    if count_nonempty_lines(FAM) != int(
        PLAN["expected_analysis_samples"]
    ):
        raise RuntimeError(
            "Final assembled FAM row count is incorrect."
        )

    integrity_prefix = (
        PREFIX.parent / "psoriasis_autosomal_step1_integrity"
    )
    update_state(
        state,
        current_stage="plink_bed_integrity_validation",
    )
    run_checked(
        [
            "plink2",
            "--bfile",
            str(PREFIX),
            "--validate",
            "--write-snplist",
            "--threads",
            "24",
            "--out",
            str(integrity_prefix),
        ],
        "Validate assembled autosomal Step 1 BED",
    )
    snplist = Path(str(integrity_prefix) + ".snplist")
    if (
        not snplist.exists()
        or count_nonempty_lines(snplist)
        != int(PLAN["expected_step1_variants"])
    ):
        raise RuntimeError(
            "Assembled BED integrity variant count is incorrect."
        )

    receipt = {
        "status": "validated",
        "version": PLAN["version"],
        "production_attempt_id": (
            PLAN["production_attempt_id"]
        ),
        "plan_sha256": PLAN["plan_sha256"],
        "chromosomes": list(PLAN["autosomes"]),
        "analysis_samples": int(
            PLAN["expected_analysis_samples"]
        ),
        "step1_variants": int(
            PLAN["expected_step1_variants"]
        ),
        "bed_bytes": BED.stat().st_size,
        "bed_sha256": sha256_file(BED),
        "bim_sha256": sha256_file(BIM),
        "fam_sha256": sha256_file(FAM),
        "plink_validation": True,
        "whole_bgen_downloaded": False,
        "minimum_protected_free_gib": float(
            PLAN["minimum_protected_free_gib"]
        ),
        "final_free_gib": round(
            current_free_bytes() / 1024**3,
            2,
        ),
        "completed_utc": now(),
    }
    write_json(RECEIPT_PATH, receipt)
    update_state(
        state,
        status="complete",
        current_chromosome=None,
        current_stage="complete",
        error=None,
    )


def main():
    state = read_json(STATE_PATH)
    if state.get("status") == "complete":
        print("Step 1 assembly is already complete.")
        return
    if LOCK_PATH.exists():
        raise RuntimeError(
            "A Step 1 assembly lock already exists."
        )
    LOCK_PATH.mkdir()
    PID_PATH.write_text(str(os.getpid()))
    try:
        TEMPORARY_ROOT.mkdir(parents=True, exist_ok=True)
        for child in TEMPORARY_ROOT.iterdir():
            if child.is_file():
                child.unlink()
        truncate_to_checkpoint(state)
        prepare_fam()
        completed = set(state["completed_chromosomes"])
        for row in PLAN["artifact_plan"]:
            chromosome = str(row["chromosome"])
            if chromosome in completed:
                continue
            append_chromosome(row, state)
            completed.add(chromosome)
        validate_final(state)
        print(
            "\nDISK-BOUNDED AUTOSOMAL STEP 1 ASSEMBLY COMPLETE",
            flush=True,
        )
    except Exception as error:
        state = read_json(STATE_PATH)
        update_state(
            state,
            status="failed",
            error=f"{type(error).__name__}: {error}",
        )
        print(f"ERROR: {error}", flush=True)
        raise
    finally:
        PID_PATH.unlink(missing_ok=True)
        if LOCK_PATH.exists():
            LOCK_PATH.rmdir()


if __name__ == "__main__":
    main()
