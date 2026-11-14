import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import zstandard as zstd


if len(sys.argv) != 2:
    raise SystemExit("Usage: runner.py STEP2_PLAN.json")

PLAN_PATH = Path(sys.argv[1])
plan = json.loads(PLAN_PATH.read_text())

ROOT = Path(plan["local_root"])
STATUS_ROOT = ROOT / "status"
RESULT_ROOT = Path(plan["local_result_root"])
SCRATCH_ROOT = Path(plan["local_scratch_root"])
CHUNK_STATUS_ROOT = STATUS_ROOT / "chunks"

STATE_PATH = STATUS_ROOT / "step2_state.json"
SUMMARY_PATH = STATUS_ROOT / "step2_summary.json"
PILOT_RECEIPT = STATUS_ROOT / "step2_pilot_validation.json"

BILLING_PROJECT = f"{os.environ['GOOGLE_PROJECT']}"
MAX_WORKERS = 4
THREADS_PER_JOB = 12
BLOCK_SIZE = 400
MAX_RETRIES = 3
MINIMUM_FREE_GIB = 100

state_lock = threading.Lock()

for directory in [
    STATUS_ROOT,
    RESULT_ROOT,
    SCRATCH_ROOT,
    CHUNK_STATUS_ROOT,
]:
    directory.mkdir(parents=True, exist_ok=True)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, payload):
    path = Path(path)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(payload, indent=2))
    os.replace(partial, path)


def sha256_file(path):
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for block in iter(
            lambda: handle.read(8 * 1024**2),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def sha256_payload(payload):
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def run_command(command):
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(command)
            + "\n\n"
            + (result.stdout or "")[-8000:]
        )

    return result


def gs_copy(source, destination):
    run_command([
        "gsutil",
        "-u", BILLING_PROJECT,
        "-o",
        "GSUtil:parallel_composite_upload_threshold=0",
        "cp",
        str(source),
        str(destination),
    ])


def decompress_pvar(source_zst, destination):
    decompressor = zstd.ZstdDecompressor()

    with Path(source_zst).open("rb") as source:
        with Path(destination).open("wb") as target:
            decompressor.copy_stream(source, target)

    if (
        not Path(destination).exists()
        or Path(destination).stat().st_size == 0
    ):
        raise RuntimeError("PVAR decompression failed.")


def task_id(task):
    return (
        f"chr{task['chromosome']}:{task['stem']}"
    )


def receipt_path(task):
    return (
        CHUNK_STATUS_ROOT
        / f"chr{task['chromosome']}__{task['stem']}.json"
    )


def validated_receipt(task):
    path = receipt_path(task)

    if not path.exists():
        return None

    try:
        receipt = json.loads(path.read_text())
    except Exception:
        return None

    if (
        receipt.get("status") == "validated"
        and receipt.get("remote_uploaded") is True
    ):
        return receipt

    return None


def write_state():
    state["updated_utc"] = utc_now()
    atomic_json(STATE_PATH, state)


def set_active(identifier, active):
    with state_lock:
        current = set(state["active_chunks"])

        if active:
            current.add(identifier)
        else:
            current.discard(identifier)

        state["active_chunks"] = sorted(current)
        write_state()


def failed_test_counts(log_text):
    match = re.search(
        r"Number of failed tests\s*:\s*\((\d+)/(\d+)\)",
        log_text,
    )

    if not match:
        return 0, 0

    return int(match.group(1)), int(match.group(2))


def process_chunk_once(task, attempt_number):
    chromosome = task["chromosome"]
    stem = task["stem"]
    identifier = task_id(task)

    free_gib = (
        shutil.disk_usage("/home/jupyter").free
        / 1024**3
    )

    if free_gib < MINIMUM_FREE_GIB:
        raise RuntimeError(
            f"Only {free_gib:.1f} GiB remain."
        )

    scratch = (
        SCRATCH_ROOT / f"chr{chromosome}" / stem
    )
    result_directory = (
        RESULT_ROOT / f"chr{chromosome}"
    )

    shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True, exist_ok=True)
    result_directory.mkdir(parents=True, exist_ok=True)

    local_prefix = scratch / stem
    local_pgen = Path(str(local_prefix) + ".pgen")
    local_pvar_zst = Path(
        str(local_prefix) + ".pvar.zst"
    )
    local_pvar = Path(str(local_prefix) + ".pvar")
    local_psam = Path(str(local_prefix) + ".psam")

    output_prefix = result_directory / stem
    console_log = Path(
        str(output_prefix) + ".console.log"
    )
    regenie_log = Path(str(output_prefix) + ".log")
    result_file = Path(
        str(output_prefix)
        + "_psoriasis_case.regenie"
    )

    try:
        gs_copy(task["pgen_uri"], local_pgen)
        gs_copy(task["pvar_uri"], local_pvar_zst)
        decompress_pvar(local_pvar_zst, local_pvar)
        gs_copy(task["psam_uri"], local_psam)

        command = [
            "regenie",
            "--step", "2",
            "--bt",
            "--pgen", str(local_prefix),
            "--phenoFile", plan["phenotype_file"],
            "--phenoCol", plan["phenotype"],
            "--covarFile", plan["covariate_file"],
            "--bsize", str(BLOCK_SIZE),
            "--firth",
            "--approx",
            "--pThresh", "0.01",
            "--pred", plan["prediction_list"],
            "--af-cc",
            "--threads", str(THREADS_PER_JOB),
        ]

        if chromosome == "X":
            command.extend([
                "--par-region", "hg38",
            ])

        command.extend([
            "--out", str(output_prefix),
        ])

        with console_log.open("w") as handle:
            completed = subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
            )

        if completed.returncode != 0:
            recent = (
                console_log.read_text(
                    errors="replace"
                )[-8000:]
            )
            raise RuntimeError(
                f"REGENIE exited with "
                f"{completed.returncode}.\n{recent}"
            )

        if (
            not result_file.exists()
            or result_file.stat().st_size == 0
        ):
            raise RuntimeError(
                "Result file is missing or empty."
            )

        pvar = pd.read_csv(
            local_pvar,
            sep="\t",
            dtype=str,
            usecols=["ID"],
        )

        results = pd.read_csv(
            result_file,
            sep=r"\s+",
            dtype=str,
            low_memory=False,
        )

        required_columns = {
            "CHROM",
            "GENPOS",
            "ID",
            "ALLELE0",
            "ALLELE1",
            "A1FREQ",
            "A1FREQ_CASES",
            "A1FREQ_CONTROLS",
            "N",
            "TEST",
            "BETA",
            "SE",
            "CHISQ",
            "LOG10P",
            "EXTRA",
        }

        missing = (
            required_columns
            - set(results.columns)
        )

        if missing:
            raise RuntimeError(
                f"Missing result columns: {sorted(missing)}"
            )

        if len(results) != len(pvar):
            raise RuntimeError(
                f"{len(results):,} results for "
                f"{len(pvar):,} PVAR variants."
            )

        if set(results["ID"]) != set(pvar["ID"]):
            raise RuntimeError(
                "Result and PVAR variant IDs differ."
            )

        log_text = regenie_log.read_text(
            errors="replace"
        )

        for evidence in [
            "number of individuals used in analysis = 392155",
            "5805 cases and 386350 controls",
            "LOCO predictions",
        ]:
            if evidence not in log_text:
                raise RuntimeError(
                    f"Missing log evidence: {evidence}"
                )

        failed_tests, firth_tests = failed_test_counts(
            log_text
        )

        remote_directory = (
            f"{plan['remote_result_root']}/"
            f"chromosomes/chr{chromosome}/{stem}"
        )

        receipt = {
            "status": "validated",
            "validated_utc": utc_now(),
            "attempt_id": plan["attempt_id"],
            "task_id": identifier,
            "chromosome": chromosome,
            "chunk_stem": stem,
            "attempt_number": attempt_number,
            "source_pgen_uri": task["pgen_uri"],
            "source_pvar_uri": task["pvar_uri"],
            "source_psam_uri": task["psam_uri"],
            "n_source_variants": len(pvar),
            "n_result_rows": len(results),
            "firth_tests": firth_tests,
            "failed_tests": failed_tests,
            "result_file": str(result_file),
            "result_sha256": sha256_file(
                result_file
            ),
            "regenie_log": str(regenie_log),
            "regenie_log_sha256": sha256_file(
                regenie_log
            ),
            "pvar_decompressed_before_regenie": True,
            "par_region": (
                "hg38" if chromosome == "X" else None
            ),
            "sex_specific_filter": None,
            "remote_directory": remote_directory,
            "remote_uploaded": False,
            "participant_rows_in_result": False,
        }

        local_receipt = receipt_path(task)
        atomic_json(local_receipt, receipt)

        for local_file in [
            result_file,
            regenie_log,
            console_log,
        ]:
            gs_copy(
                local_file,
                f"{remote_directory}/{local_file.name}",
            )

        receipt["remote_uploaded"] = True
        atomic_json(local_receipt, receipt)

        gs_copy(
            local_receipt,
            f"{remote_directory}/{local_receipt.name}",
        )

        return receipt

    finally:
        # Deletes only temporary, securely retained inputs.
        shutil.rmtree(scratch, ignore_errors=True)


def process_with_retries(task):
    identifier = task_id(task)
    set_active(identifier, True)

    last_error = None

    try:
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                receipt = process_chunk_once(
                    task,
                    attempt,
                )
                return True, receipt, None

            except Exception as error:
                last_error = (
                    f"{type(error).__name__}: {error}"
                )

                failure = {
                    "status": "retrying",
                    "updated_utc": utc_now(),
                    "task_id": identifier,
                    "chromosome": task["chromosome"],
                    "chunk_stem": task["stem"],
                    "attempt_number": attempt,
                    "error": last_error,
                    "traceback": traceback.format_exc(),
                }

                atomic_json(
                    receipt_path(task),
                    failure,
                )

                if attempt < MAX_RETRIES:
                    time.sleep(10 * attempt)

        return False, None, last_error

    finally:
        set_active(identifier, False)


# Validate the immutable plan.
stored_plan_hash = plan["plan_sha256"]

calculated_plan_hash = sha256_payload({
    key: value
    for key, value in plan.items()
    if key != "plan_sha256"
})

if stored_plan_hash != calculated_plan_hash:
    raise RuntimeError(
        "Step 2 plan checksum failed."
    )

if not PILOT_RECEIPT.exists():
    raise RuntimeError(
        "Corrected pilot receipt is missing."
    )

pilot = json.loads(PILOT_RECEIPT.read_text())

if (
    pilot.get("status") != "validated"
    or not pilot.get("autosomal_pilot_validated")
    or not pilot.get("chromosome_x_pilot_validated")
    or not pilot.get("pvar_decompression_validated")
):
    raise RuntimeError(
        "Corrected Step 2 pilots are not validated."
    )


# Construct the 1,149 chunk tasks.
tasks = []

for chromosome_record in plan["chunk_records"]:
    chromosome = chromosome_record["chromosome"]

    pvar_map = {
        Path(uri).name.removesuffix(".pvar.zst"): uri
        for uri in chromosome_record["pvar_uris"]
    }

    for pgen_uri in chromosome_record["pgen_uris"]:
        stem = (
            Path(pgen_uri).name.removesuffix(".pgen")
        )

        if stem not in pvar_map:
            raise RuntimeError(
                f"No PVAR matches {stem}."
            )

        tasks.append({
            "chromosome": chromosome,
            "stem": stem,
            "pgen_uri": pgen_uri,
            "pvar_uri": pvar_map[stem],
            "psam_uri": chromosome_record["psam_uri"],
        })

if len(tasks) != 1_149:
    raise RuntimeError(
        f"Constructed {len(tasks)} tasks."
    )

completed_ids = {
    task_id(task)
    for task in tasks
    if validated_receipt(task) is not None
}

pending_tasks = [
    task
    for task in tasks
    if task_id(task) not in completed_ids
]

started_utc = utc_now()

if STATE_PATH.exists():
    try:
        old_state = json.loads(STATE_PATH.read_text())
        started_utc = old_state.get(
            "started_utc",
            started_utc,
        )
    except Exception:
        pass

state = {
    "status": "running",
    "attempt_id": plan["attempt_id"],
    "started_utc": started_utc,
    "updated_utc": utc_now(),
    "total_chunks": 1_149,
    "completed_chunks": sorted(completed_ids),
    "failed_chunks": {},
    "active_chunks": [],
    "maximum_concurrent_jobs": MAX_WORKERS,
    "threads_per_job": THREADS_PER_JOB,
}

write_state()

print("=" * 72, flush=True)
print("STARTING/RESUMING CORRECTED STEP 2", flush=True)
print("=" * 72, flush=True)
print("Total chunks: 1,149", flush=True)
print(
    f"Already validated: {len(completed_ids):,}",
    flush=True,
)
print(
    f"Pending: {len(pending_tasks):,}",
    flush=True,
)
print("Concurrent jobs: 4", flush=True)
print("Threads per job: 12", flush=True)

with concurrent.futures.ThreadPoolExecutor(
    max_workers=MAX_WORKERS
) as executor:
    futures = {
        executor.submit(
            process_with_retries,
            task,
        ): task
        for task in pending_tasks
    }

    for future in concurrent.futures.as_completed(
        futures
    ):
        task = futures[future]
        identifier = task_id(task)

        try:
            success, receipt, error = future.result()
        except Exception as unexpected:
            success = False
            error = (
                f"{type(unexpected).__name__}: {unexpected}"
            )

        with state_lock:
            if success:
                completed_ids.add(identifier)
                state["completed_chunks"] = sorted(
                    completed_ids
                )
                state["failed_chunks"].pop(
                    identifier,
                    None,
                )
            else:
                state["failed_chunks"][identifier] = error

            write_state()

        print(
            f"[{len(completed_ids):,}/1,149] "
            f"{identifier}: "
            f"{'validated' if success else 'FAILED'}",
            flush=True,
        )

state["active_chunks"] = []
state["completed_utc"] = utc_now()
state["status"] = (
    "complete"
    if len(completed_ids) == 1_149
    and not state["failed_chunks"]
    else "failed"
)
write_state()

total_rows = 0
total_firth = 0
total_failed_tests = 0
chromosome_counts = {}

for task in tasks:
    receipt = validated_receipt(task)

    if receipt is None:
        continue

    chromosome = task["chromosome"]
    rows = int(receipt["n_result_rows"])

    total_rows += rows
    total_firth += int(receipt["firth_tests"])
    total_failed_tests += int(
        receipt["failed_tests"]
    )

    chromosome_counts.setdefault(
        chromosome,
        {
            "chunks": 0,
            "result_rows": 0,
            "failed_tests": 0,
        },
    )

    chromosome_counts[chromosome]["chunks"] += 1
    chromosome_counts[chromosome][
        "result_rows"
    ] += rows
    chromosome_counts[chromosome][
        "failed_tests"
    ] += int(receipt["failed_tests"])

summary = {
    "status": state["status"],
    "attempt_id": plan["attempt_id"],
    "completed_utc": utc_now(),
    "validated_chunks": len(completed_ids),
    "expected_chunks": 1_149,
    "total_result_rows": total_rows,
    "expected_result_rows": 13_255_796,
    "total_firth_tests": total_firth,
    "total_failed_tests": total_failed_tests,
    "failed_chunks": state["failed_chunks"],
    "chromosome_counts": chromosome_counts,
    "participant_level_results": False,
}

atomic_json(SUMMARY_PATH, summary)

remote_provenance = (
    f"{plan['remote_result_root']}/provenance"
)

for local_file in [
    STATE_PATH,
    SUMMARY_PATH,
    PLAN_PATH,
    PILOT_RECEIPT,
]:
    gs_copy(
        local_file,
        f"{remote_provenance}/{local_file.name}",
    )

print()
print("=" * 72, flush=True)
print("STEP 2 COORDINATOR FINISHED", flush=True)
print("=" * 72, flush=True)
print("Status:", summary["status"], flush=True)
print(
    "Validated chunks:",
    f"{summary['validated_chunks']:,}/1,149",
    flush=True,
)
print(
    "Result rows:",
    f"{summary['total_result_rows']:,}",
    flush=True,
)
print(
    "Failed tests:",
    f"{summary['total_failed_tests']:,}",
    flush=True,
)
