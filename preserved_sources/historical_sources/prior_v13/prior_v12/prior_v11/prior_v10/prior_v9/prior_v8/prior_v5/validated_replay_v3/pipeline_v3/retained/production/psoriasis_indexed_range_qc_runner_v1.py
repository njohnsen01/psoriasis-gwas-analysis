"""Disk-bounded indexed-range preparation for the psoriasis GWAS.

This runner derives bounded BGEN byte ranges from the controlled-v8 BGI,
applies the frozen psoriasis variant-QC rules across all source samples,
then selects the frozen analysis cohort and removes dosage information.
Validated hard-call PGEN chunks are checkpointed only to the secure
workspace bucket for later REGENIE Step 2. Autosomes also contribute
LD-pruned BED artifacts for the single autosomal REGENIE Step 1 model.
Chromosome X is sex-corrected and is excluded from Step 1.
"""

import base64
import hashlib
import json
import os
import re
import shutil
import signal
import sqlite3
import struct
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

CONFIG_PATH = Path(sys.argv[1])
CONFIG = json.loads(CONFIG_PATH.read_text())

VERSION = CONFIG["version"]
PROJECT = CONFIG["project"]
CHROMOSOME = str(CONFIG["chromosome"])
if CHROMOSOME not in {
    *{str(value) for value in range(1, 23)},
    "X",
}:
    raise RuntimeError(
        "This runner supports chromosomes 1-22 and X only."
    )

REMOTE_BGEN = CONFIG["remote_bgen"]
REMOTE_BGEN_GENERATION = str(
    CONFIG["remote_bgen_generation"]
)
REMOTE_BGEN_VERSIONED = (
    f"{REMOTE_BGEN}#{REMOTE_BGEN_GENERATION}"
)
BGI_PATH = Path(CONFIG["local_bgi"])
SAMPLE_PATH = Path(CONFIG["local_sample"])
KEEP_PATH = Path(CONFIG["local_keep"])
HEADER_PATH = Path(CONFIG["local_header"])
EXPECTED_HEADER_BYTES = int(CONFIG["expected_header_bytes"])
EXPECTED_SOURCE_SAMPLES = int(
    CONFIG["expected_source_samples"]
)
EXPECTED_ANALYSIS_SAMPLES = int(
    CONFIG["expected_analysis_samples"]
)
EXPECTED_MALES = int(CONFIG.get("expected_males", 0))
EXPECTED_FEMALES = int(CONFIG.get("expected_females", 0))
EXPECTED_KEEP_SET_SHA256 = CONFIG[
    "expected_keep_set_sha256"
]
EXPECTED_SOURCE_VARIANTS = int(
    CONFIG["expected_source_variants"]
)
EXPECTED_REMOTE_BGEN_BYTES = int(
    CONFIG["expected_remote_bgen_bytes"]
)
EXPECTED_INPUT_SHA256 = CONFIG["input_sha256"]
CHUNK_TARGET_BYTES = int(CONFIG["chunk_target_bytes"])
CHUNK_HARD_CAP_BYTES = int(
    CONFIG["chunk_hard_cap_bytes"]
)
CHUNK_TARGET_VARIANTS = int(
    CONFIG["chunk_target_variants"]
)
CHUNK_HARD_CAP_VARIANTS = int(
    CONFIG["chunk_hard_cap_variants"]
)
MINIMUM_START_FREE_BYTES = int(
    float(CONFIG["minimum_start_free_gib"]) * 1024**3
)
RUNTIME_FLOOR_BYTES = int(
    float(CONFIG["runtime_floor_gib"]) * 1024**3
)
MAXIMUM_WORKING_BYTES = int(
    float(CONFIG["maximum_working_gib"]) * 1024**3
)
PLINK_THREADS = int(CONFIG["plink_threads"])
EXPECTED_PLINK_VERSION = CONFIG["plink_version"]
DISK_CHECK_PATH = Path(CONFIG["disk_check_path"])
PILOT_ROOT = Path(CONFIG["pilot_root"])
STAGE_ROOT = Path(CONFIG["stage_root"])
STATE_PATH = Path(CONFIG["state_path"])
VALIDATION_PATH = Path(CONFIG["validation_path"])
PENDING_VALIDATION_PATH = Path(
    CONFIG["pending_validation_path"]
)
PID_PATH = Path(CONFIG["pid_path"])
LOCK_PATH = Path(CONFIG["lock_path"])
REMOTE_OUTPUT_ROOT = CONFIG["remote_output_root"]
ALGORITHM_INVARIANTS = CONFIG["algorithm_invariants"]
ALGORITHM_SHA256 = CONFIG["algorithm_sha256"]
PRODUCE_STEP1_ARTIFACTS = bool(
    CONFIG.get(
        "produce_step1_artifacts",
        CHROMOSOME != "X",
    )
)
SEX_UPDATE_PATH = (
    Path(CONFIG["local_sex_update"])
    if CONFIG.get("local_sex_update")
    else None
)
REMOTE_STEP2_CHUNK_ROOT = (
    f"{REMOTE_OUTPUT_ROOT}/secure_step2_chunks"
)

ACTIVE_CHILD = None
INITIAL_FREE_BYTES = None
MINIMUM_FREE_BYTES = None
MAXIMUM_CHUNK_BODY_BYTES = 0
MAXIMUM_CHUNK_FILE_BYTES = 0
CANONICAL_FAM_ORDER_SHA256 = None
EXPECTED_KEEP_ORDER_SHA256 = None
CANONICAL_STEP2_PSAM_SHA256 = None
CANONICAL_STEP2_PSAM_REMOTE = None
STEP2_CHUNK_RECORDS = []
STATE = {}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def stable_json_sha256(payload):
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_json_atomic(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2))
    temporary.replace(path)


def sha256_path(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(16 * 1024**2)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def file_manifest(path):
    path = Path(path)
    sha256_digest = hashlib.sha256()
    md5_digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(16 * 1024**2)
            if not chunk:
                break
            sha256_digest.update(chunk)
            md5_digest.update(chunk)
    return {
        "name": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_digest.hexdigest(),
        "md5_base64": base64.b64encode(
            md5_digest.digest()
        ).decode(),
    }


def normalized_id_pairs(path):
    pairs = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                raise RuntimeError(
                    f"ID file has fewer than two columns at "
                    f"line {line_number}: {path}"
                )
            pairs.append((fields[0], fields[1]))
    if not pairs or len(set(pairs)) != len(pairs):
        raise RuntimeError(
            f"ID pairs are empty or duplicated: {path}"
        )
    return pairs


def id_set_sha256(pairs):
    digest = hashlib.sha256()
    for fid, iid in sorted(pairs):
        digest.update(f"{fid}\t{iid}\n".encode())
    return digest.hexdigest()


def id_order_sha256(pairs):
    digest = hashlib.sha256()
    for fid, iid in pairs:
        digest.update(f"{fid}\t{iid}\n".encode())
    return digest.hexdigest()


def current_free_bytes():
    return shutil.disk_usage(DISK_CHECK_PATH).free


def update_disk_metrics():
    global MINIMUM_FREE_BYTES
    free_bytes = current_free_bytes()
    if MINIMUM_FREE_BYTES is None:
        MINIMUM_FREE_BYTES = free_bytes
    else:
        MINIMUM_FREE_BYTES = min(
            MINIMUM_FREE_BYTES,
            free_bytes,
        )

    STATE["disk_free_gib"] = round(
        free_bytes / 1024**3,
        2,
    )
    STATE["minimum_free_gib"] = round(
        MINIMUM_FREE_BYTES / 1024**3,
        2,
    )
    if INITIAL_FREE_BYTES is not None:
        STATE["maximum_working_gib"] = round(
            max(
                0,
                INITIAL_FREE_BYTES - MINIMUM_FREE_BYTES,
            )
            / 1024**3,
            2,
        )
    STATE["updated_utc"] = utc_now()
    write_json_atomic(STATE_PATH, STATE)
    return free_bytes


def stop_child_group(process):
    if process is None or process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=15)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=15)


def handle_termination(signum, frame):
    stop_child_group(ACTIVE_CHILD)
    raise RuntimeError(
        f"Pilot received termination signal {signum}."
    )


signal.signal(signal.SIGTERM, handle_termination)
signal.signal(signal.SIGINT, handle_termination)


def assert_runtime_limits(stage_label):
    free_bytes = update_disk_metrics()
    working_bytes = max(
        0,
        INITIAL_FREE_BYTES - free_bytes,
    )
    if free_bytes < RUNTIME_FLOOR_BYTES:
        raise RuntimeError(
            f"Disk safety stop during {stage_label}: "
            f"{free_bytes / 1024**3:.1f} GiB free is below "
            f"the {RUNTIME_FLOOR_BYTES / 1024**3:.1f} GiB "
            "runtime floor."
        )
    if working_bytes > MAXIMUM_WORKING_BYTES:
        raise RuntimeError(
            f"Scratch safety stop during {stage_label}: "
            f"{working_bytes / 1024**3:.1f} GiB was added, "
            f"above the {MAXIMUM_WORKING_BYTES / 1024**3:.1f} "
            "GiB pilot limit."
        )


def run_monitored(
    command,
    stage_label,
    output_path=None,
    output_mode="wb",
    maximum_output_bytes=None,
):
    global ACTIVE_CHILD
    assert_runtime_limits(stage_label)
    print(
        "$",
        " ".join(str(part) for part in command),
        flush=True,
    )

    output_handle = None
    try:
        if output_path is not None:
            output_path = Path(output_path)
            output_handle = output_path.open(output_mode)

        ACTIVE_CHILD = subprocess.Popen(
            [str(part) for part in command],
            stdout=(
                output_handle
                if output_handle is not None
                else None
            ),
            stderr=None,
            start_new_session=True,
        )

        last_state_write = 0.0
        while ACTIVE_CHILD.poll() is None:
            free_bytes = current_free_bytes()
            working_bytes = max(
                0,
                INITIAL_FREE_BYTES - free_bytes,
            )
            if (
                output_path is not None
                and maximum_output_bytes is not None
                and output_path.exists()
                and output_path.stat().st_size
                > maximum_output_bytes
            ):
                stop_child_group(ACTIVE_CHILD)
                raise RuntimeError(
                    f"{stage_label} exceeded its exact output "
                    "byte limit."
                )
            if (
                free_bytes < RUNTIME_FLOOR_BYTES
                or working_bytes > MAXIMUM_WORKING_BYTES
            ):
                stop_child_group(ACTIVE_CHILD)
                assert_runtime_limits(stage_label)

            now = time.monotonic()
            if now - last_state_write >= 5:
                update_disk_metrics()
                last_state_write = now
            time.sleep(1)

        return_code = ACTIVE_CHILD.returncode
        if (
            output_path is not None
            and maximum_output_bytes is not None
            and output_path.exists()
            and output_path.stat().st_size
            > maximum_output_bytes
        ):
            raise RuntimeError(
                f"{stage_label} exceeded its exact output byte "
                "limit at process completion."
            )
        assert_runtime_limits(stage_label)
        if return_code != 0:
            raise RuntimeError(
                f"{stage_label} failed with exit code "
                f"{return_code}."
            )
    finally:
        if output_handle is not None:
            output_handle.close()
        ACTIVE_CHILD = None
        update_disk_metrics()


def count_nonempty_lines(path):
    return sum(
        1
        for line in Path(path).open()
        if line.strip()
    )


def copy_binary_with_watchdog(
    source_path,
    target_path,
    mode,
    stage_label,
    skip_bytes=0,
):
    source_path = Path(source_path)
    target_path = Path(target_path)
    additional_bytes = (
        source_path.stat().st_size - int(skip_bytes)
    )
    if additional_bytes < 0:
        raise RuntimeError(
            f"Invalid copy offset during {stage_label}."
        )

    free_before = current_free_bytes()
    predicted_free = free_before - additional_bytes
    predicted_working = max(
        0,
        INITIAL_FREE_BYTES - predicted_free,
    )
    if (
        predicted_free < RUNTIME_FLOOR_BYTES
        or predicted_working > MAXIMUM_WORKING_BYTES
    ):
        raise RuntimeError(
            f"Preflight disk safety stop during {stage_label}."
        )

    with source_path.open("rb") as source:
        source.seek(int(skip_bytes))
        with target_path.open(mode) as target:
            while True:
                chunk = source.read(16 * 1024**2)
                if not chunk:
                    break
                target.write(chunk)
                assert_runtime_limits(stage_label)
    assert_runtime_limits(stage_label)


def remove_prefix_files(prefix):
    prefix = Path(prefix)
    for path in prefix.parent.glob(prefix.name + ".*"):
        path.unlink(missing_ok=True)
    for path in prefix.parent.glob(
        prefix.name + "-temporary.*"
    ):
        path.unlink(missing_ok=True)


def verify_remote_source_unchanged():
    remote_stat = subprocess.run(
        [
            "gsutil",
            "-u",
            PROJECT,
            "stat",
            REMOTE_BGEN,
        ],
        capture_output=True,
        text=True,
    )
    if remote_stat.returncode != 0:
        raise RuntimeError(
            "Could not recheck the source BGEN metadata."
        )
    size_match = re.search(
        r"Content-Length:\s*(\d+)",
        remote_stat.stdout,
    )
    generation_match = re.search(
        r"Generation:\s*(\d+)",
        remote_stat.stdout,
    )
    if (
        size_match is None
        or generation_match is None
        or int(size_match.group(1))
        != EXPECTED_REMOTE_BGEN_BYTES
        or generation_match.group(1)
        != REMOTE_BGEN_GENERATION
    ):
        raise RuntimeError(
            "The source BGEN object changed after Block 24."
        )


def verify_static_inputs():
    global EXPECTED_KEEP_ORDER_SHA256

    if stable_json_sha256(
        ALGORITHM_INVARIANTS
    ) != ALGORITHM_SHA256:
        raise RuntimeError(
            "The algorithm invariant hash is invalid."
        )

    plink_version_result = subprocess.run(
        ["plink2", "--version"],
        capture_output=True,
        text=True,
    )
    observed_plink_version = (
        plink_version_result.stdout.strip()
        or plink_version_result.stderr.strip()
    )
    if (
        plink_version_result.returncode != 0
        or observed_plink_version
        != EXPECTED_PLINK_VERSION
    ):
        raise RuntimeError(
            "The PLINK2 binary changed after Block 24."
        )

    observed_hashes = {
        "bgi": sha256_path(BGI_PATH),
        "sample": sha256_path(SAMPLE_PATH),
        "keep": sha256_path(KEEP_PATH),
        "header": sha256_path(HEADER_PATH),
    }
    if CHROMOSOME == "X":
        if (
            SEX_UPDATE_PATH is None
            or not SEX_UPDATE_PATH.exists()
            or SEX_UPDATE_PATH.stat().st_size == 0
        ):
            raise RuntimeError(
                "Chromosome X requires the frozen PLINK sex-update file."
            )
        observed_hashes["sex_update"] = sha256_path(
            SEX_UPDATE_PATH
        )
    if observed_hashes != EXPECTED_INPUT_SHA256:
        raise RuntimeError(
            "An indexed-range input changed after preparation."
        )

    keep_pairs = normalized_id_pairs(KEEP_PATH)
    if (
        len(keep_pairs) != EXPECTED_ANALYSIS_SAMPLES
        or id_set_sha256(keep_pairs)
        != EXPECTED_KEEP_SET_SHA256
    ):
        raise RuntimeError(
            "The keep file does not match the certified cohort."
        )
    EXPECTED_KEEP_ORDER_SHA256 = id_order_sha256(keep_pairs)

    verify_remote_source_unchanged()


def build_chunk_plan():
    chunks = []
    total_rows = 0
    total_indexed_bytes = 0
    previous_block_end = None
    previous_position = None
    current = None
    current_last_position = None

    with sqlite3.connect(str(BGI_PATH)) as connection:
        cursor = connection.execute(
            """
            SELECT
                position,
                file_start_position,
                size_in_bytes
            FROM Variant
            ORDER BY file_start_position
            """
        )

        for position, block_start, block_size in cursor:
            position = int(position)
            block_start = int(block_start)
            block_size = int(block_size)
            block_end = block_start + block_size

            if block_size <= 0:
                raise RuntimeError(
                    "The BGI contains a nonpositive block size."
                )
            if previous_block_end is None:
                if block_start != EXPECTED_HEADER_BYTES:
                    raise RuntimeError(
                        "The first BGI block does not begin "
                        "immediately after the BGEN header."
                    )
            elif block_start != previous_block_end:
                raise RuntimeError(
                    "BGI variant blocks are not byte-contiguous."
                )
            if (
                previous_position is not None
                and position < previous_position
            ):
                raise RuntimeError(
                    "BGI physical order regresses in genomic "
                    "position."
                )

            if current is None:
                current = {
                    "first_start": block_start,
                    "last_end": block_end,
                    "variant_count": 0,
                    "indexed_bytes": 0,
                    "start_position": position,
                    "end_position": position,
                }

            exceeds_target = (
                current["variant_count"] > 0
                and (
                    current["indexed_bytes"] + block_size
                    > CHUNK_TARGET_BYTES
                    or current["variant_count"] + 1
                    > CHUNK_TARGET_VARIANTS
                )
                and position != current_last_position
            )
            if exceeds_target:
                chunks.append(current)
                current = {
                    "first_start": block_start,
                    "last_end": block_end,
                    "variant_count": 0,
                    "indexed_bytes": 0,
                    "start_position": position,
                    "end_position": position,
                }

            current["last_end"] = block_end
            current["variant_count"] += 1
            current["indexed_bytes"] += block_size
            current["end_position"] = position

            if (
                current["indexed_bytes"]
                > CHUNK_HARD_CAP_BYTES
                or current["variant_count"]
                > CHUNK_HARD_CAP_VARIANTS
            ):
                raise RuntimeError(
                    "A same-position group forced a chunk above "
                    "its hard byte or variant-count cap."
                )

            current_last_position = position
            previous_position = position
            previous_block_end = block_end
            total_rows += 1
            total_indexed_bytes += block_size

    if current is not None:
        chunks.append(current)
    if total_rows != EXPECTED_SOURCE_VARIANTS:
        raise RuntimeError(
            "The chunk plan does not cover every BGI variant."
        )
    if previous_block_end != EXPECTED_REMOTE_BGEN_BYTES:
        raise RuntimeError(
            "The final BGI block does not match the BGEN size."
        )
    if not chunks:
        raise RuntimeError("The chunk plan is empty.")

    for index, chunk in enumerate(chunks, start=1):
        chunk["chunk_index"] = index
    return chunks, total_indexed_bytes


def construct_chunk_bgen(chunk, chunk_path):
    global MAXIMUM_CHUNK_BODY_BYTES
    global MAXIMUM_CHUNK_FILE_BYTES

    partial_path = chunk_path.with_suffix(".bgen.partial")
    partial_path.unlink(missing_ok=True)
    shutil.copyfile(HEADER_PATH, partial_path)
    with partial_path.open("r+b") as handle:
        handle.seek(8)
        handle.write(
            struct.pack(
                "<I",
                int(chunk["variant_count"]),
            )
        )

    expected_size = (
        EXPECTED_HEADER_BYTES
        + int(chunk["indexed_bytes"])
    )
    run_monitored(
        [
            "gsutil",
            "-u",
            PROJECT,
            "cat",
            "-r",
            (
                f"{chunk['first_start']}-"
                f"{chunk['last_end'] - 1}"
            ),
            REMOTE_BGEN_VERSIONED,
        ],
        f"streaming chunk {chunk['chunk_index']}",
        output_path=partial_path,
        output_mode="ab",
        maximum_output_bytes=expected_size,
    )
    if partial_path.stat().st_size != expected_size:
        raise RuntimeError(
            "A streamed chunk has the wrong exact byte count."
        )

    with partial_path.open("rb") as handle:
        header = handle.read(16)
    offset, _, variant_count, sample_count = struct.unpack(
        "<IIII",
        header,
    )
    if (
        offset + 4 != EXPECTED_HEADER_BYTES
        or variant_count != int(chunk["variant_count"])
        or sample_count != EXPECTED_SOURCE_SAMPLES
    ):
        raise RuntimeError(
            "A constructed chunk has an invalid BGEN header."
        )

    partial_path.replace(chunk_path)
    MAXIMUM_CHUNK_BODY_BYTES = max(
        MAXIMUM_CHUNK_BODY_BYTES,
        int(chunk["indexed_bytes"]),
    )
    MAXIMUM_CHUNK_FILE_BYTES = max(
        MAXIMUM_CHUNK_FILE_BYTES,
        chunk_path.stat().st_size,
    )


def validate_bed_prefix(
    prefix,
    expected_variants,
    expected_samples=EXPECTED_ANALYSIS_SAMPLES,
    expected_id_set_sha256=EXPECTED_KEEP_SET_SHA256,
):
    prefix = Path(prefix)
    bed_path = Path(str(prefix) + ".bed")
    bim_path = Path(str(prefix) + ".bim")
    fam_path = Path(str(prefix) + ".fam")
    for required_path in [bed_path, bim_path, fam_path]:
        if (
            not required_path.exists()
            or required_path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing BED fileset artifact: {required_path}"
            )

    fam_pairs = normalized_id_pairs(fam_path)
    if len(fam_pairs) != int(expected_samples):
        raise RuntimeError(
            f"BED fileset has the wrong cohort: {prefix}"
        )
    if (
        expected_id_set_sha256 is not None
        and id_set_sha256(fam_pairs)
        != expected_id_set_sha256
    ):
        raise RuntimeError(
            f"BED fileset has the wrong ID set: {prefix}"
        )

    bim_count = count_nonempty_lines(bim_path)
    if bim_count != int(expected_variants):
        raise RuntimeError(
            f"BED fileset has {bim_count} variants instead of "
            f"{expected_variants}: {prefix}"
        )

    bytes_per_variant = (
        int(expected_samples) + 3
    ) // 4
    expected_bed_bytes = (
        3 + bytes_per_variant * int(expected_variants)
    )
    if bed_path.stat().st_size != expected_bed_bytes:
        raise RuntimeError(
            f"BED fileset has the wrong exact byte size: {prefix}"
        )
    with bed_path.open("rb") as handle:
        if handle.read(3) != b"\x6c\x1b\x01":
            raise RuntimeError(
                f"BED fileset is not SNP-major: {prefix}"
            )

    validation_prefix = Path(str(prefix) + "_readback")
    run_monitored(
        [
            "plink2",
            "--bfile",
            str(prefix),
            "--write-snplist",
            "--out",
            str(validation_prefix),
        ],
        f"reading BED fileset {prefix.name}",
    )
    validation_snplist = Path(
        str(validation_prefix) + ".snplist"
    )
    if (
        not validation_snplist.exists()
        or count_nonempty_lines(validation_snplist)
        != int(expected_variants)
    ):
        raise RuntimeError(
            f"PLINK could not read every BED variant: {prefix}"
        )
    remove_prefix_files(validation_prefix)
    return fam_pairs


def append_candidate_bed(
    candidate_prefix,
    candidate_variant_count,
    merged_prefix,
    seen_candidate_ids,
    prior_last_position,
):
    global CANONICAL_FAM_ORDER_SHA256

    candidate_prefix = Path(candidate_prefix)
    merged_prefix = Path(merged_prefix)
    candidate_bed = Path(str(candidate_prefix) + ".bed")
    candidate_bim = Path(str(candidate_prefix) + ".bim")
    candidate_fam = Path(str(candidate_prefix) + ".fam")
    merged_bed = Path(str(merged_prefix) + ".bed")
    merged_bim = Path(str(merged_prefix) + ".bim")
    merged_fam = Path(str(merged_prefix) + ".fam")

    fam_pairs = validate_bed_prefix(
        candidate_prefix,
        candidate_variant_count,
    )
    candidate_fam_order_sha256 = id_order_sha256(fam_pairs)
    if CANONICAL_FAM_ORDER_SHA256 is None:
        CANONICAL_FAM_ORDER_SHA256 = (
            candidate_fam_order_sha256
        )
        copy_binary_with_watchdog(
            candidate_fam,
            merged_fam,
            "wb",
            "copying canonical candidate FAM",
        )
    elif (
        candidate_fam_order_sha256
        != CANONICAL_FAM_ORDER_SHA256
    ):
        raise RuntimeError(
            "Candidate BED participant order changed between "
            "chunks."
        )

    candidate_bim_lines = [
        line.rstrip()
        for line in candidate_bim.open()
        if line.strip()
    ]
    if len(candidate_bim_lines) != candidate_variant_count:
        raise RuntimeError(
            "Candidate BIM count changed during append."
        )

    last_position = prior_last_position
    for line in candidate_bim_lines:
        fields = line.split()
        if len(fields) < 6:
            raise RuntimeError("A candidate BIM row is invalid.")
        normalized_chromosome = re.sub(
            r"^chr",
            "",
            fields[0],
            flags=re.IGNORECASE,
        )
        if normalized_chromosome != CHROMOSOME:
            raise RuntimeError(
                "A candidate BIM contains the wrong chromosome."
            )
        position = int(fields[3])
        if (
            last_position is not None
            and position < last_position
        ):
            raise RuntimeError(
                "Candidate variants are not in genomic order."
            )
        variant_id = fields[1]
        if variant_id in seen_candidate_ids:
            raise RuntimeError(
                "A candidate variant ID occurs more than once."
            )
        seen_candidate_ids.add(variant_id)
        last_position = position

    if not merged_bed.exists():
        copy_binary_with_watchdog(
            candidate_bed,
            merged_bed,
            "wb",
            "starting merged candidate BED",
        )
        with merged_bim.open("w") as target:
            for line_index, line in enumerate(
                candidate_bim_lines,
                start=1,
            ):
                target.write(line + "\n")
                if line_index % 10_000 == 0:
                    assert_runtime_limits(
                        "writing merged candidate BIM"
                    )
    else:
        with candidate_bed.open("rb") as source:
            if source.read(3) != b"\x6c\x1b\x01":
                raise RuntimeError(
                    "A candidate BED is not SNP-major."
                )
        copy_binary_with_watchdog(
            candidate_bed,
            merged_bed,
            "ab",
            "appending merged candidate BED",
            skip_bytes=3,
        )
        with merged_bim.open("a") as target:
            for line_index, line in enumerate(
                candidate_bim_lines,
                start=1,
            ):
                target.write(line + "\n")
                if line_index % 10_000 == 0:
                    assert_runtime_limits(
                        "appending merged candidate BIM"
                    )

    remove_prefix_files(candidate_prefix)
    assert_runtime_limits("after candidate BED append")
    return last_position


def upload_secure_object_once(local_path, remote_uri):
    """Upload one controlled intermediate without silent overwrite."""
    local_path = Path(local_path)
    local_manifest = file_manifest(local_path)
    probe = subprocess.run(
        [
            "gsutil",
            "-u",
            PROJECT,
            "stat",
            remote_uri,
        ],
        capture_output=True,
        text=True,
    )
    if probe.returncode == 0:
        observed = remote_stat_record(remote_uri)
        if (
            observed["size_bytes"]
            != local_manifest["size_bytes"]
            or observed["md5_base64"]
            != local_manifest["md5_base64"]
        ):
            raise RuntimeError(
                "A pre-existing secure chunk object does not match "
                f"the local validated file: {remote_uri}"
            )
    else:
        run_monitored(
            [
                "gsutil",
                "-u",
                PROJECT,
                "-o",
                "GSUtil:parallel_composite_upload_threshold=0",
                "cp",
                str(local_path),
                remote_uri,
            ],
            f"uploading secure Step 2 object {local_path.name}",
        )
        observed = remote_stat_record(remote_uri)
        if (
            observed["size_bytes"]
            != local_manifest["size_bytes"]
            or observed["md5_base64"]
            != local_manifest["md5_base64"]
        ):
            raise RuntimeError(
                "Secure Step 2 object failed remote verification: "
                f"{remote_uri}"
            )
    return {
        **local_manifest,
        "uri": remote_uri,
        "generation": observed["generation"],
        "versioned_uri": observed["versioned_uri"],
    }


def validate_analysis_pgen(prefix, expected_variants):
    """Validate one dosage-free, analysis-cohort PGEN chunk."""
    prefix = Path(prefix)
    pgen = Path(str(prefix) + ".pgen")
    pvar = Path(str(prefix) + ".pvar.zst")
    psam = Path(str(prefix) + ".psam")
    for required in [pgen, pvar, psam]:
        if not required.exists() or required.stat().st_size == 0:
            raise RuntimeError(
                f"Missing analysis PGEN artifact: {required}"
            )

    pairs = normalized_id_pairs(psam)
    if (
        len(pairs) != EXPECTED_ANALYSIS_SAMPLES
        or id_set_sha256(pairs) != EXPECTED_KEEP_SET_SHA256
        or id_order_sha256(pairs) != EXPECTED_KEEP_ORDER_SHA256
    ):
        raise RuntimeError(
            f"Analysis PGEN has the wrong cohort or order: {prefix}"
        )
    if CHROMOSOME == "X":
        with psam.open() as handle:
            header = handle.readline().strip().split()
            normalized_header = [
                value.lstrip("#") for value in header
            ]
            if "SEX" not in normalized_header:
                raise RuntimeError(
                    "Chromosome X PSAM lacks the SEX field."
                )
            sex_index = normalized_header.index("SEX")
            observed_sexes = []
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split()
                observed_sexes.append(fields[sex_index])
        male_count = sum(
            value in {"1", "M", "male", "Male"}
            for value in observed_sexes
        )
        female_count = sum(
            value in {"2", "F", "female", "Female"}
            for value in observed_sexes
        )
        if (
            male_count != EXPECTED_MALES
            or female_count != EXPECTED_FEMALES
            or male_count + female_count
            != EXPECTED_ANALYSIS_SAMPLES
        ):
            raise RuntimeError(
                "Chromosome X PSAM sex counts do not match the "
                "frozen psoriasis cohort."
            )

    check_prefix = Path(str(prefix) + "_integrity")
    run_monitored(
        [
            "plink2",
            "--pfile",
            str(prefix),
            "vzs",
            "--validate",
            "--write-snplist",
            "--out",
            str(check_prefix),
        ],
        f"validating analysis PGEN {prefix.name}",
    )
    snplist = Path(str(check_prefix) + ".snplist")
    if (
        not snplist.exists()
        or count_nonempty_lines(snplist)
        != int(expected_variants)
    ):
        raise RuntimeError(
            f"Analysis PGEN variant count is invalid: {prefix}"
        )

    info = subprocess.run(
        [
            "plink2",
            "--pfile",
            str(prefix),
            "vzs",
            "--pgen-info",
            "--out",
            str(check_prefix) + "_pgen_info",
        ],
        capture_output=True,
        text=True,
    )
    info_text = (info.stdout or "") + (info.stderr or "")
    info_log = Path(str(check_prefix) + "_pgen_info.log")
    if info_log.exists():
        info_text += info_log.read_text(errors="replace")
    if (
        info.returncode != 0
        or re.search(
            r"(?im)^\s*No dosages present\s*$",
            info_text,
        )
        is None
    ):
        raise RuntimeError(
            f"Analysis PGEN is not confirmed dosage-free: {prefix}"
        )
    remove_prefix_files(check_prefix)
    remove_prefix_files(
        Path(str(check_prefix) + "_pgen_info")
    )
    return pairs


def upload_step2_chunk(prefix, chunk, variant_count):
    """Freeze one hard-call analysis chunk in the secure bucket."""
    global CANONICAL_STEP2_PSAM_SHA256
    global CANONICAL_STEP2_PSAM_REMOTE

    prefix = Path(prefix)
    chunk_index = int(chunk["chunk_index"])
    chunk_label = f"{chunk_index:04d}"
    pgen = Path(str(prefix) + ".pgen")
    pvar = Path(str(prefix) + ".pvar.zst")
    psam = Path(str(prefix) + ".psam")

    current_psam_sha256 = sha256_path(psam)
    if CANONICAL_STEP2_PSAM_SHA256 is None:
        CANONICAL_STEP2_PSAM_SHA256 = current_psam_sha256
        CANONICAL_STEP2_PSAM_REMOTE = (
            upload_secure_object_once(
                psam,
                (
                    f"{REMOTE_STEP2_CHUNK_ROOT}/"
                    f"chr{CHROMOSOME}_pso_step2.psam"
                ),
            )
        )
    elif current_psam_sha256 != CANONICAL_STEP2_PSAM_SHA256:
        raise RuntimeError(
            "Analysis PSAM changed between indexed ranges."
        )

    pgen_record = upload_secure_object_once(
        pgen,
        (
            f"{REMOTE_STEP2_CHUNK_ROOT}/"
            f"chr{CHROMOSOME}_pso_chunk_{chunk_label}.pgen"
        ),
    )
    pvar_record = upload_secure_object_once(
        pvar,
        (
            f"{REMOTE_STEP2_CHUNK_ROOT}/"
            f"chr{CHROMOSOME}_pso_chunk_{chunk_label}.pvar.zst"
        ),
    )
    record = {
        "chromosome": CHROMOSOME,
        "chunk_index": chunk_index,
        "source_start_position": int(chunk["start_position"]),
        "source_end_position": int(chunk["end_position"]),
        "source_first_byte": int(chunk["first_start"]),
        "source_last_byte_exclusive": int(chunk["last_end"]),
        "source_variant_count": int(chunk["variant_count"]),
        "qc_variant_count": int(variant_count),
        "pgen": pgen_record,
        "pvar": pvar_record,
        "psam_sha256": CANONICAL_STEP2_PSAM_SHA256,
    }
    STEP2_CHUNK_RECORDS.append(record)
    return record


def process_chunk(
    chunk,
    output_root,
    merged_candidate_prefix,
    aggregate_qc_handle,
    seen_qc_ids,
    seen_candidate_ids,
    prior_last_candidate_position,
):
    chunk_index = int(chunk["chunk_index"])
    chunk_label = f"{chunk_index:04d}"
    chunk_bgen = output_root / f"chunk_{chunk_label}.bgen"
    source_qc_prefix = (
        output_root / f"chunk_{chunk_label}_source_qc"
    )
    qc_prefix = (
        output_root / f"chunk_{chunk_label}_analysis_qc"
    )
    qc_snplist_prefix = (
        output_root / f"chunk_{chunk_label}_qc_ids"
    )
    prune_prefix = output_root / f"chunk_{chunk_label}_prune"
    candidate_prefix = (
        output_root / f"chunk_{chunk_label}_candidate"
    )

    STATE["current_stage"] = (
        f"stream_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    STATE["current_chunk"] = chunk_index
    update_disk_metrics()
    construct_chunk_bgen(chunk, chunk_bgen)

    STATE["current_stage"] = (
        f"qc_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    update_disk_metrics()
    # Reproduce the historical psoriasis variant-QC population: filters are
    # calculated across all 414,830 source samples before the 392,155-person
    # analysis cohort is selected. Writing BED here makes the historical
    # hard-call round trip explicit and prevents dosage retention.
    source_qc_command = [
        "plink2",
        "--bgen",
        str(chunk_bgen),
        "ref-first",
        "--sample",
        str(SAMPLE_PATH),
        "--geno",
        "0.02",
        "--maf",
        "0.01",
        "--mac",
        "100",
        "--hwe",
        "1e-15",
        "keep-fewhet",
        "--set-missing-var-ids",
        "@:#$1,$2",
        "--new-id-max-allele-len",
        "1000",
    ]
    if CHROMOSOME == "X":
        source_qc_command.extend(
            ["--split-par", "hg38"]
        )
    source_qc_command.extend([
        "--make-bed",
        "--threads",
        str(PLINK_THREADS),
        "--out",
        str(source_qc_prefix),
    ])
    try:
        run_monitored(
            source_qc_command,
            f"PLINK source-population QC chunk {chunk_index}",
        )
    except RuntimeError:
        qc_log = Path(str(source_qc_prefix) + ".log")
        qc_log_text = (
            qc_log.read_text(errors="replace")
            if qc_log.exists()
            else ""
        )
        if "No variants remaining" not in qc_log_text:
            raise
        chunk_bgen.unlink(missing_ok=True)
        remove_prefix_files(source_qc_prefix)
        remove_prefix_files(qc_prefix)
        remove_prefix_files(qc_snplist_prefix)
        STATE["completed_chunks"] = chunk_index
        STATE["empty_qc_chunks"] = int(
            STATE.get("empty_qc_chunks", 0)
        ) + 1
        update_disk_metrics()
        return (
            0,
            0,
            prior_last_candidate_position,
        )

    source_bim = Path(str(source_qc_prefix) + ".bim")
    source_qc_variant_count = count_nonempty_lines(source_bim)
    if source_qc_variant_count <= 0:
        raise RuntimeError(
            f"Source-population QC produced no BIM rows: {chunk_index}"
        )
    validate_bed_prefix(
        source_qc_prefix,
        source_qc_variant_count,
        expected_samples=EXPECTED_SOURCE_SAMPLES,
        expected_id_set_sha256=None,
    )

    STATE["current_stage"] = (
        f"subset_analysis_cohort_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    update_disk_metrics()
    analysis_command = [
        "plink2",
        "--bfile",
        str(source_qc_prefix),
        "--keep",
        str(KEEP_PATH),
        "--indiv-sort",
        "file",
        str(KEEP_PATH),
    ]
    if CHROMOSOME == "X":
        analysis_command.extend(
            ["--update-sex", str(SEX_UPDATE_PATH)]
        )
    analysis_command.extend([
        "--make-pgen",
        "vzs",
        "erase-dosage",
        "--threads",
        str(PLINK_THREADS),
        "--out",
        str(qc_prefix),
    ])
    run_monitored(
        analysis_command,
        f"write analysis hard-call PGEN chunk {chunk_index}",
    )

    qc_pgen = Path(str(qc_prefix) + ".pgen")
    qc_pvar = Path(str(qc_prefix) + ".pvar.zst")
    qc_psam = Path(str(qc_prefix) + ".psam")
    qc_log = Path(str(qc_prefix) + ".log")
    if (
        any(
            not required_path.exists()
            for required_path in [
                qc_pgen,
                qc_pvar,
                qc_psam,
            ]
        )
        and qc_log.exists()
        and "No variants remaining"
        in qc_log.read_text(errors="replace")
    ):
        chunk_bgen.unlink(missing_ok=True)
        remove_prefix_files(source_qc_prefix)
        remove_prefix_files(qc_prefix)
        remove_prefix_files(qc_snplist_prefix)
        STATE["completed_chunks"] = chunk_index
        STATE["empty_qc_chunks"] = int(
            STATE.get("empty_qc_chunks", 0)
        ) + 1
        update_disk_metrics()
        return (
            0,
            0,
            prior_last_candidate_position,
        )
    for required_path in [
        qc_pgen,
        qc_pvar,
        qc_psam,
    ]:
        if (
            not required_path.exists()
            or required_path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing chunk QC output: {required_path}"
            )

    validate_analysis_pgen(
        qc_prefix,
        source_qc_variant_count,
    )

    # Read the completed PGEN separately because Workbench's Jan 2025 PLINK2
    # build does not permit --write-snplist with --make-pgen erase-dosage.
    STATE["current_stage"] = (
        f"read_qc_ids_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    update_disk_metrics()
    run_monitored(
        [
            "plink2",
            "--pfile",
            str(qc_prefix),
            "vzs",
            "--write-snplist",
            "--out",
            str(qc_snplist_prefix),
        ],
        f"reading QC variant IDs for chunk {chunk_index}",
    )
    qc_snplist = Path(str(qc_snplist_prefix) + ".snplist")
    if (
        not qc_snplist.exists()
        or qc_snplist.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing chunk QC SNP list: {qc_snplist}"
        )

    qc_pairs = normalized_id_pairs(qc_psam)
    if (
        len(qc_pairs) != EXPECTED_ANALYSIS_SAMPLES
        or id_set_sha256(qc_pairs)
        != EXPECTED_KEEP_SET_SHA256
    ):
        raise RuntimeError(
            f"Chunk {chunk_index} retained the wrong cohort."
        )

    qc_validate_prefix = (
        output_root / f"chunk_{chunk_label}_qc_validate"
    )
    run_monitored(
        [
            "plink2",
            "--pfile",
            str(qc_prefix),
            "vzs",
            "--validate",
            "--out",
            str(qc_validate_prefix),
        ],
        f"validating QC PGEN chunk {chunk_index}",
    )
    remove_prefix_files(qc_validate_prefix)

    qc_ids = [
        line.strip()
        for line in qc_snplist.open()
        if line.strip()
    ]
    if not qc_ids or len(set(qc_ids)) != len(qc_ids):
        raise RuntimeError(
            f"Chunk {chunk_index} QC IDs are empty or duplicated."
        )
    for variant_id in qc_ids:
        if variant_id in seen_qc_ids:
            raise RuntimeError(
                "A QC variant ID occurs in multiple chunks."
            )
        seen_qc_ids.add(variant_id)
        aggregate_qc_handle.write(variant_id + "\n")
    aggregate_qc_handle.flush()
    qc_variant_count = len(qc_ids)
    if qc_variant_count != source_qc_variant_count:
        raise RuntimeError(
            "Analysis hard-call PGEN changed the source-QC "
            f"variant count in chunk {chunk_index}."
        )

    STATE["current_stage"] = (
        f"secure_step2_checkpoint_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    update_disk_metrics()
    upload_step2_chunk(
        qc_prefix,
        chunk,
        qc_variant_count,
    )

    if not PRODUCE_STEP1_ARTIFACTS:
        chunk_bgen.unlink(missing_ok=True)
        remove_prefix_files(source_qc_prefix)
        remove_prefix_files(qc_prefix)
        remove_prefix_files(qc_snplist_prefix)
        STATE["completed_chunks"] = chunk_index
        STATE["qc_variants_so_far"] = int(
            STATE.get("qc_variants_so_far", 0)
        ) + qc_variant_count
        update_disk_metrics()
        return (
            qc_variant_count,
            0,
            prior_last_candidate_position,
        )

    STATE["current_stage"] = (
        f"initial_ld_prune_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    update_disk_metrics()
    run_monitored(
        [
            "plink2",
            "--pfile",
            str(qc_prefix),
            "vzs",
            "--indep-pairwise",
            "200kb",
            "1",
            "0.5",
            "--threads",
            str(PLINK_THREADS),
            "--out",
            str(prune_prefix),
        ],
        f"initial LD pruning chunk {chunk_index}",
    )

    prune_in = Path(str(prune_prefix) + ".prune.in")
    candidate_variant_count = (
        count_nonempty_lines(prune_in)
        if prune_in.exists()
        else 0
    )
    if candidate_variant_count <= 0:
        chunk_bgen.unlink(missing_ok=True)
        remove_prefix_files(source_qc_prefix)
        remove_prefix_files(qc_prefix)
        remove_prefix_files(qc_snplist_prefix)
        remove_prefix_files(prune_prefix)
        STATE["completed_chunks"] = chunk_index
        STATE["qc_variants_so_far"] = int(
            STATE.get("qc_variants_so_far", 0)
        ) + qc_variant_count
        STATE["empty_candidate_chunks"] = int(
            STATE.get("empty_candidate_chunks", 0)
        ) + 1
        update_disk_metrics()
        return (
            qc_variant_count,
            0,
            prior_last_candidate_position,
        )

    STATE["current_stage"] = (
        f"write_candidate_chunk_{chunk_index}_of_"
        f"{STATE['planned_chunks']}"
    )
    update_disk_metrics()
    run_monitored(
        [
            "plink2",
            "--pfile",
            str(qc_prefix),
            "vzs",
            "--extract",
            str(prune_in),
            "--make-bed",
            "--threads",
            "1",
            "--out",
            str(candidate_prefix),
        ],
        f"writing candidate BED chunk {chunk_index}",
    )

    next_last_position = append_candidate_bed(
        candidate_prefix,
        candidate_variant_count,
        merged_candidate_prefix,
        seen_candidate_ids,
        prior_last_candidate_position,
    )

    chunk_bgen.unlink(missing_ok=True)
    remove_prefix_files(source_qc_prefix)
    remove_prefix_files(qc_prefix)
    remove_prefix_files(qc_snplist_prefix)
    remove_prefix_files(prune_prefix)

    STATE["completed_chunks"] = chunk_index
    STATE["qc_variants_so_far"] = int(
        STATE.get("qc_variants_so_far", 0)
    ) + qc_variant_count
    STATE["candidate_variants_so_far"] = int(
        STATE.get("candidate_variants_so_far", 0)
    ) + candidate_variant_count
    update_disk_metrics()
    return (
        qc_variant_count,
        candidate_variant_count,
        next_last_position,
    )


def remote_stat_record(remote_uri):
    result = subprocess.run(
        [
            "gsutil",
            "-u",
            PROJECT,
            "stat",
            remote_uri,
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Remote object verification failed: {remote_uri}"
        )
    size_match = re.search(
        r"Content-Length:\s*(\d+)",
        result.stdout,
    )
    md5_match = re.search(
        r"Hash \(md5\):\s*(\S+)",
        result.stdout,
    )
    if size_match is None or md5_match is None:
        raise RuntimeError(
            f"Remote object lacks size or MD5: {remote_uri}"
        )
    generation_match = re.search(
        r"Generation:\s*(\d+)",
        result.stdout,
    )
    if generation_match is None:
        raise RuntimeError(
            f"Remote object lacks a generation: {remote_uri}"
        )
    generation = generation_match.group(1)
    return {
        "size_bytes": int(size_match.group(1)),
        "md5_base64": md5_match.group(1),
        "generation": generation,
        "versioned_uri": f"{remote_uri}#{generation}",
    }


def main():
    global INITIAL_FREE_BYTES
    global MINIMUM_FREE_BYTES

    try:
        initial_free = current_free_bytes()
        if initial_free < MINIMUM_START_FREE_BYTES:
            raise RuntimeError(
                "The indexed-range job requires at least "
                f"{MINIMUM_START_FREE_BYTES / 1024**3:.0f} "
                "GiB free before launch."
            )

        verify_static_inputs()
        INITIAL_FREE_BYTES = initial_free
        MINIMUM_FREE_BYTES = initial_free
        PID_PATH.write_text(str(os.getpid()))
        PENDING_VALIDATION_PATH.unlink(missing_ok=True)

        if STAGE_ROOT.exists():
            shutil.rmtree(STAGE_ROOT)
        output_root = STAGE_ROOT / "outputs"
        output_root.mkdir(parents=True)

        STATE.update({
            "version": VERSION,
            "status": "running",
            "chromosome": CHROMOSOME,
            "produce_step1_artifacts": (
                PRODUCE_STEP1_ARTIFACTS
            ),
            "current_stage": "planning_byte_chunks",
            "current_chunk": None,
            "completed_chunks": 0,
            "empty_qc_chunks": 0,
            "empty_candidate_chunks": 0,
            "initial_free_gib": round(
                initial_free / 1024**3,
                2,
            ),
            "started_utc": utc_now(),
            "updated_utc": utc_now(),
        })
        update_disk_metrics()
    except Exception as error:
        stop_child_group(ACTIVE_CHILD)
        if STAGE_ROOT.exists():
            shutil.rmtree(STAGE_ROOT)
        PENDING_VALIDATION_PATH.unlink(missing_ok=True)
        PID_PATH.unlink(missing_ok=True)
        if LOCK_PATH.exists():
            LOCK_PATH.rmdir()
        STATE["status"] = "failed"
        STATE["error"] = str(error)
        STATE["staging_removed"] = not STAGE_ROOT.exists()
        STATE["current_stage"] = (
            "preflight_failed_and_staging_cleaned"
        )
        STATE["updated_utc"] = utc_now()
        if INITIAL_FREE_BYTES is not None:
            update_disk_metrics()
        else:
            write_json_atomic(STATE_PATH, STATE)
        print(f"ERROR: {error}", flush=True)
        raise

    try:
        chunks, total_indexed_bytes = build_chunk_plan()
        STATE["planned_chunks"] = len(chunks)
        STATE["source_variants"] = EXPECTED_SOURCE_VARIANTS
        STATE["source_indexed_gib"] = round(
            total_indexed_bytes / 1024**3,
            2,
        )
        write_json_atomic(
            PILOT_ROOT / "chunk_plan.json",
            {
                "chromosome": CHROMOSOME,
                "chunk_target_bytes": CHUNK_TARGET_BYTES,
                "chunk_hard_cap_bytes": (
                    CHUNK_HARD_CAP_BYTES
                ),
                "chunk_target_variants": (
                    CHUNK_TARGET_VARIANTS
                ),
                "chunk_hard_cap_variants": (
                    CHUNK_HARD_CAP_VARIANTS
                ),
                "source_variant_count": (
                    EXPECTED_SOURCE_VARIANTS
                ),
                "total_indexed_bytes": total_indexed_bytes,
                "chunks": chunks,
            },
        )
        update_disk_metrics()

        aggregate_qc_snplist = (
            output_root
            / f"chr{CHROMOSOME}_pso_qc_variants.snplist"
        )
        merged_candidate_prefix = (
            output_root
            / f"chr{CHROMOSOME}_pso_candidate_merged"
        )
        seen_qc_ids = set()
        seen_candidate_ids = set()
        total_qc_variants = 0
        total_candidate_variants = 0
        last_candidate_position = None

        with aggregate_qc_snplist.open("w") as qc_handle:
            for chunk in chunks:
                (
                    chunk_qc_count,
                    chunk_candidate_count,
                    last_candidate_position,
                ) = process_chunk(
                    chunk,
                    output_root,
                    merged_candidate_prefix,
                    qc_handle,
                    seen_qc_ids,
                    seen_candidate_ids,
                    last_candidate_position,
                )
                total_qc_variants += chunk_qc_count
                total_candidate_variants += (
                    chunk_candidate_count
                )

        if total_qc_variants <= 0:
            raise RuntimeError(
                "Chromosome-wide psoriasis QC produced no variants."
            )
        if PRODUCE_STEP1_ARTIFACTS:
            if total_candidate_variants <= 0:
                raise RuntimeError(
                    "The autosomal Step 1 candidate set is empty."
                )
        elif (
            total_candidate_variants != 0
            or seen_candidate_ids
        ):
            raise RuntimeError(
                "Chromosome X unexpectedly produced Step 1 "
                "candidate variants."
            )
        if (
            total_qc_variants != len(seen_qc_ids)
            or total_candidate_variants
            != len(seen_candidate_ids)
            or total_candidate_variants
            > total_qc_variants
            or total_qc_variants
            > EXPECTED_SOURCE_VARIANTS
        ):
            raise RuntimeError(
                "Aggregate QC/candidate count invariants failed."
            )

        expected_checkpoint_count = (
            len(chunks)
            - int(STATE.get("empty_qc_chunks", 0))
        )
        if (
            CANONICAL_STEP2_PSAM_REMOTE is None
            or CANONICAL_STEP2_PSAM_SHA256 is None
            or len(STEP2_CHUNK_RECORDS)
            != expected_checkpoint_count
            or sum(
                int(record["qc_variant_count"])
                for record in STEP2_CHUNK_RECORDS
            )
            != total_qc_variants
        ):
            raise RuntimeError(
                "Secure Step 2 checkpoint invariants failed."
            )

        verify_remote_source_unchanged()

        final_prune_variant_count = 0
        final_prune_in = None
        upload_paths = [aggregate_qc_snplist]
        if PRODUCE_STEP1_ARTIFACTS:
            merged_fam_pairs = validate_bed_prefix(
                merged_candidate_prefix,
                total_candidate_variants,
            )
            if (
                id_order_sha256(merged_fam_pairs)
                != EXPECTED_KEEP_ORDER_SHA256
                or CANONICAL_FAM_ORDER_SHA256
                != EXPECTED_KEEP_ORDER_SHA256
            ):
                raise RuntimeError(
                    "Merged candidate FAM order changed."
                )

            final_prune_prefix = (
                output_root
                / f"chr{CHROMOSOME}_pso_global_ldpruned"
            )
            STATE["current_stage"] = "global_ld_pruning"
            STATE["current_chunk"] = None
            update_disk_metrics()
            run_monitored(
                [
                    "plink2",
                    "--bfile",
                    str(merged_candidate_prefix),
                    "--indep-pairwise",
                    "200kb",
                    "1",
                    "0.5",
                    "--threads",
                    str(PLINK_THREADS),
                    "--out",
                    str(final_prune_prefix),
                ],
                (
                    f"global chromosome {CHROMOSOME} "
                    "LD pruning"
                ),
            )

            final_prune_in = Path(
                str(final_prune_prefix) + ".prune.in"
            )
            final_prune_variant_count = (
                count_nonempty_lines(final_prune_in)
                if final_prune_in.exists()
                else 0
            )
            if (
                final_prune_variant_count <= 0
                or final_prune_variant_count
                > total_candidate_variants
            ):
                raise RuntimeError(
                    "Global LD-pruned variant count is invalid."
                )

            step1_prefix = (
                output_root / f"chr{CHROMOSOME}_pso_step1"
            )
            STATE["current_stage"] = "write_final_step1_bed"
            update_disk_metrics()
            run_monitored(
                [
                    "plink2",
                    "--bfile",
                    str(merged_candidate_prefix),
                    "--extract",
                    str(final_prune_in),
                    "--make-bed",
                    "--threads",
                    "1",
                    "--out",
                    str(step1_prefix),
                ],
                (
                    f"writing chromosome {CHROMOSOME} "
                    "Step 1 BED"
                ),
            )
            step1_fam_pairs = validate_bed_prefix(
                step1_prefix,
                final_prune_variant_count,
            )
            if (
                id_set_sha256(step1_fam_pairs)
                != EXPECTED_KEEP_SET_SHA256
                or id_order_sha256(step1_fam_pairs)
                != EXPECTED_KEEP_ORDER_SHA256
            ):
                raise RuntimeError(
                    "Final Step 1 cohort set/order fingerprint "
                    "is invalid."
                )

            remove_prefix_files(merged_candidate_prefix)
            Path(
                str(final_prune_prefix) + ".prune.out"
            ).unlink(missing_ok=True)
            assert_runtime_limits(
                "after removing merged candidate BED"
            )

            upload_paths.extend([
                final_prune_in,
                Path(str(step1_prefix) + ".bed"),
                Path(str(step1_prefix) + ".bim"),
                Path(str(step1_prefix) + ".fam"),
            ])

        step2_manifest_path = (
            output_root
            / f"chr{CHROMOSOME}_pso_step2_chunk_manifest.json"
        )
        step2_manifest = {
            "version": VERSION,
            "chromosome": CHROMOSOME,
            "representation": (
                "hard-call PGEN; dosage-free"
            ),
            "analysis_sample_count": (
                EXPECTED_ANALYSIS_SAMPLES
            ),
            "analysis_id_set_sha256": (
                EXPECTED_KEEP_SET_SHA256
            ),
            "analysis_id_order_sha256": (
                EXPECTED_KEEP_ORDER_SHA256
            ),
            "sex_update_applied": CHROMOSOME == "X",
            "common_psam": CANONICAL_STEP2_PSAM_REMOTE,
            "chunk_count": len(STEP2_CHUNK_RECORDS),
            "total_qc_variants": total_qc_variants,
            "secure_chunk_root": REMOTE_STEP2_CHUNK_ROOT,
            "chunks": STEP2_CHUNK_RECORDS,
        }
        write_json_atomic(
            step2_manifest_path,
            step2_manifest,
        )
        upload_paths.append(step2_manifest_path)

        for path in upload_paths:
            if not path.exists() or path.stat().st_size == 0:
                raise RuntimeError(
                    f"Missing final artifact: {path}"
                )
        artifact_manifest = {
            path.name: file_manifest(path)
            for path in upload_paths
        }

        runner_sha256 = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        config_sha256 = hashlib.sha256(
            CONFIG_PATH.read_bytes()
        ).hexdigest()

        pending_validation = {
            "version": VERSION,
            "status": "pending_remote_verification",
            "chromosome": CHROMOSOME,
            "algorithm_sha256": ALGORITHM_SHA256,
            "algorithm_invariants": ALGORITHM_INVARIANTS,
            "plink_version": EXPECTED_PLINK_VERSION,
            "runner_sha256": runner_sha256,
            "config_sha256": config_sha256,
            "source_bgen_generation": (
                REMOTE_BGEN_GENERATION
            ),
            "source_uris": {
                "bgen": CONFIG["remote_bgen"],
                "bgi": CONFIG["remote_bgi"],
                "sample": CONFIG["remote_sample"],
            },
            "source_input_sha256": dict(
                CONFIG["input_sha256"]
            ),
            "qc_population_sample_count": (
                EXPECTED_SOURCE_SAMPLES
            ),
            "analysis_population_sample_count": (
                EXPECTED_ANALYSIS_SAMPLES
            ),
            "analysis_male_count": (
                EXPECTED_MALES
                if CHROMOSOME == "X"
                else None
            ),
            "analysis_female_count": (
                EXPECTED_FEMALES
                if CHROMOSOME == "X"
                else None
            ),
            "keep_set_sha256": EXPECTED_KEEP_SET_SHA256,
            "analysis_id_order_sha256": (
                EXPECTED_KEEP_ORDER_SHA256
            ),
            "hard_call_roundtrip": True,
            "dosage_retained": False,
            "source_variant_count": (
                EXPECTED_SOURCE_VARIANTS
            ),
            "qc_variant_count": total_qc_variants,
            "produce_step1_artifacts": (
                PRODUCE_STEP1_ARTIFACTS
            ),
            "candidate_variant_count": (
                total_candidate_variants
            ),
            "final_prune_variant_count": (
                final_prune_variant_count
            ),
            "step1_variant_count": (
                final_prune_variant_count
            ),
            "secure_step2_chunk_count": (
                len(STEP2_CHUNK_RECORDS)
            ),
            "secure_step2_common_psam": (
                CANONICAL_STEP2_PSAM_REMOTE
            ),
            "secure_step2_chunk_manifest": (
                step2_manifest_path.name
            ),
            "planned_chunks": len(chunks),
            "completed_chunks": len(chunks),
            "empty_qc_chunks": int(
                STATE.get("empty_qc_chunks", 0)
            ),
            "empty_candidate_chunks": int(
                STATE.get("empty_candidate_chunks", 0)
            ),
            "chunk_target_gib": round(
                CHUNK_TARGET_BYTES / 1024**3,
                3,
            ),
            "chunk_body_hard_cap_gib": round(
                CHUNK_HARD_CAP_BYTES / 1024**3,
                3,
            ),
            "maximum_chunk_body_gib": round(
                MAXIMUM_CHUNK_BODY_BYTES / 1024**3,
                3,
            ),
            "maximum_chunk_file_gib": round(
                MAXIMUM_CHUNK_FILE_BYTES / 1024**3,
                3,
            ),
            "header_gib": round(
                EXPECTED_HEADER_BYTES / 1024**3,
                6,
            ),
            "initial_free_gib": round(
                INITIAL_FREE_BYTES / 1024**3,
                2,
            ),
            "minimum_free_gib": round(
                MINIMUM_FREE_BYTES / 1024**3,
                2,
            ),
            "maximum_working_gib": round(
                (
                    INITIAL_FREE_BYTES
                    - MINIMUM_FREE_BYTES
                )
                / 1024**3,
                2,
            ),
            "runtime_floor_gib": round(
                RUNTIME_FLOOR_BYTES / 1024**3,
                1,
            ),
            "remote_output_root": REMOTE_OUTPUT_ROOT,
            "artifacts": artifact_manifest,
            "staging_removed": False,
            "completed_utc": utc_now(),
        }
        write_json_atomic(
            PENDING_VALIDATION_PATH,
            pending_validation,
        )

        STATE["current_stage"] = "upload_validated_outputs"
        STATE["current_chunk"] = None
        update_disk_metrics()
        run_monitored(
            [
                "gsutil",
                "-u",
                PROJECT,
                "-o",
                (
                    "GSUtil:parallel_composite_upload_"
                    "threshold=0"
                ),
                "cp",
                *[str(path) for path in upload_paths],
                REMOTE_OUTPUT_ROOT + "/",
            ],
            (
                f"uploading validated chromosome "
                f"{CHROMOSOME} outputs"
            ),
        )

        for local_name, manifest_row in (
            artifact_manifest.items()
        ):
            remote_uri = (
                f"{REMOTE_OUTPUT_ROOT}/{local_name}"
            )
            observed_remote = remote_stat_record(remote_uri)
            if (
                observed_remote["size_bytes"]
                != manifest_row["size_bytes"]
                or observed_remote["md5_base64"]
                != manifest_row["md5_base64"]
            ):
                raise RuntimeError(
                    "Remote artifact size/MD5 mismatch: "
                    f"{local_name}"
                )
            manifest_row["generation"] = (
                observed_remote["generation"]
            )
            manifest_row["versioned_uri"] = (
                observed_remote["versioned_uri"]
            )

        pending_validation["artifacts"] = artifact_manifest
        write_json_atomic(
            PENDING_VALIDATION_PATH,
            pending_validation,
        )

        shutil.rmtree(STAGE_ROOT)
        update_disk_metrics()

        final_validation = dict(pending_validation)
        final_validation["status"] = "validated"
        final_validation["staging_removed"] = (
            not STAGE_ROOT.exists()
        )
        final_validation["final_free_gib"] = round(
            current_free_bytes() / 1024**3,
            2,
        )
        final_validation["minimum_free_gib"] = round(
            MINIMUM_FREE_BYTES / 1024**3,
            2,
        )
        final_validation["maximum_working_gib"] = round(
            (
                INITIAL_FREE_BYTES - MINIMUM_FREE_BYTES
            )
            / 1024**3,
            2,
        )
        write_json_atomic(
            PENDING_VALIDATION_PATH,
            final_validation,
        )

        remote_validation_uri = (
            f"{REMOTE_OUTPUT_ROOT}/"
            f"chr{CHROMOSOME}_pso_indexed_range_validation.json"
        )
        run_monitored(
            [
                "gsutil",
                "-u",
                PROJECT,
                "-o",
                (
                    "GSUtil:parallel_composite_upload_"
                    "threshold=0"
                ),
                "cp",
                str(PENDING_VALIDATION_PATH),
                remote_validation_uri,
            ],
            "uploading final indexed-range validation",
        )
        remote_validation_readback = subprocess.run(
            [
                "gsutil",
                "-u",
                PROJECT,
                "cat",
                remote_validation_uri,
            ],
            capture_output=True,
            text=True,
        )
        if (
            remote_validation_readback.returncode != 0
            or json.loads(
                remote_validation_readback.stdout
            )
            != final_validation
        ):
            raise RuntimeError(
                "Remote indexed-range validation readback failed."
            )

        PENDING_VALIDATION_PATH.replace(VALIDATION_PATH)

        STATE["status"] = "complete"
        STATE["current_stage"] = "complete"
        STATE["current_chunk"] = None
        STATE["completed_chunks"] = len(chunks)
        update_disk_metrics()
        print(
            f"Chromosome {CHROMOSOME} indexed-range psoriasis "
            "preparation completed and validated.",
            flush=True,
        )
    except Exception as error:
        stop_child_group(ACTIVE_CHILD)
        if STAGE_ROOT.exists():
            shutil.rmtree(STAGE_ROOT)
        PENDING_VALIDATION_PATH.unlink(missing_ok=True)
        STATE["status"] = "failed"
        STATE["error"] = str(error)
        STATE["staging_removed"] = not STAGE_ROOT.exists()
        STATE["current_stage"] = (
            "failed_and_staging_cleaned"
            if not STAGE_ROOT.exists()
            else "failed_cleanup_incomplete"
        )
        STATE["updated_utc"] = utc_now()
        update_disk_metrics()
        print(f"ERROR: {error}", flush=True)
        raise
    finally:
        PID_PATH.unlink(missing_ok=True)
        if LOCK_PATH.exists():
            LOCK_PATH.rmdir()


if __name__ == "__main__":
    main()
