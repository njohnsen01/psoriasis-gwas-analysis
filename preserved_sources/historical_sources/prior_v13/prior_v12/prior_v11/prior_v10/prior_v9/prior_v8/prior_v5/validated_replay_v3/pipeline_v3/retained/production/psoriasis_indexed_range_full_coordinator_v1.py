"""Full psoriasis indexed-range preparation; resumes by chromosome."""

import hashlib
import json
import os
import re
import secrets
import shutil
import sqlite3
import struct
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

BASE = Path("/home/jupyter/psoriasis_indexed_range_production")
CURRENT = BASE / "current_full_production.json"
RUNNER = Path("/home/jupyter/psoriasis_indexed_range_qc_runner_v1.py")
ORDER = [
    "21", "20", "18", "15", "14", "13", "19", "17", "16",
    "9", "11", "12", "10", "8", "X", "5", "6", "7", "3",
    "4", "1", "2",
]


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2))
    tmp.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def command(args, label, output=None, quiet=False):
    if not quiet:
        print(f"\n[{label}]", flush=True)
        print(" ".join(str(x) for x in args), flush=True)
    if output is None:
        result = subprocess.run(
            args, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.stdout and not quiet:
            print(result.stdout, flush=True)
        content = result.stdout
    else:
        with Path(output).open("wb") as handle:
            result = subprocess.run(
                args, stdout=handle, stderr=subprocess.PIPE
            )
        content = result.stderr.decode(errors="replace")
    if result.returncode:
        raise RuntimeError(f"{label} failed:\n{content}")
    return content


def validated_pilot():
    for root in sorted(
        BASE.glob("pilot_chr22_*"),
        key=lambda x: x.stat().st_mtime,
        reverse=True,
    ):
        config = root / "indexed_range_config.json"
        receipts = list(root.rglob("indexed_range_validation.json"))
        for receipt in receipts:
            value = read(receipt)
            if (
                config.exists()
                and value.get("status") == "validated"
                and str(value.get("chromosome")) == "22"
            ):
                return root, config, receipt, value
    raise RuntimeError("No validated chr22 pilot was found.")


def prepare():
    BASE.mkdir(parents=True, exist_ok=True)
    if CURRENT.exists():
        attempt = read(CURRENT)
        state_path = Path(attempt["state"])
        if state_path.exists():
            print("Existing attempt:", attempt["id"])
            print("Status:", read(state_path).get("status"))
            print("Run the launcher block next.")
            return

    pilot_root, pilot_config_path, pilot_receipt, receipt = (
        validated_pilot()
    )
    base = read(pilot_config_path)
    if not RUNNER.exists() or sha(RUNNER) != receipt["runner_sha256"]:
        raise RuntimeError("The pilot-certified runner changed.")

    attempt_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "_" + secrets.token_hex(12)
    )
    root = BASE / f"production_{attempt_id}"
    status = root / "status"
    frozen = root / "frozen_inputs"
    status.mkdir(parents=True)
    frozen.mkdir()

    for source in Path(base["local_keep"]).parent.iterdir():
        if source.is_file() and source.stat().st_size:
            shutil.copy2(source, frozen / source.name)
    base["local_keep"] = str(
        frozen / Path(base["local_keep"]).name
    )
    if base.get("local_sex_update"):
        base["local_sex_update"] = str(
            frozen / Path(base["local_sex_update"]).name
        )
    if not Path(base["local_keep"]).exists():
        raise RuntimeError("Frozen keep file was not preserved.")

    immutable_runner = root / RUNNER.name
    shutil.copy2(RUNNER, immutable_runner)
    base_path = root / "base_config.json"
    write(base_path, base)
    pilot_remote = receipt["remote_output_root"].rstrip("/")
    remote = (
        pilot_remote.rsplit("/", 1)[0]
        + f"/production_{attempt_id}"
    )
    state_path = status / "production_state.json"
    state = {
        "status": "prepared",
        "completed": ["22"],
        "current": None,
        "error": None,
        "updated": now(),
    }
    write(state_path, state)
    attempt = {
        "id": attempt_id,
        "root": str(root),
        "state": str(state_path),
        "pid": str(status / "production.pid"),
        "log": str(status / "production.log"),
        "base": str(base_path),
        "runner": str(immutable_runner),
        "runner_sha256": sha(immutable_runner),
        "project": base["project"],
        "source": base["remote_bgen"].rsplit("/", 1)[0],
        "remote": remote,
        "pilot_root": str(pilot_root),
        "pilot_receipt": str(pilot_receipt),
        "pilot_remote": pilot_remote,
        "order": ORDER,
    }
    write(root / "production_attempt.json", attempt)
    write(CURRENT, attempt)
    print("Prepared attempt:", attempt_id)
    print("Reusing validated chr22:", pilot_remote)
    print("Remaining chromosomes:", len(ORDER))
    print("Run the launcher block next.")


def stat(project, uri):
    text = command(
        ["gsutil", "-u", project, "stat", uri],
        f"Read {uri} metadata", quiet=True,
    )
    size = re.search(r"Content-Length:\s*(\d+)", text)
    generation = re.search(r"Generation:\s*(\d+)", text)
    if not size or not generation:
        raise RuntimeError(f"Cannot parse metadata for {uri}.")
    return int(size.group(1)), generation.group(1)


def stage(project, uri, path):
    path = Path(path)
    if path.exists() and path.stat().st_size:
        return
    partial = path.with_suffix(path.suffix + ".partial")
    partial.unlink(missing_ok=True)
    command(
        ["gsutil", "-u", project, "cp", uri, str(partial)],
        f"Stage {Path(uri).name}",
    )
    partial.replace(path)


def chromosome_config(attempt, chromosome):
    root = Path(attempt["root"]) / f"chr{chromosome}"
    config_path = root / "indexed_range_config.json"
    if config_path.exists():
        return config_path
    metadata = root / "source_metadata"
    status = root / "status"
    metadata.mkdir(parents=True, exist_ok=True)
    status.mkdir(parents=True, exist_ok=True)

    project = attempt["project"]
    bgen = f"{attempt['source']}/chr{chromosome}.bgen"
    bgi_uri = bgen + ".bgi"
    sample_uri = f"{attempt['source']}/chr{chromosome}.sample"
    bgi = metadata / f"chr{chromosome}.bgen.bgi"
    sample = metadata / f"chr{chromosome}.sample"
    header = metadata / f"chr{chromosome}.bgen.header"
    stage(project, bgi_uri, bgi)
    stage(project, sample_uri, sample)
    bgen_bytes, generation = stat(project, bgen)

    with sqlite3.connect(str(bgi)) as connection:
        count, header_bytes, final_byte = connection.execute(
            "SELECT COUNT(*), MIN(file_start_position), "
            "MAX(file_start_position + size_in_bytes) FROM Variant"
        ).fetchone()
    count, header_bytes, final_byte = map(
        int, (count, header_bytes, final_byte)
    )
    if count <= 0 or header_bytes < 16 or final_byte != bgen_bytes:
        raise RuntimeError(f"chr{chromosome} BGI validation failed.")

    if not header.exists() or header.stat().st_size != header_bytes:
        header.unlink(missing_ok=True)
        command(
            [
                "gsutil", "-u", project, "cat", "-r",
                f"0-{header_bytes - 1}",
                f"{bgen}#{generation}",
            ],
            f"Stage chr{chromosome} BGEN header",
            output=header,
        )
    with header.open("rb") as handle:
        offset, _, variants, samples = struct.unpack(
            "<IIII", handle.read(16)
        )
    base = read(attempt["base"])
    sample_rows = sum(
        1 for line in sample.open(errors="replace") if line.strip()
    ) - 2
    if (
        offset + 4 != header_bytes
        or variants != count
        or samples != int(base["expected_source_samples"])
        or sample_rows != int(base["expected_source_samples"])
    ):
        raise RuntimeError(
            f"chr{chromosome} source validation failed."
        )

    remote = f"{attempt['remote']}/chr{chromosome}"
    config = dict(base)
    config.update({
        "chromosome": chromosome,
        "remote_bgen": bgen,
        "remote_bgen_generation": generation,
        "remote_bgi": bgi_uri,
        "remote_sample": sample_uri,
        "local_bgi": str(bgi),
        "local_sample": str(sample),
        "local_header": str(header),
        "expected_header_bytes": header_bytes,
        "expected_source_variants": count,
        "expected_remote_bgen_bytes": bgen_bytes,
        "pilot_root": str(root),
        "stage_root": str(root / "temporary_stage"),
        "state_path": str(status / "indexed_range_state.json"),
        "validation_path": str(
            status / "indexed_range_validation.json"
        ),
        "pending_validation_path": str(
            status / "indexed_range_validation.pending.json"
        ),
        "pid_path": str(status / "indexed_range.pid"),
        "lock_path": str(status / "indexed_range.lock"),
        "remote_output_root": remote,
        "produce_step1_artifacts": chromosome != "X",
    })
    config["input_sha256"] = {
        "bgi": sha(bgi),
        "sample": sha(sample),
        "keep": sha(config["local_keep"]),
        "header": sha(header),
    }
    if chromosome == "X":
        sex_update = config.get("local_sex_update")
        if not sex_update:
            raise RuntimeError("Frozen chrX sex file is missing.")
        config["input_sha256"]["sex_update"] = sha(sex_update)
    write(config_path, config)
    command(
        [
            "gsutil", "-u", project, "cp",
            str(config_path), remote + "/",
        ],
        f"Upload chr{chromosome} frozen configuration",
    )
    return config_path


def receipt_for(attempt, chromosome):
    path = (
        Path(attempt["root"]) / f"chr{chromosome}" / "status"
        / "indexed_range_validation.json"
    )
    if not path.exists():
        return None
    value = read(path)
    if (
        value.get("status") == "validated"
        and str(value.get("chromosome")) == chromosome
        and value.get("runner_sha256") == attempt["runner_sha256"]
        and value.get("staging_removed")
    ):
        return value
    return None


def run():
    attempt = read(CURRENT)
    if sha(attempt["runner"]) != attempt["runner_sha256"]:
        raise RuntimeError("The immutable runner changed.")
    state_path = Path(attempt["state"])
    state = read(state_path)
    state.update(status="running", error=None, updated=now())
    write(state_path, state)
    Path(attempt["pid"]).write_text(str(os.getpid()))
    records = {"22": read(attempt["pilot_receipt"])}
    try:
        for chromosome in attempt["order"]:
            existing = receipt_for(attempt, chromosome)
            if existing is not None:
                records[chromosome] = existing
                if chromosome not in state["completed"]:
                    state["completed"].append(chromosome)
                write(state_path, state)
                continue
            state.update(current=chromosome, updated=now())
            write(state_path, state)
            print(
                f"\n{'=' * 25} CHROMOSOME {chromosome} "
                f"{'=' * 25}",
                flush=True,
            )
            config_path = chromosome_config(attempt, chromosome)
            config = read(config_path)
            lock = Path(config["lock_path"])
            if lock.exists():
                lock.rmdir()
            lock.mkdir()
            code = subprocess.run(
                [sys.executable, attempt["runner"], str(config_path)]
            ).returncode
            if code:
                detail = read(config["state_path"]).get("error", "")
                raise RuntimeError(
                    f"chr{chromosome} failed ({code}): {detail}"
                )
            value = receipt_for(attempt, chromosome)
            if value is None:
                raise RuntimeError(
                    f"chr{chromosome} lacks a validated receipt."
                )
            records[chromosome] = value
            state["completed"].append(chromosome)
            state["updated"] = now()
            write(state_path, state)

        expected = {str(x) for x in range(1, 23)} | {"X"}
        if set(records) != expected:
            raise RuntimeError("Final chromosome set is incomplete.")
        manifest = {
            "status": "validated",
            "attempt_id": attempt["id"],
            "chr22_reused": True,
            "whole_bgen_downloaded": False,
            "regenie_run": False,
            "chromosomes": records,
            "completed_utc": now(),
        }
        manifest_path = (
            Path(attempt["root"])
            / "indexed_range_full_production_manifest.json"
        )
        write(manifest_path, manifest)
        command(
            [
                "gsutil", "-u", attempt["project"], "cp",
                str(manifest_path), attempt["remote"] + "/",
            ],
            "Upload full-production manifest",
        )
        state.update(
            status="complete", current=None,
            completed=sorted(
                expected,
                key=lambda x: 23 if x == "X" else int(x),
            ),
            updated=now(),
        )
        write(state_path, state)
        print("\nFULL INDEXED-RANGE PREPARATION COMPLETE", flush=True)
    except Exception as error:
        state.update(
            status="failed",
            error=f"{type(error).__name__}: {error}",
            updated=now(),
        )
        write(state_path, state)
        raise
    finally:
        Path(attempt["pid"]).unlink(missing_ok=True)


if len(sys.argv) != 2 or sys.argv[1] not in {"prepare", "run"}:
    raise SystemExit("Use: coordinator.py prepare|run")
prepare() if sys.argv[1] == "prepare" else run()
