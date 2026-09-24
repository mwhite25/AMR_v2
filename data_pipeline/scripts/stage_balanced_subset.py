#!/usr/bin/env python3
"""Stage a balanced CAMDA subset one accession at a time, prefetching the next.

Only prefetch overlaps the current assemble. fasterq-dump, fastp, and SPAdes
stay exclusive. Failed accessions are replaced from the same-species, same-class
backup pool and appended to the subset metadata.
"""

import argparse
import csv
import json
import shutil
import subprocess
import time
from pathlib import Path


BINARY_MAP = {
    "Susceptible": "S",
    "Intermediate": "R",
    "Resistant": "R",
}

MISSING_FASTA_SPECIES = (
    "streptococcus_pneumoniae",
    "escherichia_coli",
    "staphylococcus_aureus",
)
STRAINS_SCRIPT = Path("data_pipeline/scripts/make_dbgwas_strains.py")
FILL_SCRIPT = Path("data_pipeline/scripts/fill_missing_dbgwas_fastas.py")


def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        return list(reader), list(reader.fieldnames)


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def species_key(row):
    return f"{row['genus']}_{row['species']}".lower()


def binary_class(row):
    phenotype = row.get("binary_class") or BINARY_MAP.get(row["phenotype"])
    if phenotype is None:
        raise ValueError(f"Unsupported phenotype for {row.get('accession')}")
    return phenotype


def is_staged(assemblies_dir, accession):
    path = assemblies_dir / f"{accession}.fasta"
    return path.is_file() and path.stat().st_size > 0


def work_state(work_dir, accession):
    path = work_dir / accession / "state.json"
    if not path.is_file():
        return ""
    try:
        return json.loads(path.read_text()).get("state") or ""
    except (OSError, json.JSONDecodeError):
        return ""


def has_sra_cache(work_dir, accession):
    sra_dir = work_dir / accession / "sra"
    if not sra_dir.is_dir():
        return False
    return any(sra_dir.rglob("*"))


def queue_sort_key(row):
    key = species_key(row)
    try:
        priority = MISSING_FASTA_SPECIES.index(key)
    except ValueError:
        priority = len(MISSING_FASTA_SPECIES)
    return (priority, key, row["accession"])


def append_unique_line(path, value):
    existing = set()
    if path.exists():
        existing = {line.strip() for line in path.read_text().splitlines() if line.strip()}
    if value in existing:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(f"{value}\n")


class PrefetchHandle:
    def __init__(self, accession, process, log_path):
        self.accession = accession
        self.process = process
        self.log_path = log_path

    def wait(self):
        if self.process is None:
            return 0
        return self.process.wait()

    def poll(self):
        if self.process is None:
            return 0
        return self.process.poll()


def start_prefetch(accession, args, logs_dir):
    if is_staged(args.assemblies_dir, accession) or has_sra_cache(args.work_dir, accession):
        return None
    sra_dir = args.work_dir / accession / "sra"
    sra_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"{accession}.prefetch.log"
    command = [args.prefetch_tool, accession, "-O", str(sra_dir)]
    handle = log_path.open("w")
    handle.write(" ".join(command) + "\n")
    handle.flush()
    process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT)
    handle.close()
    return PrefetchHandle(accession, process, log_path)


def wait_prefetch(handle):
    if handle is None:
        return
    handle.wait()


def stage_accession(accession, args, logs_dir):
    manifest = logs_dir / f"{accession}.csv"
    command = [
        args.python,
        str(args.stage_script),
        "stage-sra",
        "--dataset",
        "training",
        "--accession",
        accession,
        "--threads",
        str(args.threads),
        "--memory",
        str(args.memory),
        "--resume",
        "--bad-accession-list",
        str(args.bad_accession_list),
        "--training-metadata",
        str(args.metadata),
        "--training-dir",
        str(args.assemblies_dir),
        "--work-dir",
        str(args.work_root),
        "--manifest",
        str(manifest),
    ]
    log_path = logs_dir / f"{accession}.stage.log"
    started = time.time()
    with log_path.open("w") as log:
        log.write(" ".join(command) + "\n")
        log.flush()
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    finished = time.time()
    status = "present" if is_staged(args.assemblies_dir, accession) else "failed"
    classification = ""
    error = ""
    if manifest.is_file():
        rows, _ = read_csv(manifest)
        if rows:
            status = rows[0].get("status") or status
            classification = rows[0].get("classification") or ""
            error = rows[0].get("error") or ""
    if result.returncode != 0 and status == "present":
        status = "failed"
        error = error or f"stage-sra exited {result.returncode}"
    return {
        "accession": accession,
        "status": status,
        "classification": classification,
        "error": error,
        "returncode": result.returncode,
        "seconds": round(finished - started, 1),
        "log": str(log_path),
    }


def species_complete(selected, species, assemblies_dir):
    rows = [row for row in selected if species_key(row) == species]
    return bool(rows) and all(is_staged(assemblies_dir, row["accession"]) for row in rows)


def maybe_fill_features(species, args, logs_dir):
    strains_dir = args.subset_dir.parent / "dbgwas" / args.subset_dir.name
    strains = [
        args.python,
        str(STRAINS_SCRIPT),
        "--metadata",
        str(args.metadata),
        "--assemblies-dir",
        str(args.assemblies_dir),
        "--output-dir",
        str(strains_dir),
        "--species",
        species,
    ]
    fill = [
        args.python,
        str(FILL_SCRIPT),
        "--metadata",
        str(args.metadata),
        "--assemblies-dir",
        str(args.assemblies_dir),
        "--strains-dir",
        str(strains_dir),
        "--species",
        species,
    ]
    strains_result = subprocess.run(strains, check=False, capture_output=True, text=True)
    fill_result = subprocess.run(fill, check=False, capture_output=True, text=True)
    log_path = logs_dir / f"{species}.features.log"
    log_path.write_text(
        f"strains exit {strains_result.returncode}\n{strains_result.stdout}\n{strains_result.stderr}\n"
        f"fill exit {fill_result.returncode}\n{fill_result.stdout}\n{fill_result.stderr}\n"
    )
    print(f"features_{species} strains={strains_result.returncode} fill={fill_result.returncode}", flush=True)


def next_unfetched(queue, start, args, inflight):
    for row in queue[start:]:
        accession = row["accession"]
        if accession == inflight:
            continue
        if is_staged(args.assemblies_dir, accession):
            continue
        if has_sra_cache(args.work_dir, accession):
            continue
        return accession
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--subset-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/balanced_3_3"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/camda_amr_2025/processed/balanced_3_3/training_metadata.csv"),
    )
    parser.add_argument(
        "--backup-pool",
        type=Path,
        default=Path("data/camda_amr_2025/processed/balanced_3_3/backup_pool.csv"),
    )
    parser.add_argument(
        "--assemblies-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/assemblies/train"),
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path("data/camda_amr_2025/work/assemblies"),
    )
    parser.add_argument(
        "--bad-accession-list",
        type=Path,
        default=Path("data/camda_amr_2025/processed/bad_accessions.txt"),
    )
    parser.add_argument(
        "--stage-script",
        type=Path,
        default=Path("data_pipeline/scripts/stage_camda_assemblies.py"),
    )
    parser.add_argument("--python", default="python")
    parser.add_argument("--prefetch-tool", default="prefetch")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory", type=int, default=16)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.work_dir = args.work_root / "training"
    args.metadata = args.metadata if args.metadata.exists() else args.subset_dir / "training_metadata.csv"
    args.backup_pool = (
        args.backup_pool if args.backup_pool.exists() else args.subset_dir / "backup_pool.csv"
    )

    if shutil.which(args.prefetch_tool) is None:
        raise SystemExit(f"prefetch tool not found: {args.prefetch_tool}")

    selected, fieldnames = read_csv(args.metadata)
    backups, backup_fields = read_csv(args.backup_pool) if args.backup_pool.exists() else ([], [])
    selected_ids = {row["accession"] for row in selected}
    backup_queue = [row for row in backups if row["accession"] not in selected_ids]

    skip_listed = set()
    if args.bad_accession_list.exists():
        skip_listed = {
            line.strip()
            for line in args.bad_accession_list.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        }
    queue = []
    skipped_bad = []
    for row in selected:
        accession = row["accession"]
        if is_staged(args.assemblies_dir, accession):
            continue
        if accession in skip_listed or work_state(args.work_dir, accession) == "bad_accession":
            skipped_bad.append(accession)
            continue
        queue.append(row)
    queue.sort(key=queue_sort_key)
    logs_dir = args.subset_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    prior = args.subset_dir / "stage_driver.json"
    driver_log = []
    if prior.exists():
        try:
            previous = json.loads(prior.read_text())
            if isinstance(previous, list):
                driver_log = previous
            elif isinstance(previous, dict) and isinstance(previous.get("runs"), list):
                driver_log = previous["runs"]
        except (OSError, json.JSONDecodeError):
            driver_log = []
    prefetch_handle = None

    print(
        json.dumps(
            {
                "selected": len(selected),
                "already_staged": sum(is_staged(args.assemblies_dir, row["accession"]) for row in selected),
                "skipped_bad_accession": skipped_bad,
                "to_stage": len(queue),
                "backups": len(backup_queue),
                "dry_run": args.dry_run,
                "queue": [row["accession"] for row in queue],
            },
            indent=2,
        )
    )
    if args.dry_run:
        return

    index = 0
    while index < len(queue):
        row = queue[index]
        accession = row["accession"]
        if prefetch_handle and prefetch_handle.accession == accession:
            wait_prefetch(prefetch_handle)
            prefetch_handle = None
        if prefetch_handle is None:
            upcoming = next_unfetched(
                queue, index + 1, args, inflight=accession
            )
            if upcoming:
                prefetch_handle = start_prefetch(upcoming, args, logs_dir)
                print(f"prefetch_start {upcoming}", flush=True)
        print(f"stage_start {accession}", flush=True)
        result = stage_accession(accession, args, logs_dir)
        driver_log.append(result)
        write_json(args.subset_dir / "stage_driver.json", driver_log)
        print(json.dumps(result), flush=True)
        if result["status"] == "present":
            finished_species = species_key(row)
            if finished_species in MISSING_FASTA_SPECIES and species_complete(
                selected, finished_species, args.assemblies_dir
            ):
                maybe_fill_features(finished_species, args, logs_dir)
        if result["status"] != "present":
            key = species_key(row)
            label = binary_class(row)
            replacement = None
            remaining = []
            for candidate in backup_queue:
                if replacement is None and species_key(candidate) == key and binary_class(candidate) == label:
                    replacement = candidate
                else:
                    remaining.append(candidate)
            backup_queue = remaining
            if replacement is None:
                print(f"no_backup {accession} {key} {label}", flush=True)
            else:
                selected.append(replacement)
                selected_ids.add(replacement["accession"])
                write_csv(args.metadata, selected, fieldnames)
                append_unique_line(args.subset_dir / "accessions.txt", replacement["accession"])
                write_csv(args.backup_pool, backup_queue, backup_fields or list(replacement))
                queue.append(replacement)
                print(f"backup_queued {replacement['accession']}", flush=True)
        index += 1

    if prefetch_handle:
        wait_prefetch(prefetch_handle)

    remaining = [
        row["accession"]
        for row in selected
        if not is_staged(args.assemblies_dir, row["accession"])
        and row["accession"] not in skip_listed
        and work_state(args.work_dir, row["accession"]) != "bad_accession"
    ]
    summary = {
        "attempted": len(driver_log),
        "present": sum(item["status"] == "present" for item in driver_log),
        "failed": sum(item["status"] != "present" for item in driver_log),
        "still_missing": remaining,
    }
    write_json(args.subset_dir / "stage_driver.json", {"runs": driver_log, "summary": summary})
    print(json.dumps(summary, indent=2))
    if remaining:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
