#!/usr/bin/env python3
"""Select a class-balanced CAMDA training subset and a backup replacement pool.

Per species, already-staged isolates are kept first (capped at the class quota),
then each binary class is filled with a seeded shuffle. Binary labels match the
pipeline: Susceptible = S (class 1), Intermediate+Resistant = R (class 0).
"""

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path


BINARY_SUSCEPTIBLE = "S"
BINARY_RESISTANT = "R"


def load_skip_list(path):
    if path is None or not path.exists():
        return set()
    skipped = set()
    for raw in path.read_text().splitlines():
        value = raw.strip()
        if value and not value.startswith("#"):
            skipped.add(value)
    return skipped


def read_rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        required = {"genus", "species", "accession", "phenotype"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")
        return list(reader), reader.fieldnames


def binary_class(phenotype):
    if phenotype == "Susceptible":
        return BINARY_SUSCEPTIBLE
    if phenotype in {"Intermediate", "Resistant"}:
        return BINARY_RESISTANT
    raise ValueError(f"Unsupported phenotype: {phenotype}")


def species_key(row):
    return f"{row['genus']}_{row['species']}".lower()


def staged_accessions(assemblies_dir):
    if not assemblies_dir.is_dir():
        return set()
    return {path.stem for path in assemblies_dir.glob("*.fasta") if path.is_file()}


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_lines(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value}\n" for value in values))


def select_for_class(rows, staged, quota, backup_quota, rng):
    staged_rows = sorted(
        (row for row in rows if row["accession"] in staged),
        key=lambda row: row["accession"],
    )
    unstaged_rows = [row for row in rows if row["accession"] not in staged]
    rng.shuffle(unstaged_rows)
    selected = staged_rows[:quota]
    need = quota - len(selected)
    selected.extend(unstaged_rows[:need])
    remaining_unstaged = unstaged_rows[max(need, 0) :]
    extra_staged = staged_rows[quota:]
    backups = remaining_unstaged[:backup_quota]
    if len(backups) < backup_quota:
        backups.extend(extra_staged[: backup_quota - len(backups)])
    return selected, backups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/camda_amr_2025/processed/training_metadata.csv"),
    )
    parser.add_argument(
        "--assemblies-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/assemblies/train"),
    )
    parser.add_argument(
        "--bad-accession-list",
        type=Path,
        default=Path("data/camda_amr_2025/processed/bad_accessions.txt"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/balanced_3_3"),
    )
    parser.add_argument("--quota", type=int, default=3)
    parser.add_argument("--backup-quota", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rows, fieldnames = read_rows(args.metadata)
    skip = load_skip_list(args.bad_accession_list)
    staged = staged_accessions(args.assemblies_dir)
    eligible = [row for row in rows if row["accession"] not in skip]

    grouped = defaultdict(lambda: defaultdict(list))
    for row in eligible:
        grouped[species_key(row)][binary_class(row["phenotype"])].append(row)

    selected = []
    backups = []
    summary_species = {}
    for species in sorted(grouped):
        species_summary = {}
        for label in (BINARY_SUSCEPTIBLE, BINARY_RESISTANT):
            rng = random.Random(f"{args.seed}:{species}:{label}")
            chosen, reserve = select_for_class(
                grouped[species][label], staged, args.quota, args.backup_quota, rng
            )
            if len(chosen) < args.quota:
                raise ValueError(
                    f"{species} class {label} has only {len(chosen)} usable isolates"
                )
            selected.extend(chosen)
            backups.extend(reserve)
            reused = sum(1 for row in chosen if row["accession"] in staged)
            species_summary[label] = {
                "selected": len(chosen),
                "reused_staged": reused,
                "to_stage": len(chosen) - reused,
                "backups": len(reserve),
                "selected_accessions": [row["accession"] for row in chosen],
                "backup_accessions": [row["accession"] for row in reserve],
            }
        summary_species[species] = species_summary

    selected.sort(key=lambda row: (species_key(row), row["accession"]))
    backups.sort(key=lambda row: (species_key(row), row["accession"]))

    output_dir = args.output_dir
    write_csv(output_dir / "training_metadata.csv", selected, fieldnames)
    write_lines(output_dir / "accessions.txt", [row["accession"] for row in selected])
    write_lines(output_dir / "backup_accessions.txt", [row["accession"] for row in backups])
    backup_fields = list(fieldnames) + ["binary_class", "species_key"]
    backup_rows = []
    for row in backups:
        item = dict(row)
        item["binary_class"] = binary_class(row["phenotype"])
        item["species_key"] = species_key(row)
        backup_rows.append(item)
    write_csv(output_dir / "backup_pool.csv", backup_rows, backup_fields)

    reused_total = sum(1 for row in selected if row["accession"] in staged)
    summary = {
        "quota": args.quota,
        "backup_quota": args.backup_quota,
        "seed": args.seed,
        "selected": len(selected),
        "backups": len(backups),
        "reused_staged": reused_total,
        "remaining_jobs": len(selected) - reused_total,
        "skipped": sorted(skip),
        "species": summary_species,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
