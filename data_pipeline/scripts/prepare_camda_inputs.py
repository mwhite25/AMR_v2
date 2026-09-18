#!/usr/bin/env python3
"""Prepare CAMDA 2025 metadata and report input prerequisites.

The raw CAMDA CSV files are kept unchanged. This script writes deterministic
derived metadata, duplicate/conflict reports, and optional assembly/DBGWAS
manifests for the AMR_v2 pipeline.
"""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


SPECIES_TO_ANTIBIOTIC = {
    ("Klebsiella", "pneumoniae"): "GEN",
    ("Streptococcus", "pneumoniae"): "ERY",
    ("Escherichia", "coli"): "GEN",
    ("Campylobacter", "jejuni"): "TET",
    ("Salmonella", "enterica"): "GEN",
    ("Neisseria", "gonorrhoeae"): "TET",
    ("Staphylococcus", "aureus"): "ERY",
    ("Pseudomonas", "aeruginosa"): "CAZ",
    ("Acinetobacter", "baumannii"): "CAZ",
}

TRAIN_REQUIRED = {"genus", "species", "accession", "phenotype", "antibiotic"}
TEST_REQUIRED = {"genus", "species", "accession"}
PHENOTYPES = {"Susceptible", "Intermediate", "Resistant"}


def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return reader.fieldnames, list(reader)


def require_columns(path, columns, required):
    missing = sorted(required - set(columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def group_by_accession(rows):
    grouped = defaultdict(list)
    for row in rows:
        accession = row["accession"].strip()
        if not accession:
            raise ValueError("Encountered a training row with an empty accession")
        grouped[accession].append(row)
    return grouped


def conflict_kind(rows):
    phenotype_values = {row["phenotype"] for row in rows}
    measurement_values = {
        (row.get("measurement_sign", ""), row.get("measurement_value", ""))
        for row in rows
    }
    species_values = {
        (row["genus"], row["species"], row["antibiotic"]) for row in rows
    }
    kinds = []
    if len(phenotype_values) > 1:
        kinds.append("phenotype")
    if len(measurement_values) > 1:
        kinds.append("measurement")
    if len(species_values) > 1:
        kinds.append("metadata")
    return ";".join(kinds)


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_training(path, output_dir):
    columns, rows = read_csv(path)
    require_columns(path, columns, TRAIN_REQUIRED)
    for row in rows:
        if row["phenotype"] not in PHENOTYPES:
            raise ValueError(
                f"Unsupported phenotype {row['phenotype']!r} for {row['accession']}"
            )
        expected = SPECIES_TO_ANTIBIOTIC.get((row["genus"], row["species"]))
        if expected is None:
            raise ValueError(f"Unsupported species: {row['genus']} {row['species']}")
        if row["antibiotic"] != expected:
            raise ValueError(
                f"Unexpected antibiotic {row['antibiotic']} for {row['accession']}; "
                f"expected {expected}"
            )

    grouped = group_by_accession(rows)
    conflicts = []
    clean_rows = []
    for accession in sorted(grouped):
        accession_rows = grouped[accession]
        kind = conflict_kind(accession_rows)
        if kind:
            for row in accession_rows:
                conflicts.append({**row, "conflict_kind": kind})
            continue
        clean_rows.append(sorted(accession_rows, key=lambda row: tuple(row.values()))[0])

    normalized_columns = ["genus", "species", "accession", "phenotype", "antibiotic"]
    normalized_columns += [column for column in columns if column not in normalized_columns]
    write_csv(output_dir / "training_metadata.csv", clean_rows, normalized_columns)
    write_csv(output_dir / "training_conflicts.csv", conflicts, normalized_columns + ["conflict_kind"])
    return {
        "raw_rows": len(rows),
        "raw_unique_accessions": len(grouped),
        "clean_rows": len(clean_rows),
        "conflict_accessions": len({row["accession"] for row in conflicts}),
        "conflict_rows": len(conflicts),
    }


def normalize_testing(path, output_dir):
    columns, rows = read_csv(path)
    require_columns(path, columns, TEST_REQUIRED)
    seen = set()
    normalized_rows = []
    for row in rows:
        accession = row["accession"].strip()
        if not accession:
            raise ValueError("Encountered a testing row with an empty accession")
        if accession in seen:
            raise ValueError(f"Duplicate testing accession: {accession}")
        seen.add(accession)
        expected = SPECIES_TO_ANTIBIOTIC.get((row["genus"], row["species"]))
        if expected is None:
            raise ValueError(f"Unsupported species: {row['genus']} {row['species']}")
        normalized_rows.append({**row, "antibiotic": expected})

    normalized_columns = ["genus", "species", "accession", "antibiotic"]
    normalized_columns += [column for column in columns if column not in normalized_columns]
    normalized_rows.sort(key=lambda row: row["accession"])
    write_csv(output_dir / "testing_metadata.csv", normalized_rows, normalized_columns)
    return {"raw_rows": len(rows), "unique_accessions": len(seen)}


def write_summary(output_dir, training_summary, testing_summary, train_path, test_path):
    with train_path.open(newline="") as handle:
        train_ids = {row["accession"] for row in csv.DictReader(handle)}
    with test_path.open(newline="") as handle:
        test_ids = {row["accession"] for row in csv.DictReader(handle)}
    summary = {
        "training": training_summary,
        "testing": testing_summary,
        "train_test_accession_overlap": len(train_ids & test_ids),
        "label_policy": {"Susceptible": 1, "Intermediate": 0, "Resistant": 0},
        "raw_files_unchanged": True,
    }
    (output_dir / "preparation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def write_assembly_manifest(metadata_path, assemblies_dir, output_path, dataset):
    with metadata_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    manifest = []
    for row in rows:
        accession = row["accession"]
        assembly_path = assemblies_dir / f"{accession}.fasta"
        manifest.append(
            {
                "dataset": dataset,
                "accession": accession,
                "assembly_path": str(assembly_path),
                "status": "present" if assembly_path.is_file() else "missing",
            }
        )
    write_csv(output_path, manifest, ["dataset", "accession", "assembly_path", "status"])
    return {
        "dataset": dataset,
        "assembly_dir": str(assemblies_dir),
        "expected": len(manifest),
        "present": sum(row["status"] == "present" for row in manifest),
        "missing": sum(row["status"] == "missing" for row in manifest),
    }


def write_dbgwas_manifest(dbgwas_dir, output_path):
    manifest = []
    for genus, species in sorted(SPECIES_TO_ANTIBIOTIC):
        species_name = f"{genus}_{species}".lower()
        feature_path = dbgwas_dir / f"{species_name}_sig_sequences.fasta"
        manifest.append(
            {
                "species": species_name,
                "feature_path": str(feature_path),
                "status": "present" if feature_path.is_file() else "missing",
            }
        )
    write_csv(output_path, manifest, ["species", "feature_path", "status"])
    return {
        "feature_dir": str(dbgwas_dir),
        "expected": len(manifest),
        "present": sum(row["status"] == "present" for row in manifest),
        "missing": sum(row["status"] == "missing" for row in manifest),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("data/camda_amr_2025"))
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--train-assemblies-dir", type=Path, default=None)
    parser.add_argument("--test-assemblies-dir", type=Path, default=None)
    parser.add_argument(
        "--dbgwas-dir", type=Path, default=Path("data/dbgwas/p0.05/per_species")
    )
    args = parser.parse_args()

    output_dir = args.output_dir or args.input_dir / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)
    training_summary = normalize_training(
        args.input_dir / "training_dataset.csv", output_dir
    )
    testing_summary = normalize_testing(
        args.input_dir / "testing_dataset.csv", output_dir
    )
    train_assemblies_dir = args.train_assemblies_dir or output_dir / "assemblies" / "train"
    test_assemblies_dir = args.test_assemblies_dir or output_dir / "assemblies" / "test"
    assembly_summary = {
        "training": write_assembly_manifest(
            output_dir / "training_metadata.csv",
            train_assemblies_dir,
            output_dir / "training_assembly_manifest.csv",
            "training",
        ),
        "testing": write_assembly_manifest(
            output_dir / "testing_metadata.csv",
            test_assemblies_dir,
            output_dir / "testing_assembly_manifest.csv",
            "testing",
        ),
    }
    dbgwas_summary = write_dbgwas_manifest(
        args.dbgwas_dir, output_dir / "dbgwas_manifest.csv"
    )
    summary = write_summary(
        output_dir,
        training_summary,
        testing_summary,
        output_dir / "training_metadata.csv",
        output_dir / "testing_metadata.csv",
    )
    summary["assemblies"] = assembly_summary
    summary["dbgwas"] = dbgwas_summary
    (output_dir / "preparation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()