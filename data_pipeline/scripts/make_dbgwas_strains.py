#!/usr/bin/env python3
"""Create per-species DBGWAS strains manifests from clean CAMDA metadata."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


PHENOTYPE_MAP = {
    "Susceptible": "1",
    "Intermediate": "0",
    "Resistant": "0",
}


def read_rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        required = {"genus", "species", "accession", "phenotype"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")
        return list(reader)


def write_manifest(rows, assembly_dir, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    species_rows = sorted(rows, key=lambda row: row["accession"])
    missing = []
    invalid = []
    class_values = set()
    with output_path.open("w", newline="") as handle:
        handle.write("ID\tPhenotype\tPath\n")
        for row in species_rows:
            phenotype = PHENOTYPE_MAP.get(row["phenotype"])
            assembly_path = assembly_dir / f"{row['accession']}.fasta"
            if phenotype is None:
                invalid.append(row["accession"])
                continue
            if not assembly_path.is_file():
                missing.append(row["accession"])
                continue
            class_values.add(phenotype)
            handle.write(f"{row['accession']}\t{phenotype}\t{assembly_path.resolve()}\n")
    if invalid:
        raise ValueError(f"Unsupported phenotypes for {output_path}: {invalid}")
    if missing:
        raise FileNotFoundError(
            f"Missing {len(missing)} assemblies for {output_path}: {missing[:5]}"
        )
    if class_values != {"0", "1"}:
        raise ValueError(
            f"Expected both phenotype classes for {output_path}; found {sorted(class_values)}"
        )
    return {"rows": len(species_rows), "classes": sorted(class_values), "path": str(output_path)}


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
        "--output-dir", type=Path, default=Path("data/camda_amr_2025/processed/dbgwas")
    )
    args = parser.parse_args()

    grouped = defaultdict(list)
    for row in read_rows(args.metadata):
        grouped[(row["genus"], row["species"])].append(row)

    reports = []
    for (genus, species), rows in sorted(grouped.items()):
        species_name = f"{genus}_{species}".lower()
        output_path = args.output_dir / f"{species_name}.strains.tsv"
        reports.append(write_manifest(rows, args.assemblies_dir, output_path))

    print(f"Created {len(reports)} DBGWAS strains manifests in {args.output_dir}")
    for report in reports:
        print(f"{report['path']}: {report['rows']} rows, classes={report['classes']}")


if __name__ == "__main__":
    main()