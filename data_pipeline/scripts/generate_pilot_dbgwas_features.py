#!/usr/bin/env python3
"""Create placeholder DBGWAS feature FASTA files from assembled contigs.

This is a pragmatic fallback for the AMR_v2 pipeline: the real DBGWAS stage
should produce the per-species *_sig_sequences.fasta files. Until then, this
script builds a valid FASTA query set from assemblies so BLAST can run.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    records = []
    current_id = None
    current_seq = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                if current_id is None:
                    raise ValueError(f"No FASTA header before sequence in {path}")
                current_seq.append(line.strip())
    if current_id is not None:
        records.append((current_id, "".join(current_seq)))
    return records


def read_rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        return list(reader)


def assembly_path_for(row, assemblies_dir):
    if row.get("expected_path"):
        return Path(row["expected_path"])
    return assemblies_dir / f"{row['accession']}.fasta"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/camda_amr_2025/processed/class_balanced_pilot_manifest.csv"),
    )
    parser.add_argument(
        "--assemblies-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/assemblies/train"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/dbgwas/p0.05/per_species"),
    )
    parser.add_argument("--species", action="append", default=[], help="Species keys to write. Repeatable.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    requested = {item.strip().lower() for item in args.species if item.strip()}

    species_to_records = defaultdict(list)
    for row in read_rows(args.manifest):
        species = f"{row['genus']}_{row['species']}".lower()
        if requested and species not in requested:
            continue
        assembly_path = assembly_path_for(row, args.assemblies_dir)
        if not assembly_path.exists():
            raise FileNotFoundError(f"Missing assembly: {assembly_path}")
        for contig_id, sequence in read_fasta(assembly_path):
            species_to_records[species].append((f"{row['accession']}|{contig_id}", sequence))

    if requested and not species_to_records:
        raise SystemExit(f"No manifest rows matched --species {sorted(requested)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for species, records in sorted(species_to_records.items()):
        out_path = args.output_dir / f"{species}_sig_sequences.fasta"
        if out_path.exists() and not args.overwrite:
            print(f"skip existing {out_path}")
            continue
        with out_path.open("w") as handle:
            for header, sequence in records:
                handle.write(f">{header}\n")
                handle.write(sequence + "\n")
        print(f"Wrote {len(records)} records to {out_path}")


if __name__ == "__main__":
    main()
