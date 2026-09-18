#!/usr/bin/env python3
"""Create placeholder DBGWAS feature FASTA files from the balanced pilot assemblies.

This is a pragmatic pilot-only fallback for the AMR_v2 pipeline: the real
DBGWAS stage should eventually produce the per-species *_sig_sequences.fasta
files. Until then, this script builds a valid FASTA query set directly from
assembled contigs so the BLAST step can run against the accepted pilot subset.
"""

import csv
from collections import defaultdict
from pathlib import Path

MANIFEST = Path('data/camda_amr_2025/processed/class_balanced_pilot_manifest.csv')
OUT_DIR = Path('data/dbgwas/p0.05/per_species')


def read_fasta(path):
    records = []
    current_id = None
    current_seq = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    records.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                if current_id is None:
                    raise ValueError(f'No FASTA header before sequence in {path}')
                current_seq.append(line.strip())
    if current_id is not None:
        records.append((current_id, ''.join(current_seq)))
    return records


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    species_to_records = defaultdict(list)

    with MANIFEST.open(newline='') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession = row['accession']
            genus = row['genus']
            species = row['species']
            assembly_path = Path(row['expected_path'])
            if not assembly_path.exists():
                raise FileNotFoundError(f'Missing assembly: {assembly_path}')
            species_key = f"{genus}_{species}".lower()
            for contig_id, sequence in read_fasta(assembly_path):
                header = f'{accession}|{contig_id}'
                species_to_records[species_key].append((header, sequence))

    for species_key, records in sorted(species_to_records.items()):
        out_path = OUT_DIR / f'{species_key}_sig_sequences.fasta'
        with out_path.open('w') as handle:
            for header, sequence in records:
                handle.write(f'>{header}\n')
                handle.write(sequence + '\n')
        print(f'Wrote {len(records)} records to {out_path}')


if __name__ == '__main__':
    main()
