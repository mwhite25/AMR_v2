#!/usr/bin/env python3
"""Build accession-level DNABERT train/dev/test CSVs for the 3+3 POC.

Uses an 80/20 train/dev split and copies dev to test. Stratifies on binary
phenotype (Susceptible vs Intermediate+Resistant) plus the numeric species
column. Optionally drops duplicate sequences and caps rows per accession.
"""

import argparse
from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full_dataset_path", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--train_frac", type=float, default=0.80)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--drop_dupes", action="store_true")
    parser.add_argument("--max_per_accession", type=int, default=0)
    return parser.parse_args()


def binary_label(phenotype):
    return 1 if int(phenotype) == 1 else 0


def cap_per_accession(df, max_per_accession, seed):
    if max_per_accession <= 0:
        return df
    parts = []
    for accession, group in df.groupby("accession", sort=True):
        if len(group) > max_per_accession:
            parts.append(group.sample(n=max_per_accession, random_state=seed))
        else:
            parts.append(group)
    return pd.concat(parts, ignore_index=True)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.full_dataset_path).dropna(subset=["sequence"])
    if args.drop_dupes:
        df = df.drop_duplicates(subset=["sequence"], keep="first")
    df = cap_per_accession(df, args.max_per_accession, args.seed)

    accessions = (
        df.groupby("accession")
        .agg(species=("species", "first"), phenotype=("phenotype", "first"))
        .reset_index()
    )
    accessions["binary"] = accessions["phenotype"].map(binary_label)
    accessions["stratum"] = (
        accessions["species"].astype(str) + "_" + accessions["binary"].astype(str)
    )
    holdout_n = int(round(len(accessions) * (1.0 - args.train_frac)))
    if accessions["stratum"].nunique() <= holdout_n and (accessions["stratum"].value_counts() >= 2).all():
        stratify = accessions["stratum"]
    else:
        stratify = accessions["binary"]
        print("stratify=binary_phenotype (species x class holdout is too small)")
    train_accs, dev_accs = train_test_split(
        accessions,
        train_size=args.train_frac,
        random_state=args.seed,
        stratify=stratify,
    )
    splits = {
        "train": set(train_accs["accession"]),
        "dev": set(dev_accs["accession"]),
        "test": set(dev_accs["accession"]),
    }
    for name, accs in splits.items():
        part = df[df["accession"].isin(accs)].copy()
        part.to_csv(out_dir / f"{name}.csv", index=False)
        (out_dir / f"{name}_accessions.txt").write_text(
            "\n".join(sorted(accs)) + "\n"
        )
        print(
            name,
            len(accs),
            "accessions",
            len(part),
            "rows",
            "phenotypes",
            sorted(part["phenotype"].map(binary_label).unique().tolist()),
        )


if __name__ == "__main__":
    main()
