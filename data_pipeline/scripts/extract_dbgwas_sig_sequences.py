#!/usr/bin/env python3
"""Convert DBGWAS textual output into a per-species significant-sequence FASTA.

Prefers unitigs with p <= 0.05. If that set is empty, takes the top unitigs by
p-value (then q-value) up to --top-n. Writes a JSON sidecar describing the rule.
"""

import argparse
import csv
import json
from pathlib import Path


P_COLUMNS = ("p-value", "pvalue", "p_value", "Pvalue", "P-value")
Q_COLUMNS = ("q-value", "qvalue", "q_value", "Qvalue", "Q-value")
SEQ_COLUMNS = ("Sequence", "sequence", "SequenceAllele", "unitig", "Unitig", "seq")
ID_COLUMNS = ("NodeId", "node_id", "NodeID", "Feature", "Id", "ID", "CompId")


def first_column(fieldnames, candidates):
    lookup = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def parse_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"na", "nan", "none", "."}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def find_table(path):
    if path.is_file():
        return path
    if path.is_dir():
        matches = sorted(path.rglob("all_comps_nodes_info.tsv"))
        if matches:
            return matches[0]
    raise FileNotFoundError(f"No all_comps_nodes_info.tsv under {path}")


def read_rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader), reader.fieldnames


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            handle.write(f"{sequence}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="DBGWAS output dir or all_comps_nodes_info.tsv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.output.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite {args.output}")

    table = find_table(args.input)
    rows, fieldnames = read_rows(table)
    p_col = first_column(fieldnames, P_COLUMNS)
    q_col = first_column(fieldnames, Q_COLUMNS)
    seq_col = first_column(fieldnames, SEQ_COLUMNS)
    id_col = first_column(fieldnames, ID_COLUMNS)
    if seq_col is None:
        raise ValueError(f"No sequence column in {table}: {fieldnames}")

    annotated = []
    for index, row in enumerate(rows):
        sequence = (row.get(seq_col) or "").strip()
        if not sequence:
            continue
        annotated.append(
            {
                "header": str(row.get(id_col) or f"{args.species}_unitig_{index}").split()[0],
                "sequence": sequence,
                "p": parse_float(row.get(p_col)) if p_col else None,
                "q": parse_float(row.get(q_col)) if q_col else None,
            }
        )

    significant = [item for item in annotated if item["p"] is not None and item["p"] <= args.p_cutoff]
    rule = f"p<={args.p_cutoff}"
    chosen = significant
    if not chosen:
        ranked = sorted(
            annotated,
            key=lambda item: (
                item["p"] is None,
                item["p"] if item["p"] is not None else 1.0,
                item["q"] is None,
                item["q"] if item["q"] is not None else 1.0,
            ),
        )
        chosen = ranked[: args.top_n]
        rule = f"top_{args.top_n}_by_p"

    records = []
    seen = set()
    for item in chosen:
        header = item["header"]
        suffix = 1
        unique = header
        while unique in seen:
            unique = f"{header}_{suffix}"
            suffix += 1
        seen.add(unique)
        records.append((unique, item["sequence"]))

    write_fasta(args.output, records)
    sidecar = {
        "species": args.species,
        "source": str(table),
        "rule": rule,
        "input_rows": len(rows),
        "with_sequence": len(annotated),
        "significant": len(significant),
        "written": len(records),
        "output": str(args.output),
    }
    args.output.with_suffix(args.output.suffix + ".extract.json").write_text(
        json.dumps(sidecar, indent=2) + "\n"
    )
    print(json.dumps(sidecar, indent=2))
    if not records:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
