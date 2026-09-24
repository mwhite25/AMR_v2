#!/usr/bin/env python3
"""Fill missing per-species DBGWAS FASTAs: try real DBGWAS, then placeholders.

DBGWAS stays outside Pixi. This wrapper looks for a prebuilt binary, runs it on
a strains TSV, extracts unitigs, and falls back to contig placeholders when
install fails, DBGWAS errors, or the extracted FASTA is empty.
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MISSING = (
    "staphylococcus_aureus",
    "escherichia_coli",
    "streptococcus_pneumoniae",
)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def fasta_records(path):
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    return sum(1 for line in path.read_text().splitlines() if line.startswith(">"))


def find_dbgwas(explicit):
    if explicit:
        path = Path(explicit)
        if path.is_file():
            return str(path)
        return None
    for name in ("DBGWAS", "dbgwas"):
        found = shutil.which(name)
        if found:
            return found
    return None


def run_logged(command, log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write(" ".join(str(part) for part in command) + "\n")
        log.flush()
        return subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)


def assemblies_ready(strains_path):
    if not strains_path.is_file():
        return False
    lines = [line for line in strains_path.read_text().splitlines()[1:] if line.strip()]
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 3 or not Path(parts[2]).is_file():
            return False
    return bool(lines)


def fill_species(species, args, python):
    fasta = args.output_dir / f"{species}_sig_sequences.fasta"
    strains = args.strains_dir / f"{species}.strains.tsv"
    work_dir = args.work_dir / species
    log_dir = args.log_dir / species
    report = {"species": species, "fasta": str(fasta), "method": None}

    if fasta.exists() and not args.overwrite:
        report.update(method="existing", records=fasta_records(fasta))
        return report

    if not assemblies_ready(strains):
        report.update(method="skipped_missing_assemblies", error=f"incomplete strains file: {strains}")
        return report

    binary = find_dbgwas(args.dbgwas)
    if binary:
        work_dir.mkdir(parents=True, exist_ok=True)
        command = [binary, "-strains", str(strains.resolve()), "-output", str(work_dir), "-nb-cores", str(args.cores)]
        result = run_logged(command, log_dir / "dbgwas.log")
        report["dbgwas_returncode"] = result.returncode
        report["dbgwas_binary"] = binary
        if result.returncode == 0:
            extract = [
                python,
                str(SCRIPTS / "extract_dbgwas_sig_sequences.py"),
                "--input",
                str(work_dir),
                "--output",
                str(fasta),
                "--species",
                species,
                "--p-cutoff",
                str(args.p_cutoff),
                "--top-n",
                str(args.top_n),
            ]
            if args.overwrite:
                extract.append("--overwrite")
            extracted = run_logged(extract, log_dir / "extract.log")
            report["extract_returncode"] = extracted.returncode
            if extracted.returncode == 0 and fasta_records(fasta) > 0:
                report.update(method="dbgwas", records=fasta_records(fasta))
                return report
            report["extract_error"] = "empty or failed extract"
        else:
            report["dbgwas_error"] = f"DBGWAS exited {result.returncode}"
    else:
        report["dbgwas_error"] = "DBGWAS binary not found"

    placeholder = [
        python,
        str(SCRIPTS / "generate_pilot_dbgwas_features.py"),
        "--manifest",
        str(args.metadata),
        "--assemblies-dir",
        str(args.assemblies_dir),
        "--output-dir",
        str(args.output_dir),
        "--species",
        species,
    ]
    if args.overwrite:
        placeholder.append("--overwrite")
    placed = run_logged(placeholder, log_dir / "placeholder.log")
    report["placeholder_returncode"] = placed.returncode
    if placed.returncode == 0 and fasta_records(fasta) > 0:
        report.update(method="placeholder", records=fasta_records(fasta))
        return report
    report.update(method="failed", error="placeholder fallback failed")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/camda_amr_2025/processed/balanced_3_3/training_metadata.csv"),
    )
    parser.add_argument(
        "--assemblies-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/assemblies/train"),
    )
    parser.add_argument(
        "--strains-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/dbgwas/balanced_3_3"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/dbgwas/p0.05/per_species"),
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("data/camda_amr_2025/work/dbgwas/balanced_3_3"),
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path("data/camda_amr_2025/processed/balanced_3_3/dbgwas_logs"),
    )
    parser.add_argument("--species", action="append", default=[])
    parser.add_argument("--dbgwas", help="Path to a prebuilt DBGWAS binary, if available")
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    species_list = args.species or list(DEFAULT_MISSING)
    python = sys.executable

    reports = [fill_species(species, args, python) for species in species_list]
    summary = {"species": reports}
    write_json(args.log_dir / "fill_summary.json", summary)
    print(json.dumps(summary, indent=2))
    if any(item.get("method") in {None, "failed"} for item in reports):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
