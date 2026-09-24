#!/usr/bin/env python3
"""Stage and validate CAMDA assemblies for the AMR_v2 workflow.

This utility handles prebuilt assemblies without downloading or assembling
reads. Each final file is named ``<accession>.fasta`` as expected by
``data_pipeline/pipeline.py``. SRA/fastp/SPAdes integration can consume the
same manifest in a later stage.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import time
from collections import Counter
from pathlib import Path


REQUIRED_METADATA = {"genus", "species", "accession"}
DNA = set("ACGTNRYKMSWBDHV")
SRA_TOOLS = {
    "prefetch": {"executable": "prefetch", "version_args": ["--version"]},
    "fasterq_dump": {"executable": "fasterq-dump", "version_args": ["--version"]},
    "fastp": {"executable": "fastp", "version_args": ["--version"]},
    "spades": {"executable": "spades.py", "version_args": ["--version"]},
}
PAIRED_LAYOUTS = frozenset({"paired", "paired_sra_1_3"})
FASTQ_EXTENSIONS = (".fastq", ".fastq.gz")


class AccessionFailure(ValueError):
    def __init__(self, classification, message):
        self.classification = classification
        self.message = message
        super().__init__(message)


def read_metadata(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        missing = sorted(REQUIRED_METADATA - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")
        rows = list(reader)
    accessions = [row["accession"].strip() for row in rows]
    duplicates = sorted(accession for accession, count in Counter(accessions).items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate metadata accessions: {duplicates[:5]}")
    return rows


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_fasta(path):
    """Return basic FASTA statistics without requiring Biopython."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("missing or empty file")
    record_count = 0
    total_bases = 0
    contig_ids = set()
    current_id = None
    current_length = 0
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    if current_length == 0:
                        raise ValueError(f"empty sequence for {current_id}")
                    total_bases += current_length
                current_id = line[1:].split()[0]
                if not current_id:
                    raise ValueError(f"empty FASTA header at line {line_number}")
                if current_id in contig_ids:
                    raise ValueError(f"duplicate contig ID: {current_id}")
                contig_ids.add(current_id)
                current_length = 0
                record_count += 1
            else:
                if current_id is None:
                    raise ValueError(f"sequence before header at line {line_number}")
                sequence = line.upper()
                invalid = sorted(set(sequence) - DNA)
                if invalid:
                    raise ValueError(
                        f"invalid nucleotide(s) {invalid} at line {line_number}"
                    )
                current_length += len(sequence)
    if current_id is None or current_length == 0:
        raise ValueError("no non-empty FASTA records")
    total_bases += current_length
    return {"contig_count": record_count, "total_bases": total_bases}


def resolve_tool(name, configured=None):
    executable = SRA_TOOLS[name]["executable"]
    path = Path(configured) if configured else shutil.which(executable)
    if path is None:
        return {"path": "", "version": "", "error": f"{executable} not found"}
    try:
        result = subprocess.run(
            [str(path), *SRA_TOOLS[name]["version_args"]],
            check=True,
            capture_output=True,
            text=True,
        )
        version = (result.stdout or result.stderr).strip().splitlines()[0]
        return {"path": str(path), "version": version, "error": ""}
    except (OSError, subprocess.CalledProcessError) as error:
        return {"path": str(path), "version": "", "error": str(error)}


def discover_tools(args):
    configured = {
        "prefetch": args.prefetch_tool,
        "fasterq_dump": args.fasterq_tool,
        "fastp": args.fastp_tool,
        "spades": args.spades_tool,
    }
    tools = {name: resolve_tool(name, configured[name]) for name in SRA_TOOLS}
    print(json.dumps(tools, indent=2))
    return tools


def command_text(command):
    return " ".join(shlex.quote(str(part)) for part in command)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def nonempty_fastq(path):
    return path.is_file() and path.stat().st_size > 0


def fastq_member(work_dir, accession, member=None):
    suffix = "" if member is None else f"_{member}"
    for extension in FASTQ_EXTENSIONS:
        path = work_dir / f"{accession}{suffix}{extension}"
        if nonempty_fastq(path):
            return path
    return None


def existing_fastq_paths(work_dir, accession):
    paths = []
    for member in (1, 2, 3, None):
        suffix = "" if member is None else f"_{member}"
        for extension in FASTQ_EXTENSIONS:
            path = work_dir / f"{accession}{suffix}{extension}"
            if path.exists():
                paths.append(path)
    return paths


def detect_fastq_layout(work_dir, accession):
    """Detect paired or single FASTQ layout after fasterq-dump --split-files.

    SRA suffixes are read-member numbers, not biological mate numbers. Some
    Illumina runs store biological reads as ``_1`` and ``_3``, with a short
    technical/index ``_2`` that is empty or only a subset of spots. When both
    biological files are present and nonempty, that pair is used and any
    leftover ``_2`` is ignored.
    """
    read1 = fastq_member(work_dir, accession, 1)
    read2 = fastq_member(work_dir, accession, 2)
    read3 = fastq_member(work_dir, accession, 3)
    single = fastq_member(work_dir, accession)
    if read1 is not None and read3 is not None:
        return "paired_sra_1_3", [read1, read3]
    if read1 is not None and read2 is not None:
        return "paired", [read1, read2]
    if single is not None and read1 is None and read2 is None and read3 is None:
        return "single", [single]
    if existing_fastq_paths(work_dir, accession):
        raise AccessionFailure(
            "ambiguous_fastq_layout",
            f"incomplete or ambiguous FASTQ files for {accession}",
        )
    return None, []


def count_fastq_reads(path):
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        lines = sum(1 for _ in handle)
    if lines == 0 or lines % 4 != 0:
        return 0
    return lines // 4


def fastq_is_valid(path):
    return count_fastq_reads(path) > 0


def run_command(command, log_path, dry_run):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        log_path.write_text(command_text(command) + "\n")
        return
    with log_path.open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {command_text(command)}")


def stage_sra_row(row, args, tools):
    accession = row["accession"]
    work_dir = args.work_dir / row["dataset"] / accession
    assembly_path = Path(row["expected_path"])
    state_path = work_dir / "state.json"
    provenance_path = work_dir / "provenance.json"
    commands_dir = work_dir / "commands"
    state = {"accession": accession, "dataset": row["dataset"], "state": "planned"}
    if args.resume and state_path.is_file():
        state.update(json.loads(state_path.read_text()))
    if args.resume and assembly_path.is_file():
        try:
            stats = validate_fasta(assembly_path)
            row.update(status="present", sha256=sha256(assembly_path), **stats)
            state["state"] = "staged"
            state.pop("error", None)
            state.pop("classification", None)
            write_json(state_path, state)
            return
        except ValueError:
            pass

    row["source_path"] = f"sra:{accession}"
    row["error"] = ""
    tool_paths = {
        name: tools[name]["path"] or SRA_TOOLS[name]["executable"]
        for name in tools
    }
    provenance = {
        "accession": accession,
        "dataset": row["dataset"],
        "tools": tools,
        "started_at": time.time(),
    }
    write_json(provenance_path, provenance)
    write_json(state_path, state)
    try:
        prefetch_dir = args.sra_cache or work_dir / "sra"
        fastq_dir = work_dir / "fastq"
        trimmed_dir = work_dir / "trimmed"
        spades_dir = work_dir / "spades"
        layout, reads = None, []
        if args.resume:
            try:
                layout, reads = detect_fastq_layout(fastq_dir, accession)
            except AccessionFailure:
                layout, reads = None, []
        if layout is None:
            prefetch = [tool_paths["prefetch"], accession, "-O", str(prefetch_dir)]
            run_command(prefetch, commands_dir / "prefetch.log", args.dry_run)
            if not args.dry_run:
                state["state"] = "downloaded"
            fasterq = [tool_paths["fasterq_dump"], accession, "--split-files", "-O", str(fastq_dir)]
            run_command(fasterq, commands_dir / "fasterq-dump.log", args.dry_run)
            if args.dry_run:
                row.update(status="planned")
                write_json(state_path, state)
                return
            layout, reads = detect_fastq_layout(fastq_dir, accession)
        elif args.dry_run:
            row.update(status="planned")
            write_json(state_path, state)
            return
        if layout is None or not all(fastq_is_valid(path) for path in reads):
            raise AccessionFailure(
                "ambiguous_fastq_layout",
                f"missing or invalid FASTQ output for {accession}",
            )
        state.update(state="fastq_ready", reads_layout=layout)
        state.pop("error", None)
        state.pop("classification", None)
        trimmed_dir.mkdir(parents=True, exist_ok=True)
        if layout in PAIRED_LAYOUTS:
            trimmed_reads = [trimmed_dir / f"{accession}_1.fastq.gz", trimmed_dir / f"{accession}_2.fastq.gz"]
            fastp = [tool_paths["fastp"], "-i", str(reads[0]), "-I", str(reads[1]), "-o", str(trimmed_reads[0]), "-O", str(trimmed_reads[1]), "-j", str(trimmed_dir / "fastp.json"), "-h", str(trimmed_dir / "fastp.html"), "-w", str(args.threads)]
        else:
            trimmed_reads = [trimmed_dir / f"{accession}.fastq.gz"]
            fastp = [tool_paths["fastp"], "-i", str(reads[0]), "-o", str(trimmed_reads[0]), "-j", str(trimmed_dir / "fastp.json"), "-h", str(trimmed_dir / "fastp.html"), "-w", str(args.threads)]
        run_command(fastp, commands_dir / "fastp.log", args.dry_run)
        input_reads = sum(count_fastq_reads(path) for path in reads)
        trimmed_counts = [count_fastq_reads(path) for path in trimmed_reads]
        if not any(trimmed_counts):
            raise AccessionFailure(
                "low_quality_reads",
                f"fastp removed all reads for {accession} ({input_reads} input reads)",
            )
        if input_reads > 0 and sum(trimmed_counts) / input_reads < 0.01:
            raise AccessionFailure(
                "low_quality_reads",
                f"fastp retained only {sum(trimmed_counts)}/{input_reads} reads for {accession}",
            )
        if not all(path.is_file() and path.stat().st_size > 0 for path in trimmed_reads):
            raise AccessionFailure(
                "low_quality_reads",
                f"missing fastp output for {accession}",
            )
        state["state"] = "trimmed"
        spades = [tool_paths["spades"], "--isolate", "-o", str(spades_dir), "-t", str(args.threads)]
        if args.memory:
            spades.extend(["-m", str(args.memory)])
        if layout in PAIRED_LAYOUTS:
            spades.extend(["-1", str(trimmed_reads[0]), "-2", str(trimmed_reads[1])])
        else:
            spades.extend(["-s", str(trimmed_reads[0])])
        run_command(spades, commands_dir / "spades.log", args.dry_run)
        contigs = spades_dir / "contigs.fasta"
        stats = validate_fasta(contigs)
        assembly_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = assembly_path.with_suffix(".fasta.tmp")
        shutil.copy2(contigs, temporary)
        temporary.replace(assembly_path)
        row.update(status="present", sha256=sha256(assembly_path), **stats)
        state["state"] = "staged"
        state.pop("error", None)
        state.pop("classification", None)
    except AccessionFailure as error:
        row.update(status="bad_accession", error=str(error), classification=error.classification)
        state.update(state="bad_accession", error=str(error), classification=error.classification)
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        row.update(status="failed", error=str(error))
        state.update(state="failed", error=str(error))
    finally:
        write_json(state_path, state)
        provenance["finished_at"] = time.time()
        write_json(provenance_path, provenance)


def find_source(source_dir, accession):
    candidates = sorted(
        path
        for path in source_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".fa", ".fna", ".fasta"}
        and accession in path.name
    )
    exact = [path for path in candidates if path.stem == accession]
    if len(exact) == 1:
        return exact[0]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        return None
    raise ValueError(f"Ambiguous assembly sources for {accession}: {candidates[:5]}")


def build_manifest(metadata_paths, assembly_roots):
    rows = []
    for dataset, metadata_path in metadata_paths.items():
        for row in read_metadata(metadata_path):
            accession = row["accession"].strip()
            rows.append(
                {
                    "dataset": dataset,
                    "genus": row["genus"],
                    "species": row["species"],
                    "accession": accession,
                    "expected_path": str(assembly_roots[dataset] / f"{accession}.fasta"),
                    "status": "missing",
                    "source_path": "",
                    "sha256": "",
                    "contig_count": "",
                    "total_bases": "",
                    "error": "",
                    "classification": "",
                }
            )
    return rows


def stage_prebuilt(rows, source_dir, dry_run, resume):
    for row in rows:
        destination = Path(row["expected_path"])
        if resume and destination.is_file():
            try:
                stats = validate_fasta(destination)
                row.update(status="present", sha256=sha256(destination), **stats)
            except ValueError as error:
                row.update(status="invalid", error=str(error))
            continue
        try:
            source = find_source(source_dir, row["accession"])
            if source is None:
                row.update(status="missing", error="source assembly not found")
                continue
            stats = validate_fasta(source)
            row["source_path"] = str(source.resolve())
            row.update(status="available", **stats)
            if dry_run:
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = destination.with_suffix(destination.suffix + ".tmp")
            shutil.copy2(source, temporary)
            temporary.replace(destination)
            row.update(status="present", sha256=sha256(destination))
        except (OSError, ValueError) as error:
            row.update(status="invalid", error=str(error))


def write_manifest(path, rows):
    fields = [
        "dataset", "genus", "species", "accession", "expected_path", "status",
        "source_path", "sha256", "contig_count", "total_bases", "error",
        "classification",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_bad_accession_list(path):
    if path is None:
        return set()
    if not path.exists():
        return set()
    entries = set()
    for raw in path.read_text().splitlines():
        value = raw.strip()
        if not value or value.startswith("#"):
            continue
        entries.add(value)
    return entries


def filter_rows_for_processing(rows, args):
    skip_accessions = set(load_bad_accession_list(args.bad_accession_list))
    skip_accessions.update(
        row["accession"] for row in rows if args.skip_bad_accessions and row.get("status") == "bad_accession"
    )
    if not args.skip_bad_accessions and not skip_accessions:
        return rows
    return [row for row in rows if row["accession"] not in skip_accessions]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["manifest", "stage-prebuilt", "stage-sra", "check-tools", "validate"])
    parser.add_argument("--training-metadata", type=Path, default=Path("data/camda_amr_2025/processed/training_metadata.csv"))
    parser.add_argument("--testing-metadata", type=Path, default=Path("data/camda_amr_2025/processed/testing_metadata.csv"))
    parser.add_argument("--training-dir", type=Path, default=Path("data/camda_amr_2025/processed/assemblies/train"))
    parser.add_argument("--testing-dir", type=Path, default=Path("data/camda_amr_2025/processed/assemblies/test"))
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--work-dir", type=Path, default=Path("data/camda_amr_2025/work/assemblies"))
    parser.add_argument("--sra-cache", type=Path)
    parser.add_argument("--prefetch-tool")
    parser.add_argument("--fasterq-tool")
    parser.add_argument("--fastp-tool")
    parser.add_argument("--spades-tool")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory", type=int)
    parser.add_argument("--dataset", choices=["training", "testing"])
    parser.add_argument("--species")
    parser.add_argument("--accession")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--manifest", type=Path, default=Path("data/camda_amr_2025/processed/assembly_manifest.csv"))
    parser.add_argument("--bad-accession-list", type=Path, help="Optional file with one accession per line to skip from processing")
    parser.add_argument("--skip-bad-accessions", action="store_true", default=True, help="Skip accessions already classified as bad_accession from the manifest")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    if args.mode == "check-tools":
        tools = discover_tools(args)
        raise SystemExit(0 if all(not item["error"] for item in tools.values()) else 1)

    metadata_paths = {"training": args.training_metadata, "testing": args.testing_metadata}
    assembly_roots = {"training": args.training_dir, "testing": args.testing_dir}
    rows = build_manifest(metadata_paths, assembly_roots)
    if args.dataset:
        rows = [row for row in rows if row["dataset"] == args.dataset]
    if args.species:
        rows = [row for row in rows if f"{row['genus']}_{row['species']}".lower() == args.species]
    if args.accession:
        rows = [row for row in rows if row["accession"] == args.accession]
    if args.limit:
        rows = rows[:args.limit]
    processing_rows = filter_rows_for_processing(rows, args)
    if args.mode == "stage-prebuilt":
        if args.source_dir is None:
            parser.error("stage-prebuilt requires --source-dir")
        stage_prebuilt(processing_rows, args.source_dir, args.dry_run, args.resume)
    elif args.mode == "stage-sra":
        tools = discover_tools(args)
        if not args.dry_run and any(item["error"] for item in tools.values()):
            parser.error("stage-sra requires working prefetch, fasterq-dump, fastp, and spades.py")
        for row in processing_rows:
            stage_sra_row(row, args, tools)
    else:
        for row in processing_rows:
            destination = Path(row["expected_path"])
            if destination.is_file():
                try:
                    row.update(status="present", sha256=sha256(destination), **validate_fasta(destination))
                except ValueError as error:
                    row.update(status="invalid", error=str(error))
    write_manifest(args.manifest, rows)
    summary = {
        dataset: {status: sum(row["dataset"] == dataset and row["status"] == status for row in rows)
                  for status in sorted({item["status"] for item in rows if item["dataset"] == dataset})}
        for dataset in metadata_paths
    }
    for dataset in metadata_paths:
        status_names = sorted({item["status"] for item in rows if item["dataset"] == dataset})
        summary[dataset] = {status: sum(row["dataset"] == dataset and row["status"] == status for row in rows)
                            for status in status_names}
    (args.manifest.with_suffix(".json")).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.mode == "validate" and any(row["status"] not in {"present", "bad_accession"} for row in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()