# AMR_v2 Pipeline Guide

This guide explains the AMR_v2 workflow from the beginning. It assumes that
you do not yet know what the input files mean or what each program produces.

The repository was written for the CAMDA AMR 2025 challenge. Its intended
prediction is binary antimicrobial susceptibility:

- `Susceptible`
- `Resistant`

The challenge contains nine bacterial species and four antibiotics. The
workflow is computationally large: assembly, BLAST, DBGWAS, DNABERT, and
random-forest steps can each require substantial time, memory, and storage.
Start with one accession before processing a complete dataset.

The original CAMDA work ran on an **HPC** cluster (high-performance
computing: a shared machine with a scheduler such as SLURM and a shared
filesystem). That is why several scripts still default to paths like
`/gpfs/scratch/jvaska/CAMDA_AMR/AMR_v2` and why
`dnabert/finetune/finetune_scripts/` are SLURM wrappers. A local Pixi run
must pass `--base_dir`, `--metadata_dir`, and similar flags so those
defaults are not used. HPC is the original compute environment, not a
service this repository talks to.

## 1. The Big Picture

Each row in the metadata describes one bacterial isolate. An accession such as
`DRR124693` identifies sequencing data in the Sequence Read Archive (SRA); it
does not identify the species by itself. The species comes from the metadata.

The workflow transforms that isolate through these stages:

```text
metadata + SRA accession
        |
        v
raw sequencing reads
        |
        |  fasterq-dump, fastp
        v
trimmed reads
        |
        |  SPAdes --isolate
        v
assembly FASTA: <accession>.fasta
        |
        |  DBGWAS on training assemblies
        v
significant feature FASTA files
        |
        |  makeblastdb + blastn
        v
BLAST hit tables
        |
        |  data_pipeline/pipeline.py
        v
sequence dataset CSV
        |
        |  DNABERT
        v
per-feature sequence predictions
        |
        |  random forest / XGBoost
        v
accession-level phenotype prediction
```

The workflow has two different kinds of data:

1. **Sequence data**, such as reads, assemblies, and DNA feature sequences.
2. **Tabular data**, such as phenotype labels, hit counts, model predictions,
   and accession metadata.

The CSV metadata files alone do not contain the genome sequence needed for
assembly, BLAST, or DNABERT.

## 2. Basic Terms

| Term | Meaning in this workflow |
| --- | --- |
| Isolate | One bacterial sample used as a prediction example. |
| Accession | An identifier for sequencing data, for example `SRR10037257`. |
| SRA | NCBI Sequence Read Archive. Public store of raw sequencing runs. An accession such as `SRR1187804` names one run there. |
| prefetch | NCBI SRA Toolkit command that downloads a run’s `.sra` archive. It is not a web scraper; see section 6.2. |
| HPC | High-performance computing cluster the original CAMDA jobs used (SLURM, `/gpfs/scratch/...`). Local runs must override those paths. |
| Read | A short DNA sequence produced by a sequencing machine. |
| FASTQ | A file containing reads and their quality scores. |
| Assembly | Longer reconstructed DNA sequences made from reads. |
| Contig | One reconstructed sequence in an assembly FASTA file. |
| FASTA | A text format containing sequence records beginning with `>header`. |
| DBGWAS feature | A DNA sequence found to be statistically associated with phenotype. |
| BLAST hit | An alignment between a DBGWAS feature and an assembly contig. |
| Query ID | The identifier of a DBGWAS feature used in a BLAST query. |
| Hit count | The number of BLAST matches for one feature in one accession. |
| DNABERT | The sequence model that predicts a class for each sequence window. |
| RF | Random forest, which combines per-feature values into an accession prediction. |

## 3. Species and Antibiotics

AMR_v2 uses these fixed mappings:

| Species | Antibiotic code |
| --- | --- |
| `klebsiella_pneumoniae` | `GEN` |
| `streptococcus_pneumoniae` | `ERY` |
| `escherichia_coli` | `GEN` |
| `campylobacter_jejuni` | `TET` |
| `salmonella_enterica` | `GEN` |
| `neisseria_gonorrhoeae` | `TET` |
| `staphylococcus_aureus` | `ERY` |
| `pseudomonas_aeruginosa` | `CAZ` |
| `acinetobacter_baumannii` | `CAZ` |

The accession prefix does not determine the species:

```text
DRR124693  -> Neisseria gonorrhoeae
ERR1218574 -> Klebsiella pneumoniae
SRR10037257 -> Pseudomonas aeruginosa
```

Look up an accession in the metadata instead:

```bash
awk -F, '$3 == "SRR10037257" {print $1, $2, $3}' \
  data/camda_amr_2025/processed/training_metadata.csv
```

The repository encodes labels as follows:

```text
Susceptible -> 1
Intermediate -> 0
Resistant -> 0
```

This is the existing binary convention in `data_pipeline/pipeline.py`. Keep
the original text labels in metadata and audit files even when a model uses
the numeric values.

## 4. Important Directories

| Path | Role |
| --- | --- |
| `data/camda_amr_2025/training_dataset.csv` | Raw labeled CAMDA training metadata. |
| `data/camda_amr_2025/testing_dataset.csv` | Raw unlabeled CAMDA testing metadata. |
| `data/camda_amr_2025/processed/` | Derived metadata, reports, and manifests. |
| `data/camda_amr_2025/processed/assemblies/train/` | Final training assemblies. |
| `data/camda_amr_2025/processed/assemblies/test/` | Final testing assemblies. |
| `data/camda_amr_2025/work/assemblies/` | SRA, FASTQ, trimming, SPAdes, logs, and state files. |
| `data/dbgwas/p0.05/per_species/` | Final nine DBGWAS feature FASTA files expected by AMR_v2. |
| `work/blast/` | BLAST databases and hit tables. |
| `work/models/dnabert2-117m/` | Cached public DNABERT-2-117M base checkpoint. Smoke tests only. |
| `dnabert/finetune/train.py` | Hugging Face trainer used to fine-tune DNABERT. |
| `dnabert/finetune/finetune_scripts/` | Historical SLURM launch scripts. They contain HPC paths. |
| `<output_dir>/<run_name>/` | Sequence datasets produced by `pipeline.py`; this repository commonly uses `data/camda_amr_2025/processed/pilot_outputs/`. |
| `dnabert/finetune/data/prepare_dnabert_splits.py` | Accession-level train/dev/test split helper used by the 3+3 POC. |
| `dnabert/inference/outputs/` | DNABERT predictions and RF input datasets. |
| `rf/models/` and `rf/predictions/` | Trained models and final predictions. |
| [AMR_V2_POC_3_3.md](AMR_V2_POC_3_3.md) | Write-up of the class-balanced 3+3 tooling run on this machine. |

Files under `work/` are intermediate files. Files under `processed/` are
metadata, manifests, reports, pilot outputs, and final staged assemblies. Do
not treat a metadata file as proof that an assembly or model output exists;
inspect the corresponding directory and status fields.

### Artifact lifecycle

Keep these artifacts because later stages or audits depend on them:

- raw metadata and the prepared `training_metadata.csv` and
  `testing_metadata.csv`;
- `preparation_summary.json`, conflict reports, assembly manifests, and
  `bad_accessions.txt`;
- final staged assemblies under `processed/assemblies/`;
- the nine DBGWAS feature FASTAs under `data/dbgwas/p0.05/per_species/`;
- the command/state/provenance records under each assembly work directory;
- final model checkpoints, prediction tables, RF feature tables, and their run
  names.

These are reproducible intermediates and may be regenerated when needed:

- `work/blast/dbs/` and `work/blast/hits/`;
- per-species sequence CSVs when `full_sequence_dataset.csv` is retained;
- DNABERT prediction CSVs when the checkpoint and sequence dataset are
  retained.

These are optional pilot-only artifacts rather than production inputs:

- `processed/dbgwas_pilot/*.strains.tsv`;
- `processed/pilot_assemblies_manifest.csv` and its JSON companion;
- `processed/training_metadata_balanced_pilot.csv`;
- `work/pilot_metadata/` and the `balanced_pilot_smoke` outputs.

Do not delete pilot artifacts until the run has been reviewed. After that,
they can be archived or removed as a group. The downloaded model under
`work/models/` is also a cache: keep it if the smoke test must be reproduced,
otherwise replace it with the approved fine-tuned checkpoint before production
inference. Empty per-species CSV/RF files are harmless bookkeeping outputs;
they indicate that the current metadata subset contains no accessions for that
species and should not be used for training.

## 5. Prepare Metadata

Run:

```bash
pixi run python data_pipeline/scripts/prepare_camda_inputs.py
```

The script reads the two raw CAMDA CSV files and writes:

| Output | Meaning |
| --- | --- |
| `training_metadata.csv` | One deterministic row per clean training accession. |
| `testing_metadata.csv` | Testing metadata with antibiotic derived from species. |
| `training_conflicts.csv` | Source rows excluded because one accession has conflicting records. |
| `training_assembly_manifest.csv` | Expected training assembly paths and presence status. |
| `testing_assembly_manifest.csv` | Expected testing assembly paths and presence status. |
| `dbgwas_manifest.csv` | Presence check for the nine DBGWAS feature FASTAs. |
| `preparation_summary.json` | Counts, label policy, and train/test overlap. |

The raw files are not modified. The clean training metadata is deliberately
smaller than the raw table because duplicate accessions with contradictory
phenotypes or measurements are excluded instead of being silently overwritten.
On this checkout that is 6,144 raw training rows → 5,367 clean accessions.

These two CSVs are the **canonical isolate roster**. Later steps need a
shared answer to “which accessions, which species, which antibiotic, and
(for training) which label.” An assembly FASTA or SRA run does not carry
that. The accession is only an SRA key; species and phenotype live here.

What each processed file keeps:

- **Training:** one row per clean accession, with `genus`, `species`,
  `accession`, `phenotype`, `antibiotic`, plus the extra CAMDA columns.
- **Testing:** one row per accession and no phenotype (the challenge
  labels are hidden). `antibiotic` is filled from the species map because
  the raw test file does not have it.

What they enable:

| Need | Why the metadata file is required |
| --- | --- |
| SRA / assembly staging | `stage_camda_assemblies.py` looks up the accession and whether it is train or test, then writes `<accession>.fasta`. |
| Balanced subsets | `select_balanced_subset.py` samples by species × binary phenotype from `training_metadata.csv`. |
| DBGWAS strains | `make_dbgwas_strains.py` needs ID, binary phenotype, and the path to that isolate’s FASTA. |
| BLAST / sequence CSV | `pipeline.py` maps accession → species, antibiotic, and (if `--train`) phenotype onto every window row. |
| DNABERT / RF | Those tools inherit species and labels from that sequence table, not from the genome file. |

Why both files exist:

- `training_metadata.csv` is the labeled set used to train and to evaluate
  with ground truth.
- `testing_metadata.csv` is the unlabeled CAMDA holdout: same identity and
  species/antibiotic, no ground truth. Official `infer_rf.py` scoring is
  meant to fill a testing template for those accessions.

They do not contain sequence. Assemblies, DBGWAS FASTAs, and BLAST windows
are separate artifacts keyed off `accession`.

Inspect the result:

```bash
cat data/camda_amr_2025/processed/preparation_summary.json
head data/camda_amr_2025/processed/training_metadata.csv
head data/camda_amr_2025/processed/testing_metadata.csv
```

Do not continue if the required species or antibiotic validation fails.

## 6. Obtain and Stage Assemblies

AMR_v2 expects one final FASTA file per accession:

```text
data/camda_amr_2025/processed/assemblies/train/<accession>.fasta
```

There are two ways to create it.

### 6.1 Prebuilt assemblies

If an official assembly directory is available:

```bash
pixi run python data_pipeline/scripts/stage_camda_assemblies.py \
  stage-prebuilt \
  --dataset training \
  --source-dir /path/to/prebuilt/assemblies \
  --resume
```

The utility accepts `.fa`, `.fna`, and `.fasta` files. It validates each file,
copies it to the accession-based filename, and records a checksum and basic
FASTA statistics.

### 6.2 SRA reads to assembly

The Pixi environment must contain:

```text
prefetch
fasterq-dump
fastp
spades.py
```

Check them:

```bash
pixi run python data_pipeline/scripts/stage_camda_assemblies.py check-tools
```

Run one accession first:

```bash
pixi run python data_pipeline/scripts/stage_camda_assemblies.py \
  stage-sra \
  --dataset training \
  --accession DRR124693 \
  --threads 4 \
  --memory 16 \
  --resume
```

The adapter performs these operations:

1. `prefetch` downloads the SRA object.
2. `fasterq-dump --split-files` extracts FASTQ reads.
3. It detects paired-end or single-end reads.
4. `fastp` trims reads and writes HTML/JSON reports.
5. `spades.py --isolate` assembles the trimmed reads.
6. The resulting `contigs.fasta` is validated and staged as
   `<accession>.fasta`.

#### Where `prefetch` gets the data

`prefetch` is the NCBI SRA Toolkit client. The command is:

```text
prefetch <accession> -O <work-dir>/sra
```

`<accession>` is an SRA run ID (`SRR…`, `ERR…`, `DRR…`). The tool asks
NCBI’s Sequence Read Archive for that run and writes a binary `.sra`
archive. `fasterq-dump` then expands that archive into FASTQ.

This is not a webpage download. There is a human catalog page for each
run, for example
[https://www.ncbi.nlm.nih.gov/sra/SRR1187804](https://www.ncbi.nlm.nih.gov/sra/SRR1187804),
but that page is metadata. The sequence bytes come from NCBI’s SRA data
service (HTTPS, or a configured mirror such as ENA). The accession in
the CAMDA CSV is the catalog key, not a URL.

The 3+3 POC driver (`stage_balanced_subset.py`) overlaps **prefetch-next**:
while accession *N* is in fastp/SPAdes, it starts `prefetch` for accession
*N+1*. That only hides download latency. It does not change the source.
See [AMR_V2_POC_3_3.md](AMR_V2_POC_3_3.md) section 2 for how that driver
was run.

#### SRA read-member numbering

`fasterq-dump` names files according to the read-member number stored in the
SRA record. The suffix is not guaranteed to mean "biological mate 1" or
"biological mate 2". An Illumina run can contain two biological reads plus a
short technical/index read used for sample identification or demultiplexing.
In that layout, the files can be:

```text
<accession>_1.fastq    biological read 1, usually about 84 bp
<accession>_2.fastq    technical/index read, often about 7-8 bp
<accession>_3.fastq    biological read 2, usually about 76 bp
```

This is the layout observed for the cached `ERR0293xx` accessions. Most of
them have no non-empty `_2.fastq` because the technical read is empty in the
SRA table; `ERR029367` and `ERR029380` retain only a small `_2.fastq` subset.
The `_1` and `_3` files have matching spot identifiers and are the useful
paired biological reads. The technical/index read is not input to genome
assembly.

The staging script recognizes both conventional `<accession>_1.fastq` plus
`<accession>_2.fastq` pairs and confirmed SRA biological
`<accession>_1.fastq` plus `<accession>_3.fastq` pairs. For the latter it
records `paired_sra_1_3` in `state.json`, passes only `_1` and `_3` to
`fastp`, and supplies the resulting pair to SPAdes. The technical `_2` file
is retained for provenance but is not assembled. Do not rename `_3` blindly:
the selection is appropriate only after confirming the read roles from the
SRA metadata. To inspect or reproduce this case:

```bash
pixi run fasterq-dump \
  data/camda_amr_2025/work/assemblies/training/ERR029343/sra/ERR029343 \
  --split-files --include-technical \
  -O data/camda_amr_2025/work/retry_ERR029343
```

`--include-technical` is useful for diagnosis because it exposes the index
read when it exists. It does not turn an empty SRA read member into usable
sequence. The adapter ignores that technical read after confirming the
biological `_1` and `_3` pair; `--skip-technical` can also be used for future
manual extractions once those roles have been confirmed.

For one accession, inspect:

```bash
find data/camda_amr_2025/work/assemblies/training/DRR124693 \
  -maxdepth 4 -type f -printf '%p\t%s bytes\n' | sort
cat data/camda_amr_2025/work/assemblies/training/DRR124693/state.json
ls -lh data/camda_amr_2025/processed/assemblies/train/DRR124693.fasta
```

The `state.json` file records states such as `planned`, `downloaded`,
`fastq_ready`, `trimmed`, `failed`, and `staged`. The work directory also
contains command logs and provenance. `--resume` reuses valid completed
assemblies and intermediate data where possible.

Monitor a pilot in another terminal:

```bash
pgrep -af 'stage_camda|prefetch|fasterq-dump|fastp|spades.py'
```

Do not start the full dataset until one accession per species has succeeded.
Assembly can be slow and consumes substantial temporary disk space.

## 7. Create DBGWAS Inputs

DBGWAS is a separate bacterial genome-wide association tool. It uses training
assemblies and phenotype labels to identify DNA features associated with
resistance or susceptibility.

After the required training assemblies exist:

```bash
pixi run python data_pipeline/scripts/make_dbgwas_strains.py
```

For each species, this should create a tab-separated file containing:

```text
ID    Phenotype    Path
```

Example:

```text
DRR124693    1    /absolute/path/to/DRR124693.fasta
```

The script refuses to create incomplete manifests when assemblies are missing
or a species does not contain both binary classes.

DBGWAS itself must then be run separately. Convert its output with
`data_pipeline/scripts/extract_dbgwas_sig_sequences.py` (prefers unitigs
with p ≤ 0.05; if that set is empty it keeps the top `--top-n` by
p-value). The files BLAST expects are:

```text
data/dbgwas/p0.05/per_species/<species>_sig_sequences.fasta
```

On this checkout all nine of those files exist. They are **contig
placeholders** (older pilot dumps plus three filled during the 3+3 POC),
not GWAS-significant unitigs. `fill_missing_dbgwas_fastas.py` tries a
real `DBGWAS` binary and falls back to
`generate_pilot_dbgwas_features.py`. Real DBGWAS was not runnable here
(no binary, R/`bugwas`, or container). Treat a real feature FASTA as
still missing for science even though the path is present for tooling.

The current repository does not automatically download official CAMDA
feature files. Treat genuine unitig FASTAs as a prerequisite for any
BLAST run whose scores will be interpreted.

## 8. Generate BLAST and Sequence Datasets

The main script is:

```bash
data_pipeline/pipeline.py
```

It requires:

- a metadata CSV;
- an assembly directory containing `<accession>.fasta` files;
- the species-specific DBGWAS feature FASTAs;
- BLAST tools, including `makeblastdb` and `blastn`.

The pipeline first loads metadata and creates mappings from accession to
species, antibiotic, and (for training) phenotype. It then performs the
following work for each accession:

1. Builds a BLAST nucleotide database from the assembly.
2. Aligns DBGWAS feature sequences to the assembly.
3. Writes a ten-column hit table:

   ```text
   qseqid sseqid pident length qstart qend sstart send sstrand bitscore
   ```

4. Extracts a sequence window around each hit, normally using `seq_length=1000`.
5. Writes one row for each extracted feature sequence.
6. Adds a fallback sequence row if no feature sequence was found.

The generated full sequence dataset has these columns:

| Column | Meaning |
| --- | --- |
| `sequence` | DNA window passed to DNABERT. |
| `accession` | Isolate/SRA accession. |
| `query_id` | DBGWAS feature identifier. |
| `hit_count` | Number of hits for that feature in the accession. |
| `species` | Numeric species ID used internally. |
| `antibiotic` | Numeric antibiotic ID used internally. |
| `phenotype` | Numeric training label; omitted for test data. |

Typical training invocation:

```bash
pixi run python data_pipeline/pipeline.py \
  --assemblies_dir data/camda_amr_2025/processed/assemblies/train \
  --metadata_path data/camda_amr_2025/processed/training_metadata.csv \
  --model_type sequence_based \
  --grouping per_species \
  --perspecies_dbgwas_dir data/dbgwas/p0.05/per_species \
  --output_dir data/camda_amr_2025/processed/datasets \
  --run_name camda_train \
  --train
```

Expected output includes:

```text
data/camda_amr_2025/processed/datasets/camda_train/
  sequence_based/per_species/train/full_sequence_dataset.csv
  sequence_based/per_species/train/<species>/<species>_sequence_dataset.csv
```

The `grouping` value controls how sequence datasets and models are organized:

- `per_species`: one dataset/model for each species.
- `per_antibiotic`: one dataset/model for each antibiotic.
- `full`: one combined dataset/model.

Downstream RF output is ultimately organized per species.

## 9. Train and Run DNABERT

DNABERT is a binary sequence classifier. It reads the `sequence` column from
the CSV produced by `pipeline.py`, tokenizes and truncates each DNA window,
and writes a resistant/susceptible score for that window. It does not make
the isolate-level call. The random forest later combines those scores with
BLAST hit counts.

The README names this step as `dnabert/finetune/run_finetune.sh` and
`dnabert/finetune/data/dnabert_preprocessing.py`. Those paths are stale:
the launch scripts live under `dnabert/finetune/finetune_scripts/`, and the
preprocessing script is not in this checkout. Links to the microbiome
pretrained model and a finished fine-tuned checkpoint are also `TODO` in
the README. This section documents what the code actually expects and how
to fill those gaps.

### 9.1 What the model sees

`pipeline.py` extracts a window around each DBGWAS BLAST hit, normally
`seq_length=1000`. Fine-tuning and inference both consume that window.
Training labels follow the same encoding as the rest of the pipeline:

```text
Susceptible -> 1
Intermediate -> 0
Resistant   -> 0
```

`dnabert/finetune/train.py` requires a 7-column CSV whose last column is
that integer label:

```text
sequence, accession, query_id, hit_count, species, antibiotic, phenotype
```

That is the training schema written by `pipeline.py`. The trainer reads
`sequence` from column 0 and `phenotype` from column 6. Optional extras
`--use_hit_count` and `--use_species_antibiotic` read columns 3-5; the
default launch scripts leave those off and let the random forest use them.

#### Tokenization of BLAST windows

The model never sees the raw 1000-character DNA string as text. Both
`train.py` and `inference.py` load the DNABERT-2 tokenizer
(`AutoTokenizer`, `trust_remote_code=True`) and convert each `sequence`
into integer token IDs.

This is **BPE** (byte-pair encoding), not the overlapping k-mer scheme
from DNABERT-1. The launch scripts set `--kmer -1`, which means “do not
split the DNA into 3-mers / 6-mers first.” The tokenizer’s own vocabulary
(4096 DNA pieces in `zhihan1996/DNABERT-2-117M`) is applied directly to
the A/C/G/T string.

What happens to one BLAST window:

1. `pipeline.py` has already cut ~1000 bp of contig sequence around the
   hit (`seq_length=1000`).
2. The tokenizer merges frequent nucleotide substrings into tokens.
   DNABERT-2 BPE averages about **four nucleotides per token**.
3. Special tokens are added (classification / padding). Those count
   against the length budget.
4. The result is truncated or padded to `--model_max_length 250`.
   Inference hard-codes the same limit (`MODEL_MAX_LENGTH = 250`).
5. The model forward pass receives only `input_ids` and
   `attention_mask`. By default it does **not** receive `hit_count`,
   `species`, or `antibiotic`.

The 250-token cap is why a 1000 bp window is used: 1000 ÷ 4 ≈ 250. A
longer window would be truncated and the ends of the BLAST context would
be dropped. A much shorter window would waste the context the model was
trained to use.

The tokenizer output is intermediate. It is not written to disk and is
not an RF feature. The only DNABERT values that later stages keep are
the class probabilities described below.

### 9.2 What this repository contains

Present:

- `dnabert/finetune/train.py`: Hugging Face `Trainer` that loads
  `AutoModelForSequenceClassification` with `trust_remote_code=True`.
- `dnabert/finetune/finetune_scripts/`: SLURM wrappers for one full
  model, nine per-species models, four per-antibiotic models, 3-fold
  out-of-fold training, and an experimental Nucleotide Transformer run.
- `dnabert/inference/inference.py`: sequence-level inference and RF table
  construction.
- `work/models/dnabert2-117m/`: a local copy of the public
  `zhihan1996/DNABERT-2-117M` base checkpoint (117M parameters, 12
  layers, hidden size 768, vocab 4096).
- `dnabert/finetune/data/prepare_dnabert_splits.py`: accession-level
  split used for the 3+3 POC (80/20, test = dev, optional per-accession
  cap).
- All nine paths under `data/dbgwas/p0.05/per_species/`. They are
  placeholder contig FASTAs, not official or DBGWAS-significant
  sequences. See section 7.
- A plumbing fine-tune at
  `dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/best`
  (1 CPU epoch on a capped 3+3 set). That is not a production AMR
  checkpoint. Details are in [AMR_V2_POC_3_3.md](AMR_V2_POC_3_3.md).

Missing from this checkout:

- `dnabert/finetune/pretrained_models/bacteria_model`, the
  microbiome-pretrained DNABERT-MB weights used for the CAMDA runs.
- A finished challenge-scale AMR `best` checkpoint.
- `dnabert/finetune/data/dnabert_preprocessing.py` (historical splitter;
  use `prepare_dnabert_splits.py` instead).
- Official CAMDA DBGWAS unitig FASTAs, or a runnable DBGWAS install.
- Published download URLs. The README still marks those as `TODO`.

Do not treat `work/models/dnabert2-117m` as a trained AMR model. Its
classifier head is newly initialized. Predictions from that directory
are plumbing checks only.

### 9.3 Intended checkpoint and acceptable substitutes

The production starting point was a microbiome-pretrained DNABERT-MB
directory at `dnabert/finetune/pretrained_models/bacteria_model`. The
launch scripts then copy that directory's custom Python files
(`bert_layers.py`, `bert_padding.py`, `configuration_bert.py`,
`flash_attn_triton.py`) into `<output_dir>/best` so inference can reload
the architecture.

If DNABERT-MB is still unavailable, the following substitutes are
compatible with `train.py`. They are not equivalent to the CAMDA
training setup.

| Starting weights | Use | Architecture match | Cons |
| --- | --- | --- | --- |
| DNABERT-MB / `bacteria_model` | Preferred when the files exist | Exact | Not in this repo; download link is `TODO`. |
| `work/models/dnabert2-117m` or Hugging Face `zhihan1996/DNABERT-2-117M` | Recommended public substitute | Same DNABERT-2 code, BPE, 117M size, `trust_remote_code` path | Pretrained on general genomes, not the microbiome corpus used for CAMDA. Expect some domain shift and do not claim to reproduce the published scores. |
| `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` | Experimental only; see `run_finetune_nt.sh` | Different model; `train.py` already special-cases InstaDeep tokenizers | 500M parameters, much higher GPU memory. Its tokens are 6-mers, so `model_max_length=250` means about 1500 bp, not a 1000 bp / 4 DNABERT-2 window. Not the CAMDA production model. |
| DNABERT-1, for example `zhihan1996/DNA_bert_6` | Avoid | Different tokenizer and weights | Requires `--kmer 6` and k-merized inference. `inference.py` feeds raw DNA to the tokenizer and will not match a k-mer-trained model. |

Use the Hugging Face DNABERT-2-117M weights unless you recover the
original `bacteria_model`. Keep the substitution in the run name and in
any later write-up.

### 9.4 Fill the missing pieces

Work through these steps in order. Later steps assume the earlier
artifacts exist.

`train.py` and the split helper below import `sklearn`. That package is
listed in `pixi.toml` as `scikit-learn`. If this environment was created
before that dependency was added, run `pixi install` once so the
environment picks it up.

#### Step 1. Confirm the sequence dataset exists

Fine-tuning starts from the training CSV written in section 8:

```text
data/camda_amr_2025/processed/datasets/camda_train/sequence_based/per_species/train/full_sequence_dataset.csv
```

For a single full model, that combined file is enough. For per-species
or per-antibiotic models, also keep the grouped CSVs from `pipeline.py`.
Drop rows with a null `sequence` before splitting. Confirm the seven
training columns listed in section 9.1 and that both phenotype classes
are present.

If BLAST/DBGWAS is still incomplete, stop here. Fine-tuning on a
six-species subset is only a plumbing exercise.

#### Step 2. Obtain starting weights

If `work/models/dnabert2-117m/pytorch_model.bin` already exists, point
`--model_name_or_path` at that directory. To download a fresh copy:

```bash
mkdir -p dnabert/finetune/pretrained_models
pixi run python - <<'PY'
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id="zhihan1996/DNABERT-2-117M",
    local_dir="dnabert/finetune/pretrained_models/dnabert2-117m",
)
print("downloaded DNABERT-2-117M")
PY
```

That snapshot includes the custom BERT Python files. If you later
recover `bacteria_model`, place it at
`dnabert/finetune/pretrained_models/bacteria_model` and use that path
instead.

#### Step 3. Split by accession

`dnabert/finetune/data/dnabert_preprocessing.py` is still absent. Use
`dnabert/finetune/data/prepare_dnabert_splits.py`, which is in this
checkout. It splits **by accession** so every window from one isolate
stays in one split. The 3+3 POC used 80/20 train/dev and copied dev to
`test.csv`. `--drop_dupes` removes exact duplicate `sequence` strings.
`--max_per_accession N` keeps at most N windows per isolate (the POC
used 50 so a CPU epoch was feasible).

It first tries to stratify on species × binary phenotype. If that
holdout is too small (the 3+3 set has 11 dev isolates and 18 strata),
it falls back to binary phenotype only and prints
`stratify=binary_phenotype`.

```bash
pixi run python dnabert/finetune/data/prepare_dnabert_splits.py \
  --full_dataset_path data/camda_amr_2025/processed/datasets/camda_train/sequence_based/per_species/train/full_sequence_dataset.csv \
  --out_dir dnabert/finetune/data/camda_full \
  --train_frac 0.80 \
  --drop_dupes
```

There is no `--dev_frac` on this script. Historical OOF used 90/10 with
test = dev; a larger run can still do that by editing the script or
writing a third split by hand.

Inspect the printed counts and the `*_accessions.txt` files. Every split
should contain both phenotype classes. If you later train per-species or
per-antibiotic models, filter these CSVs by the numeric `species` or
`antibiotic` column and write one directory of `train.csv` / `dev.csv` /
`test.csv` per group.

#### Step 4. Write the inference accession map

`inference.py` does not read the metadata CSV. It requires
`<metadata_dir>/train_species_to_accession.json` (or
`test_species_to_accession.json`). All nine species keys must be
present; use an empty list for a species that is absent from the
current subset.

```bash
mkdir -p work/dnabert_metadata
pixi run python - <<'PY'
import json
from pathlib import Path
import pandas as pd

species = [
    "neisseria_gonorrhoeae", "staphylococcus_aureus",
    "streptococcus_pneumoniae", "salmonella_enterica",
    "klebsiella_pneumoniae", "escherichia_coli",
    "pseudomonas_aeruginosa", "acinetobacter_baumannii",
    "campylobacter_jejuni",
]
md = pd.read_csv("data/camda_amr_2025/processed/training_metadata.csv")
key = md["genus"].str.lower() + "_" + md["species"]
grouped = md.groupby(key)["accession"].apply(lambda s: sorted(set(s))).to_dict()
mapping = {name: grouped.get(name, []) for name in species}
path = Path("work/dnabert_metadata/train_species_to_accession.json")
path.write_text(json.dumps(mapping, indent=2) + "\n")
print(path, {k: len(v) for k, v in mapping.items()})
PY
```

A pilot must use its own file, not this full-dataset map. The
six-accession smoke test used `work/pilot_metadata/`.

#### Step 5. Fine-tune

Do not run the SLURM scripts unchanged. They point at
`/gpfs/scratch/...`, activate a `conda` environment named `dna`, and
contain a hardcoded Weights & Biases login. Use the Pixi environment
from this repository and log in to W&B yourself, or set
`WANDB_MODE=offline`.

`train.py` writes the selected checkpoint to `<output_dir>/best`. It
chooses that checkpoint by macro F1 because the dev and test sets are
class-imbalanced. Historical full-model settings from
`run_finetune.sh` were 24 epochs, learning rate `3e-5`, per-device
batch size 8, gradient accumulation 4, and 4 GPUs (effective batch
128). Per-species and per-antibiotic runs used fewer epochs.

Single-GPU equivalent of that effective batch size:

```bash
export WANDB_MODE=offline
export WANDB_NAME=camda_dnabert2_117m_full
mkdir -p dnabert/finetune/finetuned_models/camda_dnabert2_117m_full

pixi run python dnabert/finetune/train.py \
  --model_name_or_path work/models/dnabert2-117m \
  --data_path dnabert/finetune/data/camda_full \
  --kmer -1 \
  --run_name camda_dnabert2_117m_full \
  --model_max_length 250 \
  --per_device_train_batch_size 8 \
  --per_device_eval_batch_size 16 \
  --gradient_accumulation_steps 16 \
  --learning_rate 3e-5 \
  --num_train_epochs 24 \
  --fp16 \
  --eval_strategy steps \
  --eval_steps 700 \
  --save_steps 700 \
  --warmup_steps 675 \
  --logging_steps 100 \
  --output_dir dnabert/finetune/finetuned_models/camda_dnabert2_117m_full \
  --overwrite_output_dir True \
  --log_level info \
  --find_unused_parameters False
```

After training, copy the DNABERT-2 runtime files into `best` so
inference can load the custom architecture:

```bash
BEST=dnabert/finetune/finetuned_models/camda_dnabert2_117m_full/best
SRC=work/models/dnabert2-117m
cp "$SRC"/bert_layers.py "$SRC"/bert_padding.py \
   "$SRC"/configuration_bert.py "$SRC"/flash_attn_triton.py \
   "$BEST"
```

This Pixi transformers build ignores `--evaluation_strategy` on the
CLI. Use `--eval_strategy`. `train.py` also used to fail on CPU with
`torch.distributed.get_rank()` (`ValueError`); that is caught now. After
a successful save, confirm `best/training_args.bin` loads: a trailing
comma once left `evaluation_strategy` as a tuple and made the
checkpoint unpickleable even though the weights had been written.

If you used `torchrun` on several GPUs, keep the same flags as the
historical scripts and reduce `--gradient_accumulation_steps` so
`batch_size * accumulation * num_gpus` stays near 128.

Optional LoRA (`--use_lora`) exists in `train.py` but was not used for
the CAMDA launch scripts. Leave it off unless you are deliberately
changing the training recipe.

#### Step 6. Run inference with the fine-tuned checkpoint

`--grouping` must match how the model was trained. `--model_path` for a
full model is the `best` directory itself. For per-species models it is
the parent directory, and the script appends `<species>/best`.

```bash
pixi run python dnabert/inference/inference.py \
  --output_format random_forest \
  --dataset_dir data/camda_amr_2025/processed/datasets/camda_train/sequence_based/per_species/train/full_sequence_dataset.csv \
  --model_path dnabert/finetune/finetuned_models/camda_dnabert2_117m_full/best \
  --grouping full \
  --run_name camda_dnabert2_117m_full \
  --sig_seqs_dir data/dbgwas/p0.05/per_species \
  --metadata_dir work/dnabert_metadata \
  --base_dir "$PWD" \
  --train True \
  --return_logits sum
```

Default values inside `inference.py` still point at the old HPC tree.
Always pass `--base_dir`, `--metadata_dir`, and `--sig_seqs_dir`.

The script writes sequence rows plus two new columns:

```text
pred_res
pred_sus
```

`pred_res` is the softmax probability of class 0 (resistant /
intermediate). `pred_sus` is class 1 (susceptible). They sum to 1 on
each row. Null sequences are dropped. This per-window file (for example
`dnabert/inference/outputs/preds/<run>/full/train/full_preds.csv`) is
**not** what the random forest trains on.

#### What the RF stage actually uses

`get_rf_dataset` collapses those window rows to **one row per
accession**. The feature list is the headers in the species DBGWAS
FASTA, not the token IDs. For each feature `Q` it writes:

```text
Q_hit_count
Q_pred_resistant
```

`Q_hit_count` is copied from the BLAST `hit_count` column. It is not a
DNABERT output; it is passed through so the forest can use alignment
abundance.

`Q_pred_resistant` is the only DNABERT-derived value in the RF table.
How windows are combined is `--return_logits`:

| `--return_logits` | `Q_pred_resistant` |
| --- | --- |
| `sum` (POC and recommended) | sum of `pred_res` over windows with that query ID |
| `average` | mean of `pred_res` |
| omitted (`None`) | count of windows whose argmax is resistant (labels flipped) |
| `all` | six stats from both `pred_res` and `pred_sus` (sum/mean/std) |

The 3+3 POC used `sum`, so each RF DNABERT column is a **sum of
resistant-class probabilities** for that feature in that isolate.
`pred_sus` is unused in that mode. Token IDs, attention masks, and the
raw `sequence` strings are discarded.

The isolate-level S/R call is made by the forest from
`*_hit_count` + `*_pred_resistant` (or a subset if
`--feature_type` is `hits` or `dnabert`). DNABERT itself never writes
that accession label.

Hugging Face caches a module named `best`. If you change the model
architecture and keep that directory name, clear the cache or inference
may reload the previous code.

### 9.5 Grouping

| Grouping | DNABERT models | Dataset consumed | RF output |
| --- | --- | --- | --- |
| `full` | 1 | `full_sequence_dataset.csv` | still one RF table per species |
| `per_species` | 9 | `<species>/<species>_sequence_dataset.csv` | one RF table per species |
| `per_antibiotic` | 4 | `<antibiotic>/<antibiotic>_sequence_dataset.csv` | one RF table per species |

Downstream RF is always per species. Train DNABERT with the same
grouping that `pipeline.py` used to write the sequence CSVs.

### 9.6 Smoke test versus real inference

The six-accession smoke test used the untuned base checkpoint on
purpose:

```bash
pixi run python dnabert/inference/inference.py \
  --output_format random_forest \
  --dataset_dir data/camda_amr_2025/processed/pilot_outputs/balanced_pilot/sequence_based/per_species/train/full_sequence_dataset.csv \
  --model_path work/models/dnabert2-117m \
  --grouping full \
  --run_name balanced_pilot_smoke \
  --sig_seqs_dir work/pilot_metadata/per_species_sig_sequences \
  --metadata_dir work/pilot_metadata \
  --base_dir "$PWD" \
  --train True \
  --return_logits sum
```

That command checks file formats and data flow. Its classifier head is
random. Use an approved fine-tuned `best` checkpoint for any result that
will be interpreted as AMR.

## 10. Build RF Features and Predictions

`get_rf_dataset` uses the headers in the species DBGWAS FASTA as the canonical
feature list. How window-level `pred_res` is aggregated into
`<feature>_pred_resistant` is described in section 9.4. For each accession
it creates columns such as:

```text
<feature>_hit_count
<feature>_pred_resistant
```

`*_hit_count` is BLAST abundance. `*_pred_resistant` is the only DNABERT
output the forest sees (with `--return_logits sum`, a sum of resistant
probabilities). Token IDs and per-window `sequence` / `pred_sus` rows are
not RF inputs.

Absent features are filled with zero so training and testing have the same
column set. Training RF datasets also contain:

```text
ground_truth_phenotype
```

The RF code supports these feature choices:

- `hits`: BLAST hit counts only.
- `dnabert`: DNABERT predictions only.
- `both`: hit counts and DNABERT predictions.

RF models are normally trained per species. The resulting model files are
usually `.joblib` files. For final testing, `rf/infer_rf.py` loads the per-
species models, predicts each accession, maps numeric predictions back to
`Susceptible` or `Resistant`, and writes a submission-style CSV.

Local runs must override HPC defaults and match the table layout
`train_rf.py` expects:

- `--base_dir "$PWD"` (default is `/gpfs/scratch/jvaska/CAMDA_AMR/AMR_v2`).
- `--flip_phenotype` pointing at empty per-species lists, or the HPC
  mismatches directory will be opened and fail. The 3+3 POC used
  `work/rf_poc/no_flip/`.
- `--train_on all` reads
  `dataset_dir/fold_*/<dnabert_grouping>/train/<species>/`. Inference
  writes a flat `.../full/train/<species>/` tree; symlink that under a
  fake `fold_0` (see [AMR_V2_POC_3_3.md](AMR_V2_POC_3_3.md) section 6)
  or build real folds.
- `infer_rf.py` must drop `ground_truth_phenotype` before `predict` if
  that column is still on the frame. That drop is in the current script.
- Python 3.14 `argparse` rejects a raw `%` in help text; the
  `--interleave` help string is already escaped.

Do not pass `--eval` or `--oof_stack` until real held-out folds exist.
In-sample accuracy on a handful of isolates is not a metric.

## 11. Beginner Checklists

### Before metadata preparation

- Raw training and testing CSVs exist.
- Raw files have not been edited manually.
- Pixi environment is available.

### After metadata preparation

- `training_metadata.csv` exists.
- `testing_metadata.csv` exists.
- `training_conflicts.csv` was reviewed.
- `preparation_summary.json` was inspected.

### After assembly staging

- Every requested accession has `<accession>.fasta`.
- `state.json` says `staged`.
- FASTA files are nonempty and contain nucleotide sequences.
- Failed accessions have been investigated before scaling.

### Before DBGWAS

- All training assemblies exist.
- Nine `strains` manifests exist.
- Each manifest contains both phenotype classes.

### Before `pipeline.py`

- All nine `*_sig_sequences.fasta` files exist. A file on that path can
  still be a contig placeholder; confirm origin before treating scores
  as GWAS features.
- BLAST tools are installed.
- Assembly paths match accession names.
- The metadata and assembly datasets refer to the same accessions.

### Before DNABERT

- Sequence CSVs contain nonempty `sequence` values and seven training
  columns.
- The train/dev/test split is defined by accession, not by sequence row.
- Starting weights exist: recovered `bacteria_model`, or the public
  DNABERT-2-117M substitute documented in section 9.3.
- `train.csv`, `dev.csv`, and `test.csv` exist for the chosen grouping
  (`prepare_dnabert_splits.py` writes these).
- `train_species_to_accession.json` covers all nine species.
- For real inference, `<output_dir>/best` is a fine-tuned checkpoint,
  not `work/models/dnabert2-117m`.

### Before RF

- DNABERT prediction CSVs exist.
- RF datasets contain the expected feature columns.
- Training includes `ground_truth_phenotype`.
- Testing uses the same feature headers as training.
- RF input is the output of `dnabert/inference/inference.py`, not the raw
  sequence dataset from `pipeline.py`.
- Local `train_rf.py` has `--base_dir` set and a `fold_*` table layout
  (section 10).
- Each species has enough accessions and both phenotype classes for the
  intended model; a six-accession pilot is suitable for plumbing checks, not
  model assessment.

## 12. Common Problems

| Symptom | Likely cause | What to inspect |
| --- | --- | --- |
| Missing assembly | SRA download or staging failed. | `state.json`, SRA log, final assembly directory. |
| Incomplete FASTQ | Single/paired layout was ambiguous or incomplete; an SRA technical read may also place the biological mate in `_3` instead of `_2`. | `fastq/` filenames, SRA metadata, and `fasterq-dump.log`. |
| SPAdes exit failure | Insufficient memory, bad reads, or incomplete output. | `spades/spades.log`, memory, disk space. |
| Missing DBGWAS FASTA | DBGWAS has not been run or output was not converted. | `data/dbgwas/p0.05/per_species/`. |
| No BLAST hits | Feature/assembly mismatch or genuinely absent feature. | `<accession>_hits.tsv`, feature FASTA, assembly. |
| Null sequence rows | No usable flanking sequence was extracted. | Sequence dataset and pipeline log. |
| `dnabert_preprocessing.py` missing | That script is not in this checkout. | Use `prepare_dnabert_splits.py` (section 9.4, step 3). |
| `--evaluation_strategy` ignored | This transformers build uses `eval_strategy`. | Pass `--eval_strategy`. |
| `train.py` `ValueError` on `get_rank()` | CPU / non-DDP process group. | Current `train.py` catches that; pull the fix if you still see it. |
| Untuned AMR scores | Inference used `work/models/dnabert2-117m` or another base checkpoint. | Confirm `--model_path` ends in a fine-tuned `best` directory. |
| `train.py` cannot load the model | Custom DNABERT-2 files were not next to the weights, or `trust_remote_code` failed. | Check `bert_layers.py` and related files in the model directory. |
| Inference reads HPC paths | `--base_dir`, `--metadata_dir`, or `--sig_seqs_dir` was omitted. | Pass local paths; defaults still point at `/gpfs/scratch/...`. |
| Hugging Face reused an old `best` | The cache kept a previous module named `best`. | Clear the cache or rename the checkpoint directory. |
| Stratified split fails | A species/phenotype pair has only one accession. | Merge that stratum or split by species only. |
| Feature mismatch in RF | Training and testing used different feature FASTAs. | FASTA headers and RF CSV columns. |
| Old error in `state.json` | A previous retry failed, then a later retry succeeded without clearing the old field. | Trust `state: staged` plus final FASTA, but review provenance. |
| Summary replaced | Multiple runs reused the same `--manifest` path. | Use a unique manifest per run/species. |

## 13. Pilot Before Scale-Up

Start with one accession:

```bash
pixi run python data_pipeline/scripts/stage_camda_assemblies.py \
  stage-sra --dataset training --accession DRR124693 \
  --threads 4 --memory 16 --resume
```

Then run one accession per species using
`one_per_species_pilot.sh`. Keep a separate manifest or log per species so
later runs do not overwrite the previous summary.

Only after the pilot succeeds should you stage all training assemblies. Do not
run the expensive downstream stages until the required inputs are complete.

The six-accession pilot used during development produced these useful checks:

- `full_sequence_dataset.csv` contained 14,904 sequence rows from six
  accessions;
- DNABERT inference produced one prediction row per sequence;
- RF feature tables were produced for the six populated species;
- the three absent species produced empty tables, as expected.

The pilot metadata had one accession per populated species and was strongly
imbalanced: five accessions had encoded phenotype `0` and one had encoded
phenotype `1`. It therefore validates file formats and data flow only. Do not
use it to claim model accuracy or train production RF models.

The 3+3 tooling POC on this machine is documented in
[AMR_V2_POC_3_3.md](AMR_V2_POC_3_3.md) (53 staged isolates, placeholder
FASTAs, one CPU DNABERT epoch, nine RF smoke models). That run is a
plumbing check, not a scale-up.

The repository still lacks official CAMDA assembly dumps, real DBGWAS
unitig FASTAs, DNABERT-MB (`bacteria_model`), and a challenge-scale AMR
checkpoint. Section 9 documents public substitutes. Treat those gaps as
explicit stop points rather than assuming that metadata preparation or
assembly creates them.

## 14. Current Limitations

This guide describes the current repository behavior, including limitations:

- Species and antibiotic mappings are hard-coded.
- Assembly staging depends on external SRA data and heavy tools.
- DBGWAS execution and feature extraction are not fully automated here.
- Some scripts contain HPC-specific absolute paths.
- `dnabert/finetune/data/dnabert_preprocessing.py` is not in this
  checkout; use `prepare_dnabert_splits.py` as in section 9.4.
- The documented split option in `pipeline.py` is not fully implemented.
- The matrix-based generator path is incomplete.
- RF and DNABERT models are not included automatically. The local
  `work/models/dnabert2-117m` directory is an untuned public base
  checkpoint, not the CAMDA fine-tuned model.
- Using Hugging Face `zhihan1996/DNABERT-2-117M` in place of DNABERT-MB
  is architecturally compatible but is a domain-shift substitute.
- Full assembly and model workflows require substantial compute and storage.
- Challenge testing overlap and independent evaluation must be interpreted
  carefully.

Read the generated manifests and logs at every stage. A file existing on disk
is not enough by itself; confirm its schema, row count, status, and provenance.

## Note to User
The SRA-to-assembly stage (Section 6.2) reconstructs inputs this repo was written to **consume**, not to generate from scratch every time.

The DNABERT + RF tooling does not need 5,000 SRA downloads. It needs:

1. Assemblies named `<accession>.fasta`
2. Nine DBGWAS feature FASTAs
3. `pipeline.py` BLAST windows
4. A fine-tuned (or at least loadable) DNABERT checkpoint
5. RF tables from `inference.py`

The README states that directly: the expected input is CAMDA isolates **already** QCed with fastp and assembled with SPAdes `--isolate`. The original project almost certainly started from those assemblies on HPC, plus feature FASTAs they already had. Downloads for both are still `TODO`. `stage_camda_assemblies.py` is the fallback for when those files are not here.

So the months of wall time for the assembly staging are **data recovery**, not model training. Fine-tuning 117M and fitting nine small RFs is cheap compared with SPAdes.

What you can skip or shrink:

- **Not all 5,367 isolates.** The 3+3 POC already walked the tools
  ([AMR_V2_POC_3_3.md](AMR_V2_POC_3_3.md)). A 5+5 or 10+10 subset is
  the next “use the tooling” size. Full staging is only for a
  challenge-scale replay.
- **Real DBGWAS** if the goal is only to exercise BLAST → DNABERT → RF.
  All nine `*_sig_sequences.fasta` files on disk are contig
  placeholders. That is scientifically weaker, but it unblocks the
  tools.
- **SRA entirely** if you can get a prebuilt assembly directory. Then section 6 is `stage-prebuilt`, not 160 SPAdes jobs.
- **DNABERT** if you only want RF. The code can train on hit counts alone (`--feature_type hits`).

What you cannot skip if you want *this* pipeline’s sequence model: some assemblies, some feature FASTAs, and BLAST windows. There is no script that trains DNABERT on raw SRA or on metadata CSVs.

The honest read: AMR_v2 is an inference/training stack sitting on a heavy CAMDA data product. This checkout is missing that product, so it feels like the repo “requires” assembly at scale. It requires those **files**. Building them yourself is optional, expensive, and only necessary because the official dumps were never added.

## 15. Optional: DNABERT embeddings (not implemented)

The public checkpoint used here (`zhihan1996/DNABERT-2-117M`, cached at
`work/models/dnabert2-117m`) can produce sequence embeddings. AMR_v2 does
not do that. This section is scope only: what the weights allow, and what
would have to change if someone used them that way.

### What the current path does

`train.py` and `inference.py` load
`AutoModelForSequenceClassification`. The encoder still computes 768-
dimensional hidden states internally, then a linear head maps them to two
logits. Softmax yields `pred_res` and `pred_sus`. Those scalars are
aggregated into `*_pred_resistant` for the random forest (section 9.4).
Hidden states, CLS vectors, and token embeddings are discarded.

### What the same checkpoint can do

The model card registers both:

- `AutoModel` → `BertModel` (encoder only)
- `AutoModelForSequenceClassification` → encoder plus a classification head

Loading `AutoModel` (or reading `hidden_states` from the classification
model) would give a 768-d vector per BLAST window after the same BPE
tokenization (`kmer=-1`, max length 250). Typical choices:

- the first-token / CLS vector, or
- the mean of last-layer token states, masked so padding is ignored

That is an embedding of the window, not a phenotype probability.

### What would still have to be designed

Embeddings are **per window**. The RF tables in this repo are **per
accession**, with one scalar per DBGWAS feature. An embedding path would
need a new aggregation, for example:

- mean of window vectors for each feature, or
- mean of all windows in the isolate, or
- concatenation / a different isolate-level model

The forest (or a replacement model) would then see 768-d blocks instead
of `*_pred_resistant`. `*_hit_count` could still be attached. That is a
different feature schema, different training code, and it would not
reproduce the CAMDA DNABERT→RF recipe.

### Scope

- **In scope for this guide:** the 117M weights are an encoder; embeddings
  are possible in principle.
- **Out of scope for this repository:** there is no embedding export in
  `inference.py`, no 768-d RF table builder, and no script that trains on
  pooled hidden states.
- **Not a substitute for the current POC:** swapping in embeddings would
  be a new experiment, not a drop-in flag on the existing `--return_logits`
  path.

If that experiment is worth doing, keep it in a separate run name and
leave the classification RF tables intact.