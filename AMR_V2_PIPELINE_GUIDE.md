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
| `<output_dir>/<run_name>/` | Sequence datasets produced by `pipeline.py`; this repository commonly uses `data/camda_amr_2025/processed/pilot_outputs/`. |
| `dnabert/inference/outputs/` | DNABERT predictions and RF input datasets. |
| `rf/models/` and `rf/predictions/` | Trained models and final predictions. |

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

DBGWAS itself must then be run separately. Its output must eventually be
converted into these nine files:

```text
data/dbgwas/p0.05/per_species/<species>_sig_sequences.fasta
```

The current repository does not automatically download DBGWAS features or
fully extract them from every DBGWAS output version. Treat the feature FASTAs
as a prerequisite for the BLAST stage.

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

DNABERT reads the `sequence` column from the generated sequence dataset. It
tokenizes and truncates each DNA window, then writes sequence-level
predictions. It does not directly produce one final prediction per accession.

The exact preprocessing script should be checked in the local repository before
running because some historical paths are referenced by older OOF scripts. In
general, training preprocessing creates `train.csv`, `dev.csv`, and `test.csv`
from accession lists, then DNABERT finetuning writes a model directory with a
`best` checkpoint.

DNABERT inference uses a command shaped like:

```bash
pixi run python dnabert/inference/inference.py \
  --dataset_dir <sequence-dataset-directory> \
  --model_path <finetuned-model> \
  --grouping per_species \
  --output_format random_forest \
  --run_name camda_train \
  --train True \
  --sig_seqs_dir data/dbgwas/p0.05/per_species
```

The inference output contains sequence rows plus:

```text
pred_res
pred_sus
```

Rows with null sequences are dropped. The next function aggregates these
sequence-level predictions by accession and DBGWAS feature.

For a pilot subset, inference also needs a JSON mapping from species to the
accessions included in that subset. The standard CAMDA workflow obtains this
from the prepared metadata; a pilot must use a separately named derived file
so that it does not accidentally use the full-dataset mapping. The mapping
must cover all nine species, using an empty list for species absent from the
pilot.

For example, the session's six-accession smoke test used:

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

This command is a pipeline smoke test only when `work/models/dnabert2-117m`
is the untuned base checkpoint. Its classifier head is newly initialized, so
the resulting predictions must not be interpreted as AMR results. Use an
approved fine-tuned `best` checkpoint for real inference.

## 10. Build RF Features and Predictions

`get_rf_dataset` uses the headers in the species DBGWAS FASTA as the canonical
feature list. For each accession it creates columns such as:

```text
<feature>_hit_count
<feature>_pred_resistant
```

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

- All nine `*_sig_sequences.fasta` files exist.
- BLAST tools are installed.
- Assembly paths match accession names.
- The metadata and assembly datasets refer to the same accessions.

### Before DNABERT

- Sequence CSVs contain nonempty `sequence` values.
- The training/dev/test split is defined by accession.
- The expected model checkpoint exists.

### Before RF

- DNABERT prediction CSVs exist.
- RF datasets contain the expected feature columns.
- Training includes `ground_truth_phenotype`.
- Testing uses the same feature headers as training.
- RF input is the output of `dnabert/inference/inference.py`, not the raw
  sequence dataset from `pipeline.py`.
- Each species has enough accessions and both phenotype classes for the
  intended model; a six-accession pilot is suitable for plumbing checks, not
  model assessment.

## 12. Common Problems

| Symptom | Likely cause | What to inspect |
| --- | --- | --- |
| Missing assembly | SRA download or staging failed. | `state.json`, SRA log, final assembly directory. |
| Incomplete FASTQ | Single/paired layout was ambiguous or incomplete. | `fastq/` filenames and `fasterq-dump.log`. |
| SPAdes exit failure | Insufficient memory, bad reads, or incomplete output. | `spades/spades.log`, memory, disk space. |
| Missing DBGWAS FASTA | DBGWAS has not been run or output was not converted. | `data/dbgwas/p0.05/per_species/`. |
| No BLAST hits | Feature/assembly mismatch or genuinely absent feature. | `<accession>_hits.tsv`, feature FASTA, assembly. |
| Null sequence rows | No usable flanking sequence was extracted. | Sequence dataset and pipeline log. |
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

The repository currently requires external or unavailable inputs for parts of
the workflow, especially the DBGWAS feature FASTAs and pretrained model
checkpoints. Treat those as explicit stop points rather than assuming that
metadata preparation or assembly automatically creates them.

## 14. Current Limitations

This guide describes the current repository behavior, including limitations:

- Species and antibiotic mappings are hard-coded.
- Assembly staging depends on external SRA data and heavy tools.
- DBGWAS execution and feature extraction are not fully automated here.
- Some scripts contain HPC-specific absolute paths.
- The documented split option in `pipeline.py` is not fully implemented.
- The matrix-based generator path is incomplete.
- RF and DNABERT models are not included automatically.
- Full assembly and model workflows require substantial compute and storage.
- Challenge testing overlap and independent evaluation must be interpreted
  carefully.

Read the generated manifests and logs at every stage. A file existing on disk
is not enough by itself; confirm its schema, row count, status, and provenance.