# AMR_v2 3+3 Tooling POC

This note records the class-balanced **3 susceptible + 3 resistant** proof of
concept that was run on this machine. The goal was to walk the real AMR_v2
tooling — subset selection, SRA staging, DBGWAS inputs, BLAST windows,
DNABERT-2-117M, and per-species random forests — not to produce a usable AMR
model.

Do not treat DNABERT or RF scores from this run as antimicrobial-resistance
performance. Features are contig placeholders, *S. pneumoniae* is short one
resistant isolate, DNABERT was fine-tuned for one CPU epoch on a capped
subset, and the forests were fit on 5–6 isolates with hundreds to thousands
of features.

The full pipeline language lives in [AMR_V2_PIPELINE_GUIDE.md](AMR_V2_PIPELINE_GUIDE.md).
Commands below were run from the repository root with Pixi unless noted.

Binary labels match the pipeline:

- `Susceptible` → class `S` / phenotype `1`
- `Intermediate` and `Resistant` → class `R` / phenotype `0`

---

## What this POC produced

| Stage | Result |
| --- | --- |
| Target | 9 species × 3 S × 3 R = 54 isolates |
| Staged assemblies | **53** (`streptococcus_pneumoniae` is 3 S + 2 R) |
| DBGWAS | Real DBGWAS was not runnable; all nine `*_sig_sequences.fasta` files are contig placeholders |
| BLAST + sequence dataset | 53 hit tables; `full_sequence_dataset.csv` has **426,394** rows |
| DNABERT | One full `DNABERT-2-117M` model, 1 epoch, CPU, 50 windows/accession |
| RF | Nine per-species `.joblib` models; isolate-level S/R predictions for all 53 accessions |

---

## End-to-end path

```text
CAMDA training CSV
        |
        |  prepare_camda_inputs.py          (already done)
        v
training_metadata.csv (5,367 accessions)
        |
        |  select_balanced_subset.py
        v
balanced_3_3/training_metadata.csv (54 + later backups)
        |
        |  stage_balanced_subset.py         (prefetch + SPAdes)
        v
assemblies/train/<accession>.fasta          (53 present)
        |
        |  filter to present FASTAs
        v
training_metadata_present.csv               (53 rows)
        |
        |  make_dbgwas_strains.py
        |  fill_missing_dbgwas_fastas.py    (real DBGWAS failed; placeholders)
        v
*.strains.tsv  +  data/dbgwas/p0.05/per_species/*_sig_sequences.fasta
        |
        |  pipeline.py                      (makeblastdb + blastn + windows)
        v
full_sequence_dataset.csv                   (426,394 rows)
        |
        |  prepare_dnabert_splits.py
        |  train.py                         (DNABERT-2-117M, 1 epoch)
        |  inference.py                     (capped_all.csv)
        v
per-species RF tables
        |
        |  train_rf.py  +  infer_rf.py
        v
nine *.joblib models  +  isolate-level predictions
```

---

## 0. Prerequisite: clean CAMDA metadata

This was already on disk before the 3+3 work. It is the source of every later
accession list.

**Script:** `data_pipeline/scripts/prepare_camda_inputs.py`

**Inputs**

- `data/camda_amr_2025/training_dataset.csv`
- `data/camda_amr_2025/testing_dataset.csv`

**Outputs**

- `data/camda_amr_2025/processed/training_metadata.csv` (5,367 unique accessions)
- `data/camda_amr_2025/processed/testing_metadata.csv`
- `data/camda_amr_2025/processed/training_conflicts.csv`
- `data/camda_amr_2025/processed/dbgwas_manifest.csv` (inventory only; `missing` is valid here)
- `data/camda_amr_2025/processed/preparation_summary.json`

```bash
pixi run python data_pipeline/scripts/prepare_camda_inputs.py
```

---

## 1. Select the 3+3 toy dataset

**Script:** `data_pipeline/scripts/select_balanced_subset.py`

Per species, already-staged assemblies are kept first (capped at the class
quota), then each binary class is filled with a seeded shuffle. A second
reserve of three backups per class is written so failed SRA jobs can be
replaced without re-sampling the whole set.

**Inputs**

- `data/camda_amr_2025/processed/training_metadata.csv`
- `data/camda_amr_2025/processed/assemblies/train/*.fasta` (reuse already-staged isolates)
- `data/camda_amr_2025/processed/bad_accessions.txt`

**Outputs** (all under `data/camda_amr_2025/processed/balanced_3_3/`)

| File | Role |
| --- | --- |
| `training_metadata.csv` | Selected rows (starts at 54; backups are appended on replace) |
| `accessions.txt` | Selected accession list |
| `backup_pool.csv` | Reserve replacements, with `binary_class` and `species_key` |
| `backup_accessions.txt` | Backup accession list |
| `summary.json` | Quota, seed, reused vs remaining jobs, per-class IDs |

```bash
pixi run python data_pipeline/scripts/select_balanced_subset.py \
  --metadata data/camda_amr_2025/processed/training_metadata.csv \
  --assemblies-dir data/camda_amr_2025/processed/assemblies/train \
  --bad-accession-list data/camda_amr_2025/processed/bad_accessions.txt \
  --output-dir data/camda_amr_2025/processed/balanced_3_3 \
  --quota 3 \
  --backup-quota 3 \
  --seed 42
```

This run reported **54 selected**, **54 backups**, **13 reused** staged
assemblies, and **41 remaining** SRA jobs. *S. aureus* reused all six
`ERR0293xx` assemblies from earlier staging.

Skip list at selection time (confirmed unusable):

```text
DRR170693
DRR170697
DRR170701
ERR016635
```

`ERR1218638` (*E. coli* S) was on an earlier skip list, later removed, and
was staged successfully. It is in the present 3+3 set.

---

## 2. Stage remaining assemblies

**Script:** `data_pipeline/scripts/stage_balanced_subset.py`

The driver calls `data_pipeline/scripts/stage_camda_assemblies.py stage-sra`
one accession at a time, with **prefetch-next** overlap. While the current
accession is in fastp/SPAdes, `prefetch` starts downloading the next
accession’s `.sra` archive from NCBI’s Sequence Read Archive (not from a
catalog webpage). Where those bytes come from, and how `prefetch` differs
from the human SRA page, is in [AMR_V2_PIPELINE_GUIDE.md](AMR_V2_PIPELINE_GUIDE.md)
section 6.2.

On failure it pulls the next same-species, same-class row from
`backup_pool.csv` and appends it to `training_metadata.csv`. Resume is
safe: already-present FASTAs and known `bad_accession` work states are
skipped.

**Inputs**

- `data/camda_amr_2025/processed/balanced_3_3/training_metadata.csv`
- `data/camda_amr_2025/processed/balanced_3_3/backup_pool.csv`
- `data/camda_amr_2025/processed/bad_accessions.txt`
- `data_pipeline/scripts/stage_camda_assemblies.py`

**Outputs**

- `data/camda_amr_2025/processed/assemblies/train/<accession>.fasta`
- `data/camda_amr_2025/work/assemblies/training/<accession>/` (SRA/fastp/SPAdes work; some later deleted)
- `data/camda_amr_2025/processed/balanced_3_3/logs/<accession>.stage.log`
- `data/camda_amr_2025/processed/balanced_3_3/stage_driver.json`

```bash
pixi run python data_pipeline/scripts/stage_balanced_subset.py \
  --subset-dir data/camda_amr_2025/processed/balanced_3_3 \
  --threads 4 \
  --memory 16
```

Final driver summary: **45 attempted**, **40 present**, **5 failed**,
`still_missing: []`. Combined with the 13 already-staged isolates, that is
the 53-assembly POC set.

### 2.1 The pneumococcus shortfall (the failed 3+3)

Eight species reached 3 S + 3 R. *S. pneumoniae* finished **3 S + 2 R**.

Original pneumococcus draw from `summary.json`:

| Class | Selected | Backups |
| --- | --- | --- |
| S | `ERR016736`, `SRR4417290`, `ERR124308` | `SRR4417353`, `ERR124227`, `ERR129158` |
| R | `SRR4417363`, `ERR016650`, `ERR016769` | `ERR016777`, `ERR019734`, `SRR2072233` |

`ERR016650`, `ERR016736`, and `ERR016769` failed as `ambiguous_fastq_layout`
(incomplete or unusable FASTQ layout; not the `_1`/`_3` biological-pair case
that was fixed for `ERR0293xx`). Their backups `ERR016777` and `ERR019734`
failed the same way. Successful replacements were `ERR124227` (S) and
`SRR2072233` (R). No further R backup remained for `ERR016769`.

Those five IDs were dropped when the present table was built:

```text
ERR016650
ERR016736
ERR016769
ERR016777
ERR019734
```

They were added under “Possible bad” in
`data/camda_amr_2025/processed/bad_accessions.txt`. The staging driver does
**not** append `bad_accession` IDs to that file automatically; it only
records them in `stage_driver.json` and the per-accession work state.

### 2.2 Other staging problems

- **WSL disconnect / idle shutdown** killed the long driver. Re-running the
  same command resumes.
- **Disk full** left `ERR124227` with a 0-byte `state.json`. The state file
  was repaired, unused work dirs were deleted (resume dirs and failed/skip
  work dirs were kept), and the driver was resumed. Compact the Windows
  `ext4.vhdx` from Windows after `wsl --shutdown` if C: is still tight.
- **`sklearn` was missing** from Pixi until it was added for later split/RF
  steps (`scikit-learn` in `pixi.toml`).

### 2.3 Present-only metadata

`make_dbgwas_strains.py` refuses a species table if any listed assembly is
missing, so the 58-row metadata (original 54 plus appended backups) was
filtered to FASTAs that actually exist.

**Inputs**

- `data/camda_amr_2025/processed/balanced_3_3/training_metadata.csv`
- `data/camda_amr_2025/processed/assemblies/train/<accession>.fasta`

**Outputs**

- `data/camda_amr_2025/processed/balanced_3_3/training_metadata_present.csv` (53 rows)
- `data/camda_amr_2025/processed/balanced_3_3/present_summary.json`

Present class counts:

| Species | S | R |
| --- | --- | --- |
| `acinetobacter_baumannii` | 3 | 3 |
| `campylobacter_jejuni` | 3 | 3 |
| `escherichia_coli` | 3 | 3 |
| `klebsiella_pneumoniae` | 3 | 3 |
| `neisseria_gonorrhoeae` | 3 | 3 |
| `pseudomonas_aeruginosa` | 3 | 3 |
| `salmonella_enterica` | 3 | 3 |
| `staphylococcus_aureus` | 3 | 3 |
| `streptococcus_pneumoniae` | 3 | **2** |

Use `training_metadata_present.csv` for every later POC step.

---

## 3. DBGWAS inputs and feature FASTAs

Real DBGWAS was attempted and skipped. There is no `DBGWAS` binary on PATH,
R/`bugwas` is not installed, and no container runtime is available. The
install attempt is recorded in:

- `data/camda_amr_2025/processed/balanced_3_3/dbgwas_logs/install_attempt.json`

Decision: `fallback_placeholder`. Expected external command, if a prebuilt
binary is added later:

```text
DBGWAS -strains <species>.strains.tsv -output <work_dir> -nb-cores 4
```

then

```bash
pixi run python data_pipeline/scripts/extract_dbgwas_sig_sequences.py \
  --input <dbgwas-output-dir> \
  --output data/dbgwas/p0.05/per_species/<species>_sig_sequences.fasta \
  --species <species> \
  --p-cutoff 0.05 \
  --top-n 50
```

### 3.1 Strains TSVs

**Script:** `data_pipeline/scripts/make_dbgwas_strains.py`

**Inputs**

- `data/camda_amr_2025/processed/balanced_3_3/training_metadata_present.csv`
- `data/camda_amr_2025/processed/assemblies/train/<accession>.fasta`

**Outputs** (nine files)

- `data/camda_amr_2025/processed/dbgwas/balanced_3_3/<species>.strains.tsv`

Columns: `ID`, `Phenotype` (`1`/`0`), `Path` (absolute FASTA). Eight species
have 6 rows and both classes. Pneumococcus has **5 rows** and both classes.

```bash
pixi run python data_pipeline/scripts/make_dbgwas_strains.py \
  --metadata data/camda_amr_2025/processed/balanced_3_3/training_metadata_present.csv \
  --assemblies-dir data/camda_amr_2025/processed/assemblies/train \
  --output-dir data/camda_amr_2025/processed/dbgwas/balanced_3_3
```

### 3.2 Placeholder feature FASTAs

**Scripts**

- `data_pipeline/scripts/fill_missing_dbgwas_fastas.py` (try DBGWAS, then placeholder)
- `data_pipeline/scripts/generate_pilot_dbgwas_features.py` (dump every contig as a “feature”)

The three species that were missing FASTAs at the start of this POC were
*S. aureus*, *E. coli*, and *S. pneumoniae*. The other six files were older
pilot contig dumps (dated 2025-09-17) and were left in place so BLAST had a
full set of nine queries.

**Inputs**

- `data/camda_amr_2025/processed/balanced_3_3/training_metadata.csv` (or the present table)
- `data/camda_amr_2025/processed/dbgwas/balanced_3_3/<species>.strains.tsv`
- `data/camda_amr_2025/processed/assemblies/train/<accession>.fasta`

**Outputs**

- `data/dbgwas/p0.05/per_species/<species>_sig_sequences.fasta`
- `data/camda_amr_2025/processed/balanced_3_3/dbgwas_logs/fill_summary.json`
- `data/camda_amr_2025/processed/balanced_3_3/dbgwas_logs/<species>/placeholder.log`

```bash
pixi run python data_pipeline/scripts/fill_missing_dbgwas_fastas.py \
  --metadata data/camda_amr_2025/processed/balanced_3_3/training_metadata.csv \
  --assemblies-dir data/camda_amr_2025/processed/assemblies/train \
  --strains-dir data/camda_amr_2025/processed/dbgwas/balanced_3_3 \
  --output-dir data/dbgwas/p0.05/per_species \
  --work-dir data/camda_amr_2025/work/dbgwas/balanced_3_3 \
  --log-dir data/camda_amr_2025/processed/balanced_3_3/dbgwas_logs \
  --species staphylococcus_aureus \
  --species escherichia_coli \
  --species streptococcus_pneumoniae
```

On-disk record counts (headers). These are **not** GWAS-significant unitigs:

| File | Records | Origin in this POC |
| --- | ---: | --- |
| `acinetobacter_baumannii_sig_sequences.fasta` | 158 | older pilot dump |
| `campylobacter_jejuni_sig_sequences.fasta` | 59 | older pilot dump |
| `escherichia_coli_sig_sequences.fasta` | 3505 | 3+3 placeholder |
| `klebsiella_pneumoniae_sig_sequences.fasta` | 492 | older pilot dump |
| `neisseria_gonorrhoeae_sig_sequences.fasta` | 192 | older pilot dump |
| `pseudomonas_aeruginosa_sig_sequences.fasta` | 216 | older pilot dump |
| `salmonella_enterica_sig_sequences.fasta` | 195 | older pilot dump |
| `staphylococcus_aureus_sig_sequences.fasta` | 1583 | 3+3 placeholder |
| `streptococcus_pneumoniae_sig_sequences.fasta` | 1569 | 3+3 placeholder |

Placeholder FASTA headers look like `<accession>|<contig_id>`. BLAST then
aligns those contigs back to the same (or other) assemblies, so hit counts
are large and not biologically selected.

---

## 4. BLAST and sequence windows

**Script:** `data_pipeline/pipeline.py`

For each present accession the script builds a BLAST nucleotide DB from the
assembly, aligns the species DBGWAS FASTA, writes a 10-column hit table, and
extracts a `seq_length=1000` window around each hit.

**Inputs**

- `data/camda_amr_2025/processed/assemblies/train/<accession>.fasta`
- `data/camda_amr_2025/processed/balanced_3_3/training_metadata_present.csv`
- `data/dbgwas/p0.05/per_species/<species>_sig_sequences.fasta`

**Outputs**

- `work/blast/dbs/train/<accession>/` (per-isolate BLAST DB)
- `work/blast/hits/train/balanced_3_3_poc/per_species/<accession>_hits.tsv` (53 files)
- `data/camda_amr_2025/processed/datasets/balanced_3_3_poc/sequence_based/per_species/train/full_sequence_dataset.csv`
- `data/camda_amr_2025/processed/datasets/balanced_3_3_poc/sequence_based/per_species/train/<species>/<species>_sequence_dataset.csv`
- `data/camda_amr_2025/processed/balanced_3_3/logs/pipeline.stdout.log`

```bash
pixi run python data_pipeline/pipeline.py \
  --assemblies_dir data/camda_amr_2025/processed/assemblies/train \
  --metadata_path data/camda_amr_2025/processed/balanced_3_3/training_metadata_present.csv \
  --model_type sequence_based \
  --grouping per_species \
  --perspecies_dbgwas_dir data/dbgwas/p0.05/per_species \
  --output_dir data/camda_amr_2025/processed/datasets \
  --run_name balanced_3_3_poc \
  --train
```

`full_sequence_dataset.csv` has **426,394** rows, **53** accessions, both
phenotypes, and no empty `sequence` values. Columns:

`sequence`, `accession`, `query_id`, `hit_count`, `species`, `antibiotic`,
`phenotype`

The full CSV was **not** used for DNABERT training or inference. A 50-row
cap per accession was applied next so a CPU epoch was feasible.

---

## 5. DNABERT

Starting weights are the public Hugging Face substitute
`zhihan1996/DNABERT-2-117M`, cached at `work/models/dnabert2-117m`. That is
not the missing CAMDA `bacteria_model` / DNABERT-MB checkpoint. BPE,
`kmer=-1`, `model_max_length=250`.

One **full** model was trained (not nine per-species models). Downstream RF
tables are still written per species.

### 5.1 Accession-level split and cap

**Script:** `dnabert/finetune/data/prepare_dnabert_splits.py`

**Inputs**

- `data/camda_amr_2025/processed/datasets/balanced_3_3_poc/sequence_based/per_species/train/full_sequence_dataset.csv`

**Outputs** (under `dnabert/finetune/data/balanced_3_3_poc/`)

| File | Contents |
| --- | --- |
| `train.csv` / `train_accessions.txt` | 42 accessions, 2,100 rows |
| `dev.csv` / `dev_accessions.txt` | 11 accessions, 550 rows |
| `test.csv` / `test_accessions.txt` | copy of dev |
| `capped_all.csv` | train+dev concat, 2,650 rows / 53 accessions |

```bash
pixi run python dnabert/finetune/data/prepare_dnabert_splits.py \
  --full_dataset_path data/camda_amr_2025/processed/datasets/balanced_3_3_poc/sequence_based/per_species/train/full_sequence_dataset.csv \
  --out_dir dnabert/finetune/data/balanced_3_3_poc \
  --train_frac 0.80 \
  --drop_dupes \
  --max_per_accession 50
```

**Split bug:** the first attempt stratified on 18 species×class cells. A 20%
holdout is 11 accessions, which is smaller than 18 strata, so
`train_test_split` failed. The script now falls back to stratifying on
binary phenotype only and prints
`stratify=binary_phenotype (species x class holdout is too small)`.

`capped_all.csv` was built after the split as `train.csv` + `dev.csv` so
inference could score every present isolate without the 426k-row CSV.

### 5.2 Species map for inference

**Output:** `work/dnabert_metadata/train_species_to_accession.json`

Built from `training_metadata_present.csv` (`genus_species` → accessions).
`inference.py` needs this map when `--train True`.

### 5.3 Fine-tune (1 epoch, CPU)

**Script:** `dnabert/finetune/train.py`

**Inputs**

- `work/models/dnabert2-117m/`
- `dnabert/finetune/data/balanced_3_3_poc/{train,dev,test}.csv`

**Outputs**

- `dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/checkpoint-132/`
- `dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/best/`
- `dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/results/`

```bash
export WANDB_MODE=offline
export WANDB_NAME=balanced_3_3_dnabert2_117m

pixi run python dnabert/finetune/train.py \
  --model_name_or_path work/models/dnabert2-117m \
  --data_path dnabert/finetune/data/balanced_3_3_poc \
  --kmer -1 \
  --run_name balanced_3_3_dnabert2_117m \
  --model_max_length 250 \
  --per_device_train_batch_size 2 \
  --per_device_eval_batch_size 4 \
  --gradient_accumulation_steps 8 \
  --learning_rate 3e-5 \
  --num_train_epochs 1 \
  --eval_strategy epoch \
  --save_strategy epoch \
  --logging_steps 20 \
  --output_dir dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m \
  --overwrite_output_dir True \
  --log_level info \
  --find_unused_parameters False
```

After weights were written, the DNABERT-2 runtime files were copied into
`best` so inference can load the custom architecture:

```bash
BEST=dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/best
SRC=work/models/dnabert2-117m
cp "$SRC"/bert_layers.py "$SRC"/bert_padding.py \
   "$SRC"/configuration_bert.py "$SRC"/flash_attn_triton.py \
   "$BEST"
```

### 5.4 Fine-tune bugs that had to be fixed

- **`--evaluation_strategy` is ignored** on this transformers build. Use
  `--eval_strategy`.
- **CPU DDP:** `torch.distributed.get_rank()` raised `ValueError`. The trainer
  now treats that as “not distributed.”
- **Unpickleable `TrainingArguments`:** a trailing comma left
  `evaluation_strategy` as a tuple, so `best` could not be loaded after the
  weights were saved. `checkpoint-132` was promoted to `best` and the comma
  was removed.

Do not interpret the 1-epoch loss or F1 as AMR signal.

### 5.5 DNABERT inference → RF tables

**Script:** `dnabert/inference/inference.py`

Scores each capped window, writes `pred_res` / `pred_sus`, then aggregates
by accession and DBGWAS feature into per-species RF CSVs.

**Inputs**

- `dnabert/finetune/data/balanced_3_3_poc/capped_all.csv`
- `dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/best/`
- `data/dbgwas/p0.05/per_species/`
- `work/dnabert_metadata/train_species_to_accession.json`

**Outputs**

- `dnabert/inference/outputs/preds/balanced_3_3_dnabert2_117m/full/train/full_preds.csv`
- `dnabert/inference/outputs/rf_datasets/balanced_3_3_dnabert2_117m/full/train/<species>/<species>_full_rf_dataset.csv`
- `data/camda_amr_2025/processed/balanced_3_3/logs/dnabert_inference.stdout.log`

```bash
pixi run python dnabert/inference/inference.py \
  --output_format random_forest \
  --dataset_dir dnabert/finetune/data/balanced_3_3_poc/capped_all.csv \
  --model_path dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/best \
  --grouping full \
  --run_name balanced_3_3_dnabert2_117m \
  --sig_seqs_dir data/dbgwas/p0.05/per_species \
  --metadata_dir work/dnabert_metadata \
  --base_dir "$PWD" \
  --train True \
  --return_logits sum
```

Always pass `--base_dir`, `--metadata_dir`, and `--sig_seqs_dir`. Script
defaults still point at the old HPC tree.

RF tables (capped inference; both classes; no NaNs):

| Species | Isolates | Feature columns (hit + pred) |
| --- | ---: | ---: |
| `acinetobacter_baumannii` | 6 | 316 |
| `campylobacter_jejuni` | 6 | 118 |
| `escherichia_coli` | 6 | 7010 |
| `klebsiella_pneumoniae` | 6 | 984 |
| `neisseria_gonorrhoeae` | 6 | 384 |
| `pseudomonas_aeruginosa` | 6 | 432 |
| `salmonella_enterica` | 6 | 390 |
| `staphylococcus_aureus` | 6 | 3166 |
| `streptococcus_pneumoniae` | 5 | 3138 |

Scoring the uncapped 426,394-row CSV was left optional and was not done.

---

## 6. Random forest training and inference

This was a plumbing smoke test: serialize nine forests, load them, and write
isolate-level labels. In-sample 53/53 agreement is expected and is **not**
AMR performance (`p >> n`).

### 6.1 Path adapter

`train_rf.py --train_on all` walks `dataset_dir/<fold>/full/train/<species>/`.
The inference tables are flat, so they were symlinked under a fake fold:

**Inputs**

- `dnabert/inference/outputs/rf_datasets/balanced_3_3_dnabert2_117m/full/train/<species>/<species>_full_rf_dataset.csv`

**Outputs**

- `work/rf_poc/folds/fold_0/full/train/<species>/<species>_full_rf_dataset.csv` (symlinks)
- `work/rf_poc/no_flip/<species>.txt` (empty; disables the HPC phenotype-flip default)
- `work/rf_poc/smoke_testing_template.csv` (53-row stand-in for the official CAMDA template)

HPC here means the original CAMDA cluster (SLURM jobs under
`/gpfs/scratch/...`), not a service used in this POC. The terms and why
those paths exist are in [AMR_V2_PIPELINE_GUIDE.md](AMR_V2_PIPELINE_GUIDE.md)
section 2. Defaults that were overridden for the local run:

- `--base_dir` defaults to `/gpfs/scratch/jvaska/CAMDA_AMR/AMR_v2`
- `--flip_phenotype` defaults to a mismatches directory that does not exist here

### 6.2 Train

**Script:** `rf/train_rf.py`

**Inputs**

- `work/rf_poc/folds/`
- `work/rf_poc/no_flip/`

**Outputs**

- `rf/models/rf/per_species/balanced_3_3_poc/<species>_rf_model.joblib` (nine files)

```bash
pixi run python rf/train_rf.py \
  --dataset_dir "$PWD/work/rf_poc/folds" \
  --grouping per_species \
  --dnabert_grouping full \
  --model_name balanced_3_3_poc \
  --train_on all \
  --feature_type both \
  --model_type rf \
  --base_dir "$PWD" \
  --flip_phenotype "$PWD/work/rf_poc/no_flip"
```

`--eval` and `--oof_stack` were not used (no real held-out CAMDA test fold).

**Train bugs**

- Python 3.14 `argparse` rejected a raw `%` in the `--interleave` help
  string. It is now escaped as `%%`.

### 6.3 Execute

**Script:** `rf/infer_rf.py`

**Inputs**

- `rf/models/rf/per_species/balanced_3_3_poc/<species>_rf_model.joblib`
- `dnabert/inference/outputs/rf_datasets/balanced_3_3_dnabert2_117m/full/train/`
- `work/rf_poc/smoke_testing_template.csv`

**Outputs**

- `work/rf_poc/preds/predictions_balanced_3_3_poc.csv` (53 rows: 27 Susceptible, 26 Resistant, 0 NaNs)
- `work/rf_poc/preds/smoke_insample_check.csv` (GT vs pred; in-sample only)

```bash
PYTHONPATH="$PWD/rf" pixi run python rf/infer_rf.py \
  --model_name balanced_3_3_poc \
  --grouping per_species \
  --testing_template "$PWD/work/rf_poc/smoke_testing_template.csv" \
  --models_dir "$PWD/rf/models/rf/per_species/balanced_3_3_poc" \
  --test_dataset_dir "$PWD/dnabert/inference/outputs/rf_datasets/balanced_3_3_dnabert2_117m/full/train" \
  --out_path "$PWD/work/rf_poc/preds/predictions_balanced_3_3_poc.csv" \
  --feature_type both \
  --model_type rf
```

**Infer bug:** `preprocess_df` keeps `ground_truth_phenotype`, and
`model.predict` was called on that frame. The extra column is now dropped
before predict.

---

## Path index

| Step | Main inputs | Main outputs |
| --- | --- | --- |
| Clean metadata | `data/camda_amr_2025/{training,testing}_dataset.csv` | `processed/training_metadata.csv` |
| Select 3+3 | that metadata, existing FASTAs, `bad_accessions.txt` | `processed/balanced_3_3/` |
| Stage | subset metadata + backup pool | `processed/assemblies/train/<acc>.fasta`, `stage_driver.json` |
| Present filter | subset metadata + FASTAs | `training_metadata_present.csv` (53) |
| Strains | present metadata + FASTAs | `processed/dbgwas/balanced_3_3/*.strains.tsv` |
| Feature FASTAs | strains + assemblies | `data/dbgwas/p0.05/per_species/*_sig_sequences.fasta` |
| BLAST / windows | FASTAs + present metadata + feature FASTAs | `work/blast/…/balanced_3_3_poc/`, `datasets/balanced_3_3_poc/…/full_sequence_dataset.csv` |
| DNABERT split | full sequence CSV | `dnabert/finetune/data/balanced_3_3_poc/` |
| DNABERT train | those CSVs + `work/models/dnabert2-117m` | `dnabert/finetune/finetuned_models/balanced_3_3_dnabert2_117m/best` |
| DNABERT infer | `capped_all.csv` + `best` | `dnabert/inference/outputs/{preds,rf_datasets}/balanced_3_3_dnabert2_117m/` |
| RF train | RF CSVs via `work/rf_poc/folds` | `rf/models/rf/per_species/balanced_3_3_poc/*.joblib` |
| RF infer | those models + RF CSVs | `work/rf_poc/preds/predictions_balanced_3_3_poc.csv` |

---

## What this POC does not show

- Real DBGWAS unitigs or published CAMDA feature files
- A complete 3+3 for *S. pneumoniae*
- DNABERT-MB / `bacteria_model` weights
- A GPU multi-epoch fine-tune, or inference on the full 426k windows
- Held-out CAMDA test-set RF evaluation (`--eval` / official testing template)
- Any claim that the isolate-level labels are correct AMR calls

Suggested next work is in the following section.

---

## Next steps for whoever takes this over

This file and [AMR_V2_PIPELINE_GUIDE.md](AMR_V2_PIPELINE_GUIDE.md) are the
starting point. The 3+3 run already proved the plumbing. The current RF and
DNABERT numbers are not a model result, and the full 5,334-accession SRA
queue is not needed to continue from here.

Work should follow this progression. Later stages overwrite earlier artifacts
(especially the shared `data/dbgwas/p0.05/per_species/` FASTAs and anything
BLAST/DNABERT/RF derived from them), so it is worth finishing isolate and
feature decisions before re-running BLAST.

### Worth considering before large jobs

On the machine that ran this POC, WSL filled the virtual disk once
(`ERR124227` left a 0-byte `state.json`). That is not a required first
step on a machine with plenty of free space. If the VHDX is tight, compact
it from **Windows** after exiting WSL:

```bat
wsl --shutdown
Optimize-VHD -Path <path-to-ext4.vhdx> -Mode Full
```

WSL `vmIdleTimeout` / idle-shutdown can also interrupt long SPAdes jobs.
The staging driver is resume-safe (`stage_balanced_subset.py` skips present
FASTAs and known `bad_accession` states), but a disconnect still wastes the
in-flight accession.

If the four SRA resume work dirs and the failed/skip-list work dirs are
still on disk, they are cheaper to keep than to re-prefetch.

### Already in place

These are done and can be reused:

- Clean CAMDA tables in `data/camda_amr_2025/processed/training_metadata.csv`
- 53 staged assemblies in `data/camda_amr_2025/processed/assemblies/train/`
- Present metadata: `data/camda_amr_2025/processed/balanced_3_3/training_metadata_present.csv`
- The train → BLAST → DNABERT → RF path, including the RF fold adapter under
  `work/rf_poc/folds/`
- Pixi `scikit-learn`; the DNABERT CPU/`eval_strategy`/trailing-comma fixes;
  `infer_rf.py` dropping `ground_truth_phenotype`

### 1. Close the pneumococcus hole, then scale the subset

Promote the five failed pneumococcus IDs from “possible bad” into the real
skip list in `data/camda_amr_2025/processed/bad_accessions.txt`. Only
`ERR016650`, `ERR016736`, and `ERR016769` are listed today. Also add
`ERR016777` and `ERR019734` (same `ambiguous_fastq_layout` in
`stage_driver.json`). The driver does **not** append these automatically.

If you want a complete 3+3 before scaling, pick a new *S. pneumoniae* R
isolate from `processed/training_metadata.csv` that is not on the skip list,
stage it with:

```bash
pixi run python data_pipeline/scripts/stage_camda_assemblies.py \
  stage-sra \
  --dataset training \
  --accession <NEW_PNEUMO_R> \
  --threads 4 \
  --memory 16 \
  --resume
```

Then rebuild `training_metadata_present.csv` the same way as in §2.3 (keep
only rows whose `assemblies/train/<accession>.fasta` exists and is nonempty)
and rerun `make_dbgwas_strains.py` on that present table.

After pneumococcus is 3+3, or if 3 S + 2 R is accepted and the work moves
on, a larger class-balanced set can be selected with the same scripts
rather than hand-edited accession lists:

```bash
pixi run python data_pipeline/scripts/select_balanced_subset.py \
  --metadata data/camda_amr_2025/processed/training_metadata.csv \
  --assemblies-dir data/camda_amr_2025/processed/assemblies/train \
  --bad-accession-list data/camda_amr_2025/processed/bad_accessions.txt \
  --output-dir data/camda_amr_2025/processed/balanced_5_5 \
  --quota 5 \
  --backup-quota 5 \
  --seed 42

pixi run python data_pipeline/scripts/stage_balanced_subset.py \
  --subset-dir data/camda_amr_2025/processed/balanced_5_5 \
  --threads 4 \
  --memory 16
```

`5+5` (about 90 isolates, many already staged) is the next useful size.
`10+10` is the step after that. Reuse already-staged FASTAs; the selector
counts them toward the quota. Official CAMDA assembly dumps are still
`TODO` in the README — if those files appear, prefer them over more SRA.

### 2. Replace placeholder DBGWAS FASTAs if DBGWAS can be run

This is the highest-leverage science step if it is available. The nine files
in `data/dbgwas/p0.05/per_species/` are contig dumps, not unitigs. BLAST,
DNABERT, and RF all inherit that.

A prebuilt DBGWAS install (R/`bugwas` plus the Linux tarball noted in
`processed/balanced_3_3/dbgwas_logs/install_attempt.json`) would look like
this, using the **present** strains TSVs:

```bash
DBGWAS \
  -strains data/camda_amr_2025/processed/dbgwas/balanced_3_3/<species>.strains.tsv \
  -output data/camda_amr_2025/work/dbgwas/balanced_3_3/<species> \
  -nb-cores 4

pixi run python data_pipeline/scripts/extract_dbgwas_sig_sequences.py \
  --input data/camda_amr_2025/work/dbgwas/balanced_3_3/<species> \
  --output data/dbgwas/p0.05/per_species/<species>_sig_sequences.fasta \
  --species <species> \
  --p-cutoff 0.05 \
  --top-n 50 \
  --overwrite
```

`fill_missing_dbgwas_fastas.py` can also take `--dbgwas <binary>` with
`--overwrite`. If DBGWAS still fails, the placeholders can stay; that
should be recorded rather than inventing significance cutoffs.

If the subset grows to 5+5, rebuild strains TSVs from the new present
metadata before running DBGWAS. After any FASTA change, BLAST and the
later model steps should be rerun. The six older Sep-17 FASTAs are not
from the 3+3 isolates; they would need replacing too if those species stay
in the new subset.

If the official CAMDA `*_sig_sequences.fasta` files are ever published,
those are preferable to both placeholders and a local DBGWAS run.

### 3. Re-run BLAST and the sequence dataset after features or isolates change

```bash
pixi run python data_pipeline/pipeline.py \
  --assemblies_dir data/camda_amr_2025/processed/assemblies/train \
  --metadata_path data/camda_amr_2025/processed/balanced_<N>_<N>/training_metadata_present.csv \
  --model_type sequence_based \
  --grouping per_species \
  --perspecies_dbgwas_dir data/dbgwas/p0.05/per_species \
  --output_dir data/camda_amr_2025/processed/datasets \
  --run_name balanced_<N>_<N>_poc \
  --train
```

A **new** `--run_name` keeps hits in `work/blast/hits/train/<run_name>/`
and avoids mixing with `balanced_3_3_poc`. BLAST DBs and the sequence CSV
are large; free disk is worth checking first if space is already tight.

### 4. Fine-tune DNABERT as a model, not a smoke test

`work/models/dnabert2-117m` is the public 117M substitute, not DNABERT-MB.
It can stay in use unless `bacteria_model` is recovered.

Suggested next training job, after a new sequence dataset exists:

1. Re-run `dnabert/finetune/data/prepare_dnabert_splits.py` on the new
   `full_sequence_dataset.csv`. Raise `--max_per_accession` above 50 if you
   have a GPU; keep a cap on CPU. If species×class stratification fails,
   the script already falls back to binary phenotype — that is expected on
   small subsets.
2. Rebuild `work/dnabert_metadata/train_species_to_accession.json` from the
   new present metadata (not the old 3+3 map).
3. Fine-tune with `--eval_strategy` (not `--evaluation_strategy`),
   `WANDB_MODE=offline` unless you intend to log, and
   `--model_name_or_path work/models/dnabert2-117m`. Historical CAMDA
   settings in the guide are 24 epochs and a much larger effective batch;
   use a GPU if you have one.
4. Copy `bert_layers.py`, `bert_padding.py`, `configuration_bert.py`, and
   `flash_attn_triton.py` from `work/models/dnabert2-117m` into the new
   `best/` directory.
5. Run `dnabert/inference/inference.py` with `--grouping full`,
   `--base_dir "$PWD"`, `--metadata_dir work/dnabert_metadata`,
   `--sig_seqs_dir data/dbgwas/p0.05/per_species`, and `--train True`.
   Point `--dataset_dir` at the new sequence CSV (or a new capped file).
   The 426,394-row 3+3 CSV was never scored; scoring it is optional and
   only useful if you still care about that exact set.

`--model_path` should be a fine-tuned `best` directory for any result that
will be interpreted; the untuned `work/models/dnabert2-117m` checkpoint is
only a smoke-test base.

### 5. Train RF on a real split, then execute

Keep `--base_dir "$PWD"` and disable `--flip_phenotype` (empty per-species
files, as in `work/rf_poc/no_flip/`). `--train_on all` still expects
`dataset_dir/fold_*/<dnabert_grouping>/train/<species>/`; reuse the symlink
layout in `work/rf_poc/folds/` or add a real `fold_0`/`fold_1`/`fold_2`
split once you have enough isolates.

A useful RF job:

1. Train per-species models on the new RF tables (`--feature_type both`,
   `--model_type rf`, no `--eval` until you have a held-out fold).
2. Hold out isolates by accession, not by sequence row.
3. For official CAMDA scoring, run `rf/infer_rf.py` against the challenge
   testing template (not `work/rf_poc/smoke_testing_template.csv`) and a
   test RF dataset built from testing assemblies. The official template was
   not on disk for the POC.

In-sample accuracy on 5–6 isolates is not a metric and does not belong
in a results table.

### 6. Optional recoveries, if they appear

- Official CAMDA training/testing assemblies named `<accession>.fasta` —
  drop SRA staging for those IDs.
- Official or recovered `data/dbgwas/p0.05/per_species/*_sig_sequences.fasta`.
- The missing DNABERT-MB `bacteria_model` directory. If it arrives, use it
  as `--model_name_or_path` instead of 117M and re-fine-tune.

None of these are required to keep using the tooling.

### A reasonable first week

1. If disk is tight, compact the VHDX and confirm free space; skip this
   if the machine already has room.
2. Add `ERR016777` and `ERR019734` to `bad_accessions.txt`.
3. Decide whether to finish pneumococcus 3+3 or jump to `balanced_5_5`.
4. Decide whether to run real DBGWAS or stay on placeholders, and record
   which in the next write-up.
5. After those isolate and feature choices, rerun BLAST → DNABERT → RF
   with a new `--run_name` / `--model_name` so the 3+3 smoke artifacts
   stay intact for comparison.

