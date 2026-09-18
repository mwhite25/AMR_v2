# AMR Prediction Pipeline using DNABERT and Random Forest

## Introduction
- Find accepted ISMB/ECCB [extended abstract here](https://github.com/jaaxk/CAMDA_AMR/blob/main/AMR_extended_abstract.pdf)
- This pipeline was developed for the CAMDA AMR Prediction Challenge 2025 with 9 pathogenic species each treated on one of 4 antibiotics (GEN, ERY, CAZ, TET)
- 9 pathogenic species: `Neisseria gonorrhoeae, Staphylococcus aureus, Streptococcus pneumoniae, Salmonella enterica, Klebsiella pneumoniae, Escherichia coli, Pseudomonas aeruginosa, Acinetobacter baumannii, Campylobacter jejuni`
- The train and test data can be found on the [official 2025 CAMDA challenges website](https://bipress.boku.ac.at/camda2025/the-camda-contest-challenges/#amr)
- The predictions are binary (resistant or susceptible)
- This work was presented orally at [ISMB/ECCB 2025 in the CAMDA COSI session](https://www.iscb.org/ismbeccb2025/programme-agenda/scientific-programme/CAMDA)
- The repository contains the code to finetune a pretrained DNABERT model and run inference using DNABERT and/or Random Forest/XGBoost
- **This work is currently under patent consideration as of September 2025**

## Usage

For a beginner-oriented walkthrough of the inputs, commands, intermediate
files, outputs, and troubleshooting checks, see
[AMR_V2_PIPELINE_GUIDE.md](AMR_V2_PIPELINE_GUIDE.md).

The expected input to the pipeline is a directory containing the CAMDA accessions QCed and trimmed with fastp, and assembled using SPAdes with the --isolate option.

Links to model and FASTA downloads are currently unavailable and will be added after publication.

To run with DNABERT and RF:

1. TODO: Make Conda env
2. TODO: (train set only) Run DBGWAS to get significant sequences: `data_pipeline/run_dbgwas.sh`
    - Or [download significant species here](TODO) and place in `data_pipeline/data/dbgwas/per_species`
3. Generate initial dataset: `data_pipeline/run_pipeline.sh` 
    - Fill out arguments including assembly path
4. (train set only) Finetune DNABERT using full dataset generated from step 3: `dnabert/finetune/run_finetune.sh`
    - Class balancing and train/dev split at `dnabert/finetune/data/dnabert_preprocessing.py`
    - Requires pretrained model path. This project used a microbiome-pretrained DNABERT model (DNABERT-MB) available [here](TODO).
    - Or [download finetuned model here](TODO) and place in `dnabert/finetune/finetuned_models/{MODEL_NAME}`
5. Generate RF dataset by getting DNABERT preds and hit counts: `dnabert/inference/run_inference.sh`
6. TODO: (train set only) Train RF/XGBoost model: `rf/train_rf.sh`
7. TODO: (test set only) Run inference on test set and fill out testing template: `rf/run_inference.sh`
8. Evaluate on [CAMDA official website](https://bipress.boku.ac.at/camda2025/competitions/camda-2025/?submissions)

## Options

- These scripts have the option to run with different groupings (per-species, per-antibiotic, and full)
    - Specify these with the --grouping option in `run_pipeline.sh`, `run_inference.sh`, and TODO: RF scripts, and finetune multiple DNABERT-MB models accordingly

- These scripts can also perform 'abalation studies' by removing DNABERT predictions and using only RF/XGBoost (hit counts)
    - Use --model_type matrix_based in `run_inference.sh` and pass directly to TODO: RF scripts

- We can also perform 'abalation studies' by removing species/antibiotic features
    - This can only be done using a full model, so use --grouping full in `run_pipeline.sh` and `run_inference.sh`
    - Also, change LEAKAGE=False in `dnabert/inference/run_inference.sh` to prevent mapping to species-specific DBGWAS features only (maps to features from all species combined)
