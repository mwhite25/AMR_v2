#!/bin/bash
mkdir -p data/camda_amr_2025/processed/pilot_logs

for species in \
  neisseria_gonorrhoeae \
  staphylococcus_aureus \
  streptococcus_pneumoniae \
  salmonella_enterica \
  klebsiella_pneumoniae \
  escherichia_coli \
  pseudomonas_aeruginosa \
  acinetobacter_baumannii \
  campylobacter_jejuni
do
  log="data/camda_amr_2025/processed/pilot_logs/${species}.log"
  manifest="data/camda_amr_2025/processed/pilot_logs/${species}.csv"

  pixi run python data_pipeline/scripts/stage_camda_assemblies.py stage-sra \
    --dataset training \
    --species "$species" \
    --limit 1 \
    --threads 4 \
    --memory 16 \
    --resume \
    --work-dir data/camda_amr_2025/work/assemblies \
    --manifest "$manifest" \
    --bad-accession-list data/camda_amr_2025/processed/bad_accessions.txt \
    --skip-bad-accessions \
    2>&1 | tee "$log"
done