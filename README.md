# SHAPE-Markov-modeling

This repository contains the scripts used to estimate the Markov transition-score matrices used by RNAstructure.

`est_trans_score.py` estimates the Markov transition-score matrices used by `Fold` (`--mc-matrix-a/b`, `--mc-pu-matrix-a/b`) in RNAstructure. SHAPE reactivities are split into low/high states by a threshold, and each transition is scored as the log-odds of its probability in NI-Paired vs. I-Paired nucleotides (`--scheme in`) or Unpaired vs. Paired nucleotides (`--scheme pu`), counted over RNAs with known structures.

## Preparation

```bash
git clone https://github.com/AviranLab/SHAPE-Markov-modeling.git
```

Download the [RNAstructure](https://rna.urmc.rochester.edu/RNAstructureDownload.html) source code and organize it into this folder as `RNAstructure/`.

```
estimate_Markov_Score/
├── data/
│   ├── weeks_set_domains.fa
│   ├── weeks_set_domains.dot
│   └── weeks_set_domains.shape
├── src/
│   └── est_trans_score.py
├── results/
└── RNAstructure/
    └── data_tables/
```

The `weeks_set_domains.*` inputs contain 34 records in the order of the original `rna_list_all.txt`: non-rRNA sequences and separate rRNA domains. They preserve the original per-record FASTA sequences, structures from `dotfile_withoutpk`, and SHAPE numeric precision. The original `weeks_set.*` files are retained for reference; use the domain inputs below to reproduce the legacy whole-dataset scores.

## Quick start

Create and activate a new environment for the steps below:

```bash
conda create -n markov_score python=3.11 numpy pandas -y
conda activate markov_score
```

### 1. Write the score matrices to `results/`

```bash
for scheme in in pu; do
  dir=MarkovScore_whole; [ "$scheme" = pu ] && dir=MarkovScore_PorU_whole
  for model in 1order 1dropend 2order; do
    for t in $(seq 10 59); do
      python src/est_trans_score.py \
        --dot data/weeks_set_domains.dot --fa data/weeks_set_domains.fa --shape data/weeks_set_domains.shape \
        --model $model --scheme $scheme --threshold 0.$t \
        --out results/$dir/score_${model/1dropend/1orddropend}_thd$t.csv
    done
  done
done
```

### 2. Move the results to `RNAstructure/data_tables/`

```bash
mkdir -p RNAstructure/data_tables
mv results/MarkovScore_whole results/MarkovScore_PorU_whole RNAstructure/data_tables/
```
