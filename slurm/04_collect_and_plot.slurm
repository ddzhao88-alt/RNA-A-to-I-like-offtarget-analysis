#!/bin/bash
#SBATCH --job-name=abe_collect
#SBATCH --partition=qcpu_23a
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=20G
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate abe_rna_editing
python scripts/collect_and_plot.py --config "${CONFIG:-config.yaml}"
