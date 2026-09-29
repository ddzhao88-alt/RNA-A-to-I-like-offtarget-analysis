#!/bin/bash
#SBATCH --job-name=abe_star_index
#SBATCH --partition=qcpu_23a
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=100G
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate abe_rna_editing

CONFIG=${CONFIG:-config.yaml}
FASTA=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['reference']['fasta']))
PY
)
INDEX=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['reference']['star_index']))
PY
)
THREADS=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['alignment']['threads'])
PY
)

mkdir -p "$INDEX"
STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$INDEX" --genomeFastaFiles "$FASTA" --genomeSAindexNbases 14
