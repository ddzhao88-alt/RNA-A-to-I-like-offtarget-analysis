#!/bin/bash
#SBATCH --job-name=abe_mpileup
#SBATCH --partition=qcpu_23a
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=30G
#SBATCH --output=./logs/%x.%A_%a.out
#SBATCH --error=./logs/%x.%A_%a.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate abe_rna_editing
CONFIG=${CONFIG:-config.yaml}

SAMPLE=$(awk -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' sample_fastq.tsv)
FASTA=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['reference']['fasta']))
PY
)
EXON_BED=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['reference']['exon_bed']))
PY
)
RESULTS=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['results']))
PY
)
MAPQ=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['filters']['mapping_quality'])
PY
)
BASEQ=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['filters']['base_quality'])
PY
)
MINCOV=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['filters']['min_coverage'])
PY
)
MINALT=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['filters']['min_alt_reads'])
PY
)
MINFRAC=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['filters']['min_alt_fraction'])
PY
)

BAM="$RESULTS/02_star_bam/${SAMPLE}.star.sorted.bam"
OUT="$RESULTS/03_mpileup_candidates/${SAMPLE}.ati_candidates.tsv"
mkdir -p "$RESULTS/03_mpileup_candidates"

samtools mpileup -f "$FASTA" -l "$EXON_BED" -q "$MAPQ" -Q "$BASEQ" -B -d 1000000 "$BAM" | python scripts/mpileup_to_candidates.py --sample "$SAMPLE" --out "$OUT" --min_cov "$MINCOV" --min_alt "$MINALT" --min_fraction "$MINFRAC"
