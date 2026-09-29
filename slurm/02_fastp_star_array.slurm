#!/bin/bash
#SBATCH --job-name=abe_align
#SBATCH --partition=qcpu_23a
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=50G
#SBATCH --output=./logs/%x.%A_%a.out
#SBATCH --error=./logs/%x.%A_%a.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate abe_rna_editing
CONFIG=${CONFIG:-config.yaml}

read SAMPLE GROUP R1 R2 <<< "$(awk -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1,$2,$3,$4}' sample_fastq.tsv)"
THREADS=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['alignment']['threads'])
PY
)
FASTP_THREADS=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['alignment']['fastp_threads'])
PY
)
MINLEN=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(cfg['alignment']['min_read_length'])
PY
)
INDEX=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['reference']['star_index']))
PY
)
TWOPASS=$(python - "$CONFIG" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(str(cfg['alignment'].get('star_twopass', False)).lower())
PY
)
RESULTS=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1])); print(Path(cfg['results']))
PY
)

FASTP_DIR="$RESULTS/01_fastp"
BAM_DIR="$RESULTS/02_star_bam"
OUTDIR="$BAM_DIR/$SAMPLE"
mkdir -p "$FASTP_DIR" "$OUTDIR"

fastp -i "$R1" -I "$R2" -o "$FASTP_DIR/${SAMPLE}.R1.clean.fq.gz" -O "$FASTP_DIR/${SAMPLE}.R2.clean.fq.gz" --detect_adapter_for_pe --qualified_quality_phred 20 --length_required "$MINLEN" --thread "$FASTP_THREADS" --html "$FASTP_DIR/${SAMPLE}.fastp.html" --json "$FASTP_DIR/${SAMPLE}.fastp.json"

STAR_ARGS=(--runThreadN "$THREADS" --genomeDir "$INDEX" --readFilesIn "$FASTP_DIR/${SAMPLE}.R1.clean.fq.gz" "$FASTP_DIR/${SAMPLE}.R2.clean.fq.gz" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$OUTDIR/")
if [[ "$TWOPASS" == "true" ]]; then
  STAR_ARGS+=(--twopassMode Basic)
fi
STAR "${STAR_ARGS[@]}"

samtools index -@ 4 "$OUTDIR/Aligned.sortedByCoord.out.bam"
ln -sfn "$PWD/$OUTDIR/Aligned.sortedByCoord.out.bam" "$BAM_DIR/${SAMPLE}.star.sorted.bam"
ln -sfn "$PWD/$OUTDIR/Aligned.sortedByCoord.out.bam.bai" "$BAM_DIR/${SAMPLE}.star.sorted.bam.bai"
