# RNA A-to-I-like off-target analysis

This repository contains the RNA-seq workflow used to identify A-to-I-like editing candidates in mouse, human and cynomolgus monkey samples.

## Workflow

1. Adapter and quality trimming with fastp.
2. Alignment to the species-specific reference genome with STAR.
3. BAM indexing with SAMtools.
4. Restriction to merged annotated exonic regions.
5. Candidate calling from `samtools mpileup` output.
6. Retention of A>G and reverse-complement T>C substitutions with coverage >=10, alternative reads >=3 and alternative-base fraction >=0.05.
7. Per-sample candidate burden summarization and figure generation.

## Reference genomes

- Mouse: Mus musculus GRCm39, Ensembl release 107.
- Human: Homo sapiens GRCh38, Ensembl release 107.
- Cynomolgus monkey: Macaca fascicularis T2T-MFA8v1.1, GCF_037993035.2.

Reference genomes and annotations are not included in this repository.

## Environment

```bash
conda env create -f environment.yml
conda activate abe_rna_editing
```

## Project setup

Copy one species configuration to the project root:

```bash
cp configs/human.yaml config.yaml
```

Create `sample_fastq.tsv` and `group_colors.tsv` using the files in `examples/` as templates.

Expected sample table format:

```text
sample  group  r1  r2
```

FASTQ files, reference files and analysis results should remain outside version control.

## Run

```bash
mkdir -p logs results
sbatch slurm/00_check_project.slurm
sbatch slurm/01_build_star_index.slurm
sbatch --array=1-N slurm/02_fastp_star_array.slurm
sbatch --array=1-N slurm/03_mpileup_candidates_array.slurm
sbatch slurm/04_collect_and_plot.slurm
```

Replace `N` with the number of samples in `sample_fastq.tsv`.

To use a non-default configuration:

```bash
CONFIG=configs/mouse.yaml sbatch slurm/00_check_project.slurm
```

The same `CONFIG` value should be used for all steps.

## Main outputs

```text
results/03_mpileup_candidates/*.ati_candidates.tsv
results/04_summary/all_A_to_I_like_candidates.tsv
results/04_summary/per_sample_candidate_burden.tsv
results/04_summary/per_sample_candidate_burden_GraphPad_wide.tsv
results/05_figures/per_sample_A_to_I_like_candidates.pdf
results/05_figures/per_sample_A_to_I_like_candidates.svg
results/05_figures/per_sample_A_to_I_like_candidates.png
```

PDF and SVG outputs retain editable text. PNG output is saved at 600 dpi with a transparent background.
