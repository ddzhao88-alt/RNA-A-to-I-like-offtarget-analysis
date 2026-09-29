#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import argparse
import shutil
import subprocess
import pandas as pd
from common import load_config, resolve


def run(cmd):
    subprocess.run([str(x) for x in cmd], check=True)


def make_exon_bed(gff3, out_bed):
    intervals = defaultdict(list)
    with open(gff3) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            if end > start:
                intervals[fields[0]].append((start, end))

    if not intervals:
        raise RuntimeError(f"No exon features found in {gff3}")

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(out_bed, "w") as out:
        for chrom in sorted(intervals):
            vals = sorted(intervals[chrom])
            s0, e0 = vals[0]
            for s, e in vals[1:]:
                if s <= e0:
                    e0 = max(e0, e)
                else:
                    out.write(f"{chrom}\t{s0}\t{e0}\n")
                    s0, e0 = s, e
            out.write(f"{chrom}\t{s0}\t{e0}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    root = Path.cwd()
    cfg = load_config(args.config)
    fasta = resolve(root, cfg["reference"]["fasta"])
    gff3 = resolve(root, cfg["reference"]["gff3"])
    exon_bed = resolve(root, cfg["reference"]["exon_bed"])
    star_index = resolve(root, cfg["reference"]["star_index"])
    sample_table = resolve(root, cfg["samples"])
    results = resolve(root, cfg["results"])

    for d in [root / "logs", results / "00_check", results / "01_fastp", results / "02_star_bam", results / "03_mpileup_candidates", results / "04_summary", results / "05_figures", star_index]:
        d.mkdir(parents=True, exist_ok=True)

    for tool in ["fastp", "STAR", "samtools", "python"]:
        if not shutil.which(tool):
            raise RuntimeError(f"Missing executable: {tool}")

    for f in [fasta, gff3, sample_table]:
        if not f.exists() or f.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty file: {f}")

    samples = pd.read_csv(sample_table, sep=r"\s+", engine="python")
    required = {"sample", "group", "r1", "r2"}
    if not required.issubset(samples.columns):
        raise RuntimeError("sample_fastq.tsv must contain sample, group, r1 and r2 columns")
    if samples["sample"].duplicated().any():
        raise RuntimeError("Duplicate sample names detected")

    for _, row in samples.iterrows():
        r1 = resolve(root, row["r1"])
        r2 = resolve(root, row["r2"])
        if not r1.exists() or r1.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty R1: {r1}")
        if not r2.exists() or r2.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty R2: {r2}")

    if not Path(str(fasta) + ".fai").exists():
        run(["samtools", "faidx", fasta])

    dict_file = fasta.with_suffix(".dict")
    if not dict_file.exists():
        run(["samtools", "dict", fasta, "-o", dict_file])

    make_exon_bed(gff3, exon_bed)
    print(f"Samples: {len(samples)}")
    print(f"Reference: {fasta}")
    print(f"Exon BED: {exon_bed}")


if __name__ == "__main__":
    main()
