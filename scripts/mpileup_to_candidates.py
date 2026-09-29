#!/usr/bin/env python3
import argparse
import re
import sys

BASES = ("A", "C", "G", "T")


def base_counts(pileup, ref):
    counts = {b: 0 for b in BASES}
    i = 0
    while i < len(pileup):
        c = pileup[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            m = re.match(r"\d+", pileup[i:])
            if m:
                n = int(m.group(0))
                i += len(m.group(0)) + n
            continue
        if c in ".,":
            counts[ref] += 1
            i += 1
            continue
        b = c.upper()
        if b in counts:
            counts[b] += 1
        i += 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_cov", type=int, default=10)
    ap.add_argument("--min_alt", type=int, default=3)
    ap.add_argument("--min_fraction", type=float, default=0.05)
    args = ap.parse_args()

    with open(args.out, "w") as out:
        out.write("sample\tchrom\tpos\tref\tsubstitution\tcoverage\tA\tC\tG\tT\talt_base\talt_reads\tfrequency\n")
        for line in sys.stdin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, ref, _, pileup = fields[:5]
            ref = ref.upper()
            if ref not in {"A", "T"}:
                continue
            counts = base_counts(pileup, ref)
            cov = sum(counts.values())
            if cov < args.min_cov:
                continue
            alt = "G" if ref == "A" else "C"
            substitution = "A>G" if ref == "A" else "T>C"
            alt_reads = counts[alt]
            fraction = alt_reads / cov if cov else 0.0
            if alt_reads < args.min_alt or fraction < args.min_fraction:
                continue
            out.write(
                f"{args.sample}\t{chrom}\t{pos}\t{ref}\t{substitution}\t{cov}\t"
                f"{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\t"
                f"{alt}\t{alt_reads}\t{fraction:.6f}\n"
            )


if __name__ == "__main__":
    main()
