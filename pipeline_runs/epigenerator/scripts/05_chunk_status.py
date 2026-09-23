#!/usr/bin/env python3
r"""
Show Part 1 progress per chunk, so you know when an account is free for the next chunk.

A lane counts as trimmed/aligned when Snakemake's benchmark file for that step exists
(written only after the job succeeded) and, for alignment, the BAM exists.

Usage (run from the run directory):
python3 ../scripts/05_chunk_status.py
python3 ../scripts/05_chunk_status.py --chunks-dir chunks

Tip: name each screen session after its chunk and account (e.g. screen -S c01_lasallegrp);
`screen -ls` then shows which chunk runs where, and `squeue -u $USER -o "%.10a %.9P %.30j %.8T"`
shows the jobs per account.
"""

import argparse
import glob
import os
import statistics
import sys


def read_samples(path):
    """Read the samples list from a chunk YAML (format written by script 03)."""
    samples, in_samples = [], False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("samples:"):
                in_samples = True
            elif in_samples and s.startswith("- "):
                samples.append(s[2:].strip().strip("'\""))
            elif s and not s.startswith("#") and ":" in s and not s.startswith("- "):
                in_samples = False
    return samples


def bench_hours(path):
    try:
        with open(path) as fh:
            fh.readline()
            return float(fh.readline().split("\t")[0]) / 3600
    except (OSError, ValueError, IndexError):
        return None


def main():
    ap = argparse.ArgumentParser(description="Per-chunk progress of CpG_Me2 Part 1 (trim + align).")
    ap.add_argument("--chunks-dir", default="chunks", help="Directory with chunk_*.yaml (default: chunks/)")
    args = ap.parse_args()

    chunk_files = sorted(glob.glob(os.path.join(args.chunks_dir, "chunk_*.yaml")))
    if not chunk_files:
        sys.exit(f"❌ No chunk_*.yaml in {args.chunks_dir}/. Run 03_split_task_samples_chunks.py first.")

    print(f"{'chunk':<10}{'lanes':>6}{'trimmed':>9}{'aligned':>9}   status")
    tot_l = tot_t = tot_a = 0
    hours = []
    for cf in chunk_files:
        lanes = read_samples(cf)
        t = sum(os.path.exists(f"00_time_logs/02_trimmed_{l}.txt") for l in lanes)
        a = 0
        for l in lanes:
            b = f"00_time_logs/04_aligned_{l}.txt"
            if os.path.exists(b) and os.path.exists(f"04_aligned/{l}_1_val_1_bismark_bt2_pe.bam"):
                a += 1
                h = bench_hours(b)
                if h is not None:
                    hours.append(h)
        if a == len(lanes):
            status = "DONE - account free for next chunk"
        elif t == 0 and a == 0:
            status = "not started"
        else:
            status = f"{100 * a / len(lanes):.0f}% aligned"
        name = os.path.splitext(os.path.basename(cf))[0]
        print(f"{name:<10}{len(lanes):>6}{t:>9}{a:>9}   {status}")
        tot_l, tot_t, tot_a = tot_l + len(lanes), tot_t + t, tot_a + a

    print(f"{'TOTAL':<10}{tot_l:>6}{tot_t:>9}{tot_a:>9}   {100 * tot_a / tot_l:.0f}% aligned")
    if hours:
        print(f"\nAlignment time so far: median {statistics.median(hours):.1f} h, "
              f"max {max(hours):.1f} h over {len(hours)} lanes")


if __name__ == "__main__":
    main()