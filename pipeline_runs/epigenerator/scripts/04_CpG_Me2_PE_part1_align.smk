#!/usr/bin/env python3
# =============================================================================================
# CpG_Me2 (paired-end) PART 1: trim + align, one job per lane per sample.
#
# Same trim_galore and bismark commands as 02_CpG_Me2_PE. Differences:
#   * The sample list comes from --configfile (a chunk from 03_split_task_samples_chunks.py, or
#     the full task_samples.yaml), so several instances can run side by side on different
#     SLURM accounts, writing into the same 02_trimmed/ and 04_aligned/.
#   * The account/partition come from --config (defaults: publicgrp/low):
#       trim_account, trim_partition, align_account, align_partition
#   * Resources tuned from 1,903 earlier lane alignments (median 9.9 h, peak RAM 60 GB for 99%,
#     max 79 GB, max 24.6 h): align asks for 70 GB / 36 h first and 100 GB / 72 h on the
#     automatic retry (profile restart-times: 1). Trimming asks for 2 GB (peak 0.34 GB).
#
# Usage (from the run directory, one screen session per instance):
#   snakemake -s ../scripts/04_CpG_Me2_PE_part1_align.smk --profile 00_slurm/ \
#     --configfile chunks/chunk_01.yaml \
#     --config align_account=lasallegrp align_partition=high
#
# Note: no f-strings in this file; Snakemake 7's Snakefile parser mangles them on Python 3.12.
# =============================================================================================

import sys
from snakemake.exceptions import WorkflowError

GB = 1024  # MB per GB

for key in ("samples", "genome"):
    if not config.get(key):
        raise WorkflowError("Config key '" + key + "' is missing. Run with --configfile "
                            "chunks/chunk_XX.yaml (or task_samples.yaml).")

TRIM_ACCOUNT = str(config.get("trim_account", "publicgrp"))
TRIM_PARTITION = str(config.get("trim_partition", "low"))
ALIGN_ACCOUNT = str(config.get("align_account", "publicgrp"))
ALIGN_PARTITION = str(config.get("align_partition", "low"))

print("CpG_Me2 Part 1: {} lanes | trim -> {}/{} | align -> {}/{}".format(
      len(config["samples"]), TRIM_ACCOUNT, TRIM_PARTITION, ALIGN_ACCOUNT, ALIGN_PARTITION), file=sys.stderr)


rule all:
    input:
        expand("02_trimmed/{sample}_1.fq.gz_trimming_report.txt", sample=config["samples"]),
        expand("02_trimmed/{sample}_1_val_1_fastqc.html", sample=config["samples"]),
        expand("02_trimmed/{sample}_1_val_1_fastqc.zip", sample=config["samples"]),
        expand("02_trimmed/{sample}_1_val_1.fq.gz", sample=config["samples"]),
        expand("02_trimmed/{sample}_2.fq.gz_trimming_report.txt", sample=config["samples"]),
        expand("02_trimmed/{sample}_2_val_2_fastqc.html", sample=config["samples"]),
        expand("02_trimmed/{sample}_2_val_2_fastqc.zip", sample=config["samples"]),
        expand("02_trimmed/{sample}_2_val_2.fq.gz", sample=config["samples"]),
        expand("04_aligned/{sample}_1_val_1_bismark_bt2_pe.bam", sample=config["samples"]),
        expand("04_aligned/{sample}_1_val_1_bismark_bt2_PE_report.txt", sample=config["samples"]),


rule trim:
    message: "Trimming samples"
    input:
        r1 = "01_raw_sequences/{sample}_1.fq.gz",
        r2 = "01_raw_sequences/{sample}_2.fq.gz",
    output:
        out1 = "02_trimmed/{sample}_1.fq.gz_trimming_report.txt",
        out2 = "02_trimmed/{sample}_1_val_1_fastqc.html",
        out3 = "02_trimmed/{sample}_1_val_1_fastqc.zip",
        out4 = "02_trimmed/{sample}_1_val_1.fq.gz",
        out5 = "02_trimmed/{sample}_2.fq.gz_trimming_report.txt",
        out6 = "02_trimmed/{sample}_2_val_2_fastqc.html",
        out7 = "02_trimmed/{sample}_2_val_2_fastqc.zip",
        out8 = "02_trimmed/{sample}_2_val_2.fq.gz"
    log: "00_std_err_logs/02_trimmed_{sample}.log"
    benchmark: "00_time_logs/02_trimmed_{sample}.txt"
    resources:
        mem_mb = 2 * GB,          # peak in earlier runs: 0.34 GB (was 15 GB)
        time = 60 * 24 * 2,       # 2 days (unchanged)
        account = TRIM_ACCOUNT,
        partition = TRIM_PARTITION,
    shell: "trim_galore --paired --cores {threads} --2colour 20 --fastqc --clip_r1 10 --clip_r2 20 --three_prime_clip_r1 10 --three_prime_clip_r2 10 --output_dir 02_trimmed/ {input.r1} {input.r2} 2> {log}"


rule align:
    message: "Aligning bisulfite reads"
    input:
        r1 = "02_trimmed/{sample}_1_val_1.fq.gz",
        r2 = "02_trimmed/{sample}_2_val_2.fq.gz"
    output:
        out1 = "04_aligned/{sample}_1_val_1_bismark_bt2_pe.bam",
        out2 = "04_aligned/{sample}_1_val_1_bismark_bt2_PE_report.txt"
    log: "00_std_err_logs/04_aligned_{sample}.txt"
    benchmark: "00_time_logs/04_aligned_{sample}.txt"
    threads: 6
    resources:
        # attempt 1: 70 GB / 36 h ; automatic retry: 100 GB / 72 h
        mem_mb = lambda wildcards, attempt: (70 if attempt == 1 else 100) * GB,
        time = lambda wildcards, attempt: (36 if attempt == 1 else 72) * 60,
        account = ALIGN_ACCOUNT,
        partition = ALIGN_PARTITION,
    shell: "bismark -n 1 --genome 01_genomes/{config[genome]}/ --multicore {threads} --dovetail --score_min L,0,-0.2 --gzip --unmapped -o 04_aligned/ -1 {input.r1} -2 {input.r2} 2> {log}"
