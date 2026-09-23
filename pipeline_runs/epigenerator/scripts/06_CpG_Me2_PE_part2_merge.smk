#!/usr/bin/env python3
# =============================================================================================
# CpG_Me2 (paired-end) PART 2: merge lanes -> name sort -> deduplicate -> extract methylation
# -> cytosine reports -> MultiQC. Same commands as the Part 2 rules of 02_CpG_Me2_PE.
#
# Differences from 02_CpG_Me2_PE:
#   * Lanes are merged per sample ID (text before the first "_", as before), but the lanes are
#     taken from task_samples.yaml instead of whatever BAMs happen to be in 04_aligned/.
#   * Completeness check: a lane counts as finished only when its BAM, its Bismark report AND
#     its Snakemake benchmark file (written only after a successful job) all exist. By default
#     the workflow STOPS if any lane of any sample is unfinished, so no sample is ever merged
#     from fewer lanes than it has.
#   * --config partial=True processes only samples whose lanes are ALL finished and skips the
#     rest (e.g. to get results for finished samples while other chunks are still aligning).
#     Re-run later without it to process the remaining samples.
#   * Account/partition: --config part2_account=... part2_partition=... (default publicgrp/low).
#
# Usage (from the run directory, after Part 1):
#   snakemake -s ../scripts/06_CpG_Me2_PE_part2_merge.smk --profile 00_slurm/
#
# Note: no f-strings in this file; Snakemake 7's Snakefile parser mangles them on Python 3.12.
# =============================================================================================

import glob
import os
import sys
from collections import defaultdict
from snakemake.common import Mode
from snakemake.exceptions import WorkflowError

configfile: "task_samples.yaml"

PARTIAL = str(config.get("partial", False)).lower() in ("true", "1", "yes")
P2_ACCOUNT = str(config.get("part2_account", "publicgrp"))
P2_PARTITION = str(config.get("part2_partition", "low"))

wildcard_constraints:
    sample_id = "[^_/]+"


def lane_files(lane):
    return ("04_aligned/" + lane + "_1_val_1_bismark_bt2_pe.bam",
            "04_aligned/" + lane + "_1_val_1_bismark_bt2_PE_report.txt",
            "00_time_logs/04_aligned_" + lane + ".txt")


LANES_BY_ID = defaultdict(list)
for _lane in config["samples"]:
    LANES_BY_ID[str(_lane).split("_")[0]].append(str(_lane))

# The sample set is decided ONLY by the main Snakemake process, which builds the DAG. Cluster jobs
# re-read this file later, when more lanes may have finished; they must not re-decide the set.
IN_JOB = workflow.mode != Mode.default

_done, _incomplete = [], {}
for _sid in ([] if IN_JOB else sorted(LANES_BY_ID)):
    _missing = [l for l in LANES_BY_ID[_sid] if not all(os.path.exists(f) for f in lane_files(l))]
    if _missing:
        _incomplete[_sid] = _missing
    else:
        _done.append(_sid)

if IN_JOB:
    pass  # validation and sample selection already done by the main process
elif _incomplete and not PARTIAL:
    _n_lanes = sum(len(v) for v in _incomplete.values())
    _lines = "\n".join("  {}: {} of {} lanes unfinished, e.g. {}".format(
                           sid, len(v), len(LANES_BY_ID[sid]), v[0])
                       for sid, v in list(_incomplete.items())[:20])
    raise WorkflowError(
        "{} sample(s) have unfinished lanes ({} lanes). Part 2 would merge them from too few "
        "lanes, so it will not start.\n{}\nFinish Part 1 for these lanes, or re-run with "
        "--config partial=True to process only the complete samples now.".format(
            len(_incomplete), _n_lanes, _lines))
elif not _done:
    raise WorkflowError("No sample has all of its lanes finished yet; nothing to do.")
elif _incomplete:
    print("CpG_Me2 Part 2 (partial): processing {} complete sample(s), skipping {} with "
          "unfinished lanes: {}{}".format(len(_done), len(_incomplete), ", ".join(list(_incomplete)[:10]),
                                         " ..." if len(_incomplete) > 10 else ""), file=sys.stderr)
else:
    print("CpG_Me2 Part 2: all {} samples complete ({} lanes).".format(
          len(_done), len(config["samples"])), file=sys.stderr)

SAMPLE_IDS = _done


# Collect BAM files for each sample (all of its lanes, from task_samples.yaml)
def get_bam_files(wildcards):
    return [lane_files(l)[0] for l in sorted(LANES_BY_ID[wildcards.sample_id])]


# Final outputs for MultiQC. In the main process: the reports of the selected samples. Inside a
# cluster job: the reports already on disk (a superset of what the main process required).
def get_final_outputs(wildcards=None):
    if IN_JOB:
        return sorted(glob.glob("08_cytosine_reports/*_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz"))
    return expand("08_cytosine_reports/{sample_id}_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
                  sample_id=SAMPLE_IDS)


rule all:
    input:
        get_final_outputs(),
        "09_multiqc/multiqc_report.html"


rule merge_lanes:
    message: "Merging aligned data"
    input: get_bam_files
    output: temp("04_merged_aligned/{sample_id}_merged.bam")
    log: "00_std_err_logs/04_merged_aligned_{sample_id}.txt"
    benchmark: "00_time_logs/04_merged_aligned_{sample_id}.txt"
    threads: 4
    resources:
        mem_mb = 1024 * 25, # Last number is memory in GB
        time = 60 * 24 * 3, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "samtools merge -n -@ {threads} {output} {input} 2> {log}"

rule name_sort:
    message: "Sorting BAM files by read name"
    input: "04_merged_aligned/{sample_id}_merged.bam"
    output: temp("04_merged_aligned/{sample_id}_merged_name_sorted.bam")
    log: "00_std_err_logs/04_name_sorted_{sample_id}.txt"
    benchmark: "00_time_logs/04_name_sorted_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 25, # Last number is memory in GB
        time = 60 * 24 * 2, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "samtools sort -n -@ {threads} -o {output} {input} 2> {log}"

rule deduplicate:
    message: "Removing PCR duplicates"
    input: "04_merged_aligned/{sample_id}_merged_name_sorted.bam"
    output:
        out_bam = "05_deduplicated/{sample_id}_merged_name_sorted.deduplicated.bam",
        out_report = "05_deduplicated/{sample_id}_merged_name_sorted.deduplication_report.txt"
    log: "00_std_err_logs/05_deduplicated_{sample_id}.txt"
    benchmark: "00_time_logs/05_deduplicated_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 50, # Last number is memory in GB
        time = 60 * 24 * 2, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell:
        """
        deduplicate_bismark --bam --paired {input} 2> {log}
        mv {wildcards.sample_id}_merged_name_sorted.deduplicated.bam {output.out_bam}
        mv {wildcards.sample_id}_merged_name_sorted.deduplication_report.txt {output.out_report}
        """

rule sort:
    message: "Sorting BAM files"
    input: "05_deduplicated/{sample_id}_merged_name_sorted.deduplicated.bam"
    output: temp("06_sorted/{sample_id}_merged_name_sorted.deduplicated.sorted.bam")
    log: "00_std_err_logs/06_sorted_{sample_id}.txt"
    benchmark: "00_time_logs/06_sorted_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 25, # Last number is memory in GB
        time = 60 * 24 * 2, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "picard SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate 2> {log}"

rule size_metrics:
    message: "Collecting insert size metrics"
    input: "06_sorted/{sample_id}_merged_name_sorted.deduplicated.sorted.bam"
    output:
        ins = "07_size_metrics/{sample_id}_merged_name_sorted.deduplicated.sorted.bam.insert.txt",
        hist = "07_size_metrics/{sample_id}_merged_name_sorted.deduplicated.sorted.bam.histogram.pdf"
    log: "00_std_err_logs/07_size_metrics_{sample_id}.txt"
    benchmark: "00_time_logs/07_size_metrics_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 25, # Last number is memory in GB
        time = 60 * 24 * 1, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "picard CollectInsertSizeMetrics INPUT={input} OUTPUT={output.ins} HISTOGRAM_FILE={output.hist} ASSUME_SORTED=TRUE 2> {log}"

rule nucleotide_coverage:
    message: "Assessing nucleotide coverage"
    input: "06_sorted/{sample_id}_merged_name_sorted.deduplicated.sorted.bam"
    output: "06_nt_coverage/{sample_id}_merged_name_sorted.deduplicated.sorted.nucleotide_stats.txt"
    log: "00_std_err_logs/06_nt_coverage_{sample_id}.txt"
    benchmark: "00_time_logs/06_nt_coverage_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 25, # Last number is memory in GB
        time = 60 * 24 * 2, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "bam2nuc --genome_folder 01_genomes/{config[genome]}/ {input} --dir 06_nt_coverage 2> {log}"

rule extract_methylation:
    message: "Extracting methylation"
    input: "05_deduplicated/{sample_id}_merged_name_sorted.deduplicated.bam"
    output:
        "06_methylation/CpG_context_{sample_id}_merged_name_sorted.deduplicated.txt.gz",
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated.bedGraph.gz",
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated.bismark.cov.gz",
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated.M-bias_R1.png",
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated.M-bias_R2.png",
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated.M-bias.txt",
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated_splitting_report.txt",
        "06_methylation/Non_CpG_context_{sample_id}_merged_name_sorted.deduplicated.txt.gz"
    log: "00_std_err_logs/06_methylation_{sample_id}.txt"
    benchmark: "00_time_logs/06_methylation_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 50, # Last number is memory in GB
        time = 60 * 24 * 3, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "bismark_methylation_extractor --paired-end --gzip --comprehensive --merge_non_CpG --bedGraph --multicore 6 --buffer_size 34G {input} --output_dir 06_methylation/ 2> {log}"

rule cytosine_reports:
    message: "Generating genome-wide cytosine reports"
    input:
        "06_methylation/{sample_id}_merged_name_sorted.deduplicated.bismark.cov.gz"
    output:
        "08_cytosine_reports/{sample_id}_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz"
    log: "00_std_err_logs/08_cytosine_reports_{sample_id}.txt"
    benchmark: "00_time_logs/08_cytosine_reports_{sample_id}.txt"
    resources:
        mem_mb = 1024 * 30, # Last number is memory in GB
        time = 60 * 24 * 3, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "coverage2cytosine --dir 08_cytosine_reports --gzip --merge_CpG {input} --genome_folder 01_genomes/{config[genome]}/ --output {wildcards.sample_id}_merged_name_sorted.deduplicated.bismark.cov.gz 2> {log}"

rule qc_report:
    message: "Generating QC reports"
    input:
        get_final_outputs()
    output:
        "09_multiqc/multiqc_report.html"
    log: "00_std_err_logs/09_multiqc.txt"
    benchmark: "00_time_logs/09_multiqc.txt"
    resources:
        mem_mb = 1024 * 50, # Last number is memory in GB
        time = 60 * 24 * 2, # Last number is days
        account = P2_ACCOUNT,
        partition = P2_PARTITION,
    shell: "multiqc . --ignore raw_sequences --ignore 01_raw_sequences --outdir 09_multiqc/ 2> {log}"
