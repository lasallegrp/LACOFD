#!/usr/bin/env python3
"""
Step 1 of 2: build a CSV that maps every raw Illumina FASTQ to a UNIQUE final name.
Read-only: nothing in the raw data folder is renamed, moved, or modified.

Illumina filenames ({sample}_S{n}_L{lane}_R{1|2}_001.fastq.gz) only encode sample, sample
index and lane, so the same sample sequenced on another run/flow cell can produce an identical
filename. The final name inserts the name of the folder each file sits in (e.g. the run/lane
folder JLLW_Nova1479P_Williams_L2) and uses the _1/_2.fq.gz ending CpG_Me2 expects:

    LA001_S1_L004_R1_001.fastq.gz   (in .../JLLW_Nova1479P_Williams_L2/)
 -> LA001_S1_L004_JLLW_Nova1479P_Williams_L2_1.fq.gz

If lane folder names repeat across runs, --tag-depth 2 adds the run folder as well:
 -> LA001_S1_L004_260819_DTSA1302_1303_1304_1305_NovaX25B_JLLW_Nova1479P_Williams_L2_1.fq.gz

02_link_fastqs_gen_yaml.py then hard-links each original file into 01_raw_sequences/ under
that final name.

Usage (run from the run directory, e.g. pipeline_runs/epigenerator/EPI_AZ_01; the shared
scripts live one level up in pipeline_runs/epigenerator/scripts/):
python3 ../scripts/01_create_transfers_list_csv.py \
  --indir 01_raw_sequences/2026_LACOFD_WGBS_cellfree_Logan \
  --output 2026_LACOFD_WGBS_cellfree_Logan_transfers.csv
  
CSV columns:
  Sample_ID, Sample_Num, Lane, Read_Dir   parsed from the Illumina name
  Parent_Dir                              folder name(s) used to make the final name unique
  Original_File_Name                      delivered name (never changed)
  File_Path                               absolute path as found (may go through a soft link)
  Real_Path                               absolute path with soft links resolved (true data location)
  New_File_Name                           final name in 01_raw_sequences/ ({sample}_1.fq.gz / _2.fq.gz)
  YAML_Sample                             entry that goes into task_samples.yaml
"""

import argparse
import os
import re
import sys

import pandas as pd

parser = argparse.ArgumentParser(
    description='Generate a CSV mapping raw FASTQs to unique final names (read-only; changes no files)'
)
parser.add_argument('--indir', required=True, type=str, metavar='<str>',
                    help='Input directory (searched recursively, soft-linked folders followed)')
parser.add_argument('--output', required=True, type=str, metavar='<str>', help='Name of output CSV')
parser.add_argument('--tag-depth', type=int, default=1, metavar='<int>',
                    help='How many enclosing folder names to add to each file name (default: 1 = the '
                         'folder the file is in). Use 2 to add run + lane folder if lane folder names repeat.')
arg = parser.parse_args()
if arg.tag_depth < 1:
    sys.exit("--tag-depth must be >= 1")

columns = [
    'Sample_ID', 'Sample_Num', 'Lane', 'Read_Dir', 'Parent_Dir',
    'Original_File_Name', 'File_Path', 'Real_Path', 'New_File_Name', 'YAML_Sample'
]

# Example matched: LA001_S1_L004_R1_001.fastq.gz
pat_regular = re.compile(
    r'^(?P<Sample_ID>[^_]+)_'       # sample ID without underscores (e.g., LA001, 016B)
    r'(?P<Sample_Num>S\d+)_'        # S number (e.g., S1, S20)
    r'(?P<Lane>L\d+)_'              # L with digits (e.g., L004)
    r'(?P<Read_Dir>R[12])_'         # R1 or R2
    r'\d+\.fastq\.gz$'              # trailing chunk number like 001
)

# Undetermined reads, if present. Example: Undetermined_S0_L004_R1_001.fastq.gz
pat_undetermined = re.compile(
    r'^(?P<Sample_ID>Undetermined)_'
    r'(?P<Sample_Num>S\d+)_'
    r'(?P<Lane>L\d+)_'
    r'(?P<Read_Dir>R[12])_'
    r'\d+\.fastq\.gz$'
)

rows = []
unmatched = []

# followlinks=True: also descend into soft-linked run/lane folders
for path, subdirs, files in os.walk(arg.indir, followlinks=True):
    subdirs.sort()
    # Tag = last <tag-depth> folder names between --indir and the file, joined with "_"
    rel_parts = [p for p in os.path.relpath(path, arg.indir).split(os.sep) if p != "."]
    if not rel_parts:  # files directly inside --indir
        rel_parts = [os.path.basename(os.path.normpath(os.path.abspath(arg.indir)))]
    parent_dir = "_".join(rel_parts[-arg.tag_depth:])
    for name in sorted(files):
        if not name.endswith(".fastq.gz"):
            continue

        filepath = os.path.abspath(os.path.join(path, name))

        m = None
        if name.startswith("Undetermined"):
            m = pat_undetermined.match(name)
        if m is None:
            m = pat_regular.match(name)
        if m is None:
            unmatched.append(filepath)
            continue

        sample_id, sample_num, lane, read_dir = (
            m.group('Sample_ID'), m.group('Sample_Num'), m.group('Lane'), m.group('Read_Dir'))

        yaml_sample = f'{sample_id}_{sample_num}_{lane}_{parent_dir}'
        new_name = f'{yaml_sample}_{read_dir[1]}.fq.gz'   # R1 -> _1.fq.gz, R2 -> _2.fq.gz

        rows.append([sample_id, sample_num, lane, read_dir, parent_dir,
                     name, filepath, os.path.realpath(filepath), new_name, yaml_sample])

master_df = pd.DataFrame(rows, columns=columns)
master_df.to_csv(arg.output, index=False)
print(f"Wrote {len(master_df)} rows to {arg.output}")

# Final names must be unique across ALL folders, because step 2 puts every file into one
# directory (01_raw_sequences/). Duplicates mean two folders share the same name.
dupes = master_df.loc[master_df['New_File_Name'].duplicated(keep=False)]
if not dupes.empty:
    sys.stderr.write(
        "⚠️ New_File_Name is not unique: these files would collide in 01_raw_sequences/.\n"
        "   Step 2 will refuse to run. Re-run this script with a larger --tag-depth (e.g. 2 adds the\n"
        "   run folder too). Do not rename folders in the raw data folder:\n"
    )
    for fp in dupes.sort_values('New_File_Name')['File_Path']:
        sys.stderr.write(f"  - {fp}\n")

# Every sample/lane should have both reads (Undetermined is skipped by step 2 by default)
pairs = master_df.loc[master_df['Sample_ID'] != 'Undetermined'].groupby('YAML_Sample')['Read_Dir'].apply(set)
unpaired = pairs[pairs.apply(lambda r: r != {'R1', 'R2'})]
if not unpaired.empty:
    sys.stderr.write("⚠️ Samples missing R1 or R2 (CpG_Me2_PE needs both):\n")
    for s, r in unpaired.items():
        sys.stderr.write(f"  - {s} (has {','.join(sorted(r))})\n")

if unmatched:
    sys.stderr.write("⚠️ Skipped files that didn't match the expected Illumina pattern "
                     "(non-standard or previously renamed names):\n")
    for u in unmatched:
        sys.stderr.write(f"  - {u}\n")