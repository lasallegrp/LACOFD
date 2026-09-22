#!/usr/bin/env python3
"""
Step 1 of 3: build a CSV that maps every raw Illumina FASTQ to a UNIQUE new name.

Illumina filenames ({sample}_S{n}_L{lane}_R{1|2}_001.fastq.gz) only encode sample, sample
index and lane, so the same sample sequenced on another run/flow cell can produce an identical
filename. This script inserts the name of the folder each file sits in (e.g. the run/lane
folder JLLW_Nova1479P_Williams_L2) between the lane and the read direction:

    LA001_S1_L004_R1_001.fastq.gz
 -> LA001_S1_L004_JLLW_Nova1479P_Williams_L2_R1_001.fastq.gz

Nothing is renamed here; 02_data_renaming.py applies the CSV.

Usage (run from the run directory, e.g. pipeline_runs/epigenerator/EPI_AZ_01; the shared
scripts live one level up in pipeline_runs/epigenerator/scripts/):
python3 ../scripts/01_create_transfers_list_csv.py \
  --indir 01_raw_sequences/2026_LACOFD_WGBS_cellfree_Logan \
  --output 2026_LACOFD_WGBS_cellfree_Logan_transfers.csv
"""

import argparse
import os
import re
import sys

import pandas as pd

parser = argparse.ArgumentParser(
    description='Generate CSV mapping raw FASTQs to unique new names (tracks file processing progress)'
)
parser.add_argument('--indir', required=True, type=str, metavar='<str>', help='Input directory (searched recursively)')
parser.add_argument('--output', required=True, type=str, metavar='<str>', help='Name of output CSV')
arg = parser.parse_args()

columns = [
    'Sample_ID', 'Sample_Num', 'Lane', 'Read_Dir',
    'File_Path', 'Parent_Dir', 'Original_File_Name', 'New_File_Name',
    'Transferred', 'Processed'
]

# Patterns
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

# followlinks=True: also descend into soft-linked run/lane folders (raw data linked from the lab data folder)
for path, subdirs, files in os.walk(arg.indir, followlinks=True):
    subdirs.sort()
    parent_dir = os.path.basename(os.path.normpath(path))
    for name in sorted(files):
        if not name.endswith(".fastq.gz"):
            continue

        # Absolute path so 02_data_renaming.py works from any working directory
        filepath = os.path.abspath(os.path.join(path, name))

        m = None
        if name.startswith("Undetermined"):
            m = pat_undetermined.match(name)
        if m is None:
            m = pat_regular.match(name)

        if m is None:
            unmatched.append(filepath)
            continue

        sample_id = m.group('Sample_ID')
        sample_num = m.group('Sample_Num')
        lane = m.group('Lane')
        read_dir = m.group('Read_Dir')

        # New name: insert parent_dir between lane and read_dir
        newname = f'{sample_id}_{sample_num}_{lane}_{parent_dir}_{read_dir}_001.fastq.gz'

        rows.append([
            sample_id, sample_num, lane, read_dir,
            filepath, parent_dir, name, newname, "NA", "NA"
        ])

master_df = pd.DataFrame(rows, columns=columns)
master_df.to_csv(arg.output, index=False)
print(f"Wrote {len(master_df)} rows to {arg.output}")

# New names must be unique across ALL folders, because step 3 puts every file into one
# directory (01_raw_sequences/). Duplicates mean two folders share the same name.
dupes = master_df.loc[master_df['New_File_Name'].duplicated(keep=False)]
if not dupes.empty:
    sys.stderr.write(
        "⚠️ New_File_Name is not unique: these files would collide in 01_raw_sequences/.\n"
        "   Give their parent folders distinct names, then re-run this script:\n"
    )
    for fp in dupes.sort_values('New_File_Name')['File_Path']:
        sys.stderr.write(f"  - {fp}\n")

# Report any files we skipped (already-renamed files end up here too)
if unmatched:
    sys.stderr.write("⚠️ Skipped files that didn't match the expected Illumina pattern "
                     "(already renamed, or non-standard names):\n")
    for u in unmatched:
        sys.stderr.write(f"  - {u}\n")