#!/usr/bin/env python3
"""
Step 2 of 3: rename raw FASTQs IN PLACE using the CSV from 01_create_transfers_list_csv.py
(Original_File_Name -> New_File_Name, in the folder given by File_Path).

NOTE: if 01_raw_sequences/{batch} is a soft link to the lab data folder, the files are renamed
in the data folder itself (the link only points there).

Usage (run from the run directory, e.g. pipeline_runs/epigenerator/EPI_AZ_01):
  # Preview
  python3 ../scripts/02_data_renaming.py --file 2026_LACOFD_WGBS_cellfree_Logan_transfers.csv --dry-run
  # Apply
  python3 ../scripts/02_data_renaming.py --file 2026_LACOFD_WGBS_cellfree_Logan_transfers.csv

Optional: --overwrite replaces an existing destination file (uses os.replace).
Keep the CSV: it is the record of original -> new names.
"""

import argparse
import pandas as pd
import os
from collections import Counter

parser = argparse.ArgumentParser(
    description='Rename files using original and new file name columns in CSV'
)
parser.add_argument('--file', required=True, type=str, metavar='<str>', help='Path to CSV')
parser.add_argument('--dry-run', action='store_true', help='Preview actions without renaming')
parser.add_argument('--overwrite', action='store_true', help='Overwrite if destination exists (uses os.replace)')
arg = parser.parse_args()

df = pd.read_csv(arg.file)

required_cols = {'File_Path', 'Original_File_Name', 'New_File_Name'}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"CSV is missing required columns: {', '.join(sorted(missing))}")

ops = []  # (old_path, new_path)
skipped_same = 0

for _, row in df.iterrows():
    file_path = row['File_Path']
    orig = row['Original_File_Name']
    new = row['New_File_Name']

    # Directory where the file currently lives
    directory_path = os.path.dirname(file_path)

    old_full = os.path.join(directory_path, orig)
    new_full = os.path.join(directory_path, new)

    # Normalize
    old_full = os.path.abspath(old_full)
    new_full = os.path.abspath(new_full)

    if old_full == new_full or orig == new:
        skipped_same += 1
        continue

    ops.append((old_full, new_full))

# Check for duplicate targets in CSV (likely a mistake)
targets = [n for _, n in ops]
dupe_targets = [t for t, c in Counter(targets).items() if c > 1]
if dupe_targets:
    print("⚠️  Duplicate target filenames found in CSV (multiple rows map to the same destination). These will be skipped:")
    for t in dupe_targets:
        print(f"   - {t}")
    # Filter out all operations with duplicate targets
    ops = [(o, n) for (o, n) in ops if n not in dupe_targets]

renamed = 0
missing_src = 0
existing_dst = 0
errors = 0

for old_full, new_full in ops:
    if not os.path.exists(old_full):
        print(f"❌ Source not found: {old_full}")
        missing_src += 1
        continue

    dst_exists = os.path.exists(new_full)

    if arg.dry_run:
        action = "would overwrite" if (dst_exists and arg.overwrite) else ("would skip (exists)" if dst_exists else "would rename")
        print(f"[DRY-RUN] {action}: {old_full} -> {new_full}")
        continue

    try:
        if dst_exists and not arg.overwrite:
            print(f"⚠️  Destination exists, skipping (use --overwrite to replace): {new_full}")
            existing_dst += 1
            continue

        # Perform rename
        if arg.overwrite:
            os.replace(old_full, new_full)  # atomic on same filesystem
        else:
            os.rename(old_full, new_full)

        print(f"✅ Renamed: {old_full} -> {new_full}")
        renamed += 1

    except PermissionError as e:
        print(f"❌ Permission error for {old_full} -> {new_full}: {e}")
        errors += 1
    except OSError as e:
        # Covers cross-device moves or other OS-level issues
        print(f"❌ OS error for {old_full} -> {new_full}: {e}")
        errors += 1

print("\n===== Summary =====")
print(f"Total rows in CSV:      {len(df)}")
print(f"Planned renames:        {len(ops)}")
print(f"Skipped (same names):   {skipped_same}")
if dupe_targets:
    print(f"Skipped (dupe targets): {len(dupe_targets)} unique targets")
print(f"Sources missing:        {missing_src}")
print(f"Skipped (dst exists):   {existing_dst} {'(use --overwrite to replace)' if not arg.overwrite else ''}")
print(f"Errors:                 {errors}")
print(f"{'Would rename (dry-run):' if arg.dry_run else 'Renamed:'} {renamed}")