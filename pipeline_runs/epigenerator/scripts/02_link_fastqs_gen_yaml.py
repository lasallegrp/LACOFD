#!/usr/bin/env python3
"""
Step 2 of 2: hard-link each raw FASTQ into 01_raw_sequences/ under its unique final name
(from the CSV made by 01_create_transfers_list_csv.py), then write task_samples.yaml.

The raw data is never renamed, moved, or modified. A hard link is a second name for the same
file: it uses no extra disk space, and the original keeps its delivered name.

Usage (run from the run directory, e.g. pipeline_runs/epigenerator/EPI_AZ_01):

  # Dry run (default): prints the plan and a YAML preview, changes nothing
  python3 ../scripts/02_link_fastqs_gen_yaml.py \\
      --csv 2026_LACOFD_WGBS_cellfree_Logan_transfers.csv \\
      --target 01_raw_sequences

  # Apply
  python3 ../scripts/02_link_fastqs_gen_yaml.py \\
      --csv 2026_LACOFD_WGBS_cellfree_Logan_transfers.csv \\
      --target 01_raw_sequences --apply

Pipeline:
  1) Read File_Path -> New_File_Name from the CSV
     (Undetermined_* rows skipped unless --include-undetermined).
  2) Check: abort if two rows share a final name (non-unique names across lanes/runs).
  3) Link:  hard-link (default) or soft-link each original file to <target>/<New_File_Name>.
            Hard links need the data and target on the same filesystem.
  4) YAML:  samples = final *.fq.gz names in --target minus _1/_2.fq.gz, written to --yaml-dir
            (default: parent of --target, i.e. the run directory).

Re-running is safe: links that already point at the right file are detected and skipped.
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, Set

import pandas as pd
import yaml

PAT_FQ = re.compile(r'^(?P<base>.+)_(?P<read>[12])\.fq\.gz$')  # final form in target


def reads_by_sample(names: Iterable[str]) -> Dict[str, Set[str]]:
    """Map sample name (filename minus _1/_2.fq.gz) -> read numbers present."""
    reads = defaultdict(set)
    for n in names:
        m = PAT_FQ.match(n)
        if m:
            reads[m.group('base')].add(m.group('read'))
    return reads


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Hard-link raw FASTQs into --target under unique final names from the CSV "
                    "(raw data untouched) and write task_samples.yaml."
    )
    ap.add_argument("--csv", required=True, help="CSV from 01_create_transfers_list_csv.py")
    ap.add_argument("--target", required=True, help="Destination directory, normally <run>/01_raw_sequences")
    ap.add_argument("--apply", action="store_true", help="Actually create links and write YAML (default: dry run)")
    ap.add_argument("--mode", choices=["hardlink", "symlink"], default="hardlink",
                    help="Link type (default: hardlink; use symlink if hard links fail across filesystems)")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace a DIFFERENT file already at the destination (only the link is replaced)")
    ap.add_argument("--include-undetermined", action="store_true", help="Also link Undetermined_* reads (default: skip)")
    ap.add_argument("--genome", default="hg38", help="Genome written to the YAML (default: hg38)")
    ap.add_argument("--yaml-out", default="task_samples.yaml", help="YAML filename (default: task_samples.yaml)")
    ap.add_argument("--yaml-dir", default=None,
                    help="Where to write the YAML (default: parent of --target, i.e. the run directory)")
    args = ap.parse_args()

    target = os.path.abspath(args.target)
    if not os.path.isdir(target):
        sys.exit(f"❌ Target dir does not exist: {target}")
    yaml_dir = os.path.abspath(args.yaml_dir) if args.yaml_dir else os.path.dirname(target)
    yaml_path = os.path.join(yaml_dir, args.yaml_out)

    # ---- Step 1: read the plan from the CSV ---------------------------------------------
    df = pd.read_csv(args.csv)
    missing_cols = {"File_Path", "New_File_Name", "Sample_ID"} - set(df.columns)
    if missing_cols:
        sys.exit(f"❌ CSV is missing columns {sorted(missing_cols)}. "
                 "Re-create it with the current 01_create_transfers_list_csv.py.")
    bad = df.loc[~df["New_File_Name"].astype(str).str.match(PAT_FQ)]
    if not bad.empty:
        sys.exit(f"❌ {len(bad)} New_File_Name value(s) don't end in _1.fq.gz/_2.fq.gz, "
                 f"e.g. {bad['New_File_Name'].iloc[0]}. Re-create the CSV with the current step 1 script.")

    undetermined = df["Sample_ID"].astype(str).eq("Undetermined")
    n_undet = int(undetermined.sum())
    if not args.include_undetermined:
        df = df.loc[~undetermined]

    planned = [(os.path.abspath(src), os.path.join(target, new))
               for src, new in zip(df["File_Path"], df["New_File_Name"])]
    if not planned:
        sys.exit("❌ No rows to link.")

    # ---- Step 2: uniqueness check -------------------------------------------------------
    by_dest = defaultdict(list)
    for src, dest in planned:
        by_dest[dest].append(src)
    clashes = {dest: srcs for dest, srcs in by_dest.items() if len(srcs) > 1}
    if clashes:
        print("❌ Non-unique final names: these source files would get the same name in the target.")
        print("   Re-run step 1 with a larger --tag-depth (e.g. --tag-depth 2). Nothing was changed.")
        for dest, srcs in sorted(clashes.items()):
            print(f"   {os.path.basename(dest)}")
            for s in srcs:
                print(f"      <- {s}")
        sys.exit(1)

    # ---- Step 3: link -------------------------------------------------------------------
    n = defaultdict(int)
    pre = "" if args.apply else "[DRY-RUN] "
    for src, dest in planned:
        if not os.path.exists(src):
            print(f"❌ Source missing: {src}")
            n["missing"] += 1
            continue
        if os.path.lexists(dest):
            if os.path.exists(dest) and os.path.samefile(src, dest):
                n["already"] += 1  # already linked by an earlier run
                continue
            if not args.overwrite:
                print(f"⚠️  {pre}Different file already at destination, skipping (use --overwrite): {dest}")
                n["exists"] += 1
                continue
            replacing = True
        else:
            replacing = False

        verb = f"{'replace with ' if replacing else ''}{args.mode}"
        if not args.apply:
            print(f"{pre}would {verb}: {src} -> {dest}")
            n["would"] += 1
            continue
        try:
            if replacing:
                os.remove(dest)  # removes only the name in --target, never the raw data
            if args.mode == "hardlink":
                os.link(src, dest)
            else:
                os.symlink(os.path.realpath(src), dest)
            print(f"✅ {verb}: {src} -> {dest}")
            n["done"] += 1
        except OSError as e:
            hint = ""
            if args.mode == "hardlink":
                hint = "\n   (hard links need the data and target on the same filesystem; otherwise use --mode symlink)"
            print(f"❌ {args.mode} failed: {src} -> {dest}: {e}{hint}")
            n["errors"] += 1

    # ---- Step 4: YAML from final files in target ----------------------------------------
    in_target = {os.path.basename(p) for p in glob.glob(os.path.join(target, "*.fq.gz"))}
    names = in_target if args.apply else in_target | {os.path.basename(d) for _, d in planned}
    reads = reads_by_sample(names)
    samples = sorted(reads)
    unpaired = sorted(s for s, r in reads.items() if r != {"1", "2"})
    this_batch = {PAT_FQ.match(os.path.basename(d)).group("base") for _, d in planned}
    other = sorted(set(samples) - this_batch)

    yaml_text = yaml.dump(
        {"sequence_data": target, "genome": args.genome, "samples": samples},
        default_flow_style=False, sort_keys=False, explicit_start=True,
    )

    if other:
        print(f"\nℹ️  {len(other)} sample(s) already in {target} from other batches will also be in the YAML, e.g.:")
        for s in other[:5]:
            print(f"   - {s}")
    if unpaired:
        print(f"\n⚠️  {len(unpaired)} sample(s) missing _1 or _2 (CpG_Me2_PE needs both):")
        for s in unpaired:
            print(f"   - {s} (has read {','.join(sorted(reads[s]))})")

    problems = n["missing"] + n["errors"]
    if args.apply:
        if problems:
            print(f"\n❌ YAML NOT written: {problems} file(s) failed. Fix the errors above and re-run "
                  "(files already linked are skipped).")
        else:
            if os.path.exists(yaml_path):
                print(f"\nℹ️  Overwriting existing {yaml_path}")
            with open(yaml_path, "w") as yh:
                yh.write(yaml_text)
            print(f"\n✅ YAML written: {yaml_path}")
    else:
        print(f"\nYAML preview (would be written to {yaml_path}):\n")
        print(yaml_text)

    # ---- Summary ------------------------------------------------------------------------
    print("===== Summary =====")
    print(f"Mode:                         {args.mode}{'' if args.apply else ' (dry run)'}")
    print(f"CSV rows used:                {len(planned)}")
    print(f"Undetermined rows:            {n_undet} {'(included)' if args.include_undetermined else '(skipped)'}")
    if args.apply:
        print(f"Linked:                       {n['done']}")
    else:
        print(f"Would link:                   {n['would']}")
    print(f"Already linked:               {n['already']}")
    print(f"Destination conflicts:        {n['exists']}")
    print(f"Sources missing:              {n['missing']}")
    print(f"Errors:                       {n['errors']}")
    print(f"YAML samples:                 {len(samples)}")


if __name__ == "__main__":
    main()