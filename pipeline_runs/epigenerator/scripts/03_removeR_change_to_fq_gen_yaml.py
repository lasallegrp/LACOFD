#!/usr/bin/env python3
"""
Step 3 of 3: give FASTQs their final CpG_Me2 names, place them in 01_raw_sequences/,
and write task_samples.yaml.

Run AFTER 01_create_transfers_list_csv.py and 02_data_renaming.py have given every
lane/run file a unique name.

Usage (run from the run directory, e.g. pipeline_runs/epigenerator/EPI_AZ_01; the shared
scripts live one level up). Pass every folder that directly contains FASTQs; the shell glob
*/*/ expands to the run/lane folders (soft-linked folders work too):

  # Dry run (default): prints the plan and a YAML preview, changes nothing
  python3 ../scripts/03_removeR_change_to_fq_gen_yaml.py \\
      01_raw_sequences/2026_LACOFD_WGBS_cellfree_Logan/*/*/ \\
      --target 01_raw_sequences

  # Apply
  python3 ../scripts/03_removeR_change_to_fq_gen_yaml.py \\
      01_raw_sequences/2026_LACOFD_WGBS_cellfree_Logan/*/*/ \\
      --target 01_raw_sequences --apply

Pipeline:
  1) Plan:  {prefix}_R1_001.fastq.gz -> {prefix}_1.fq.gz ; {prefix}_R2_001.fastq.gz -> {prefix}_2.fq.gz
            (files already named {prefix}_1/_2.fastq.gz are accepted too).
            Undetermined_* files are skipped unless --include-undetermined.
  2) Check: abort if two source files map to the same final name (non-unique names).
  3) Place: hard-link (default), symlink, or move each file into --target under its final name.
            Hard links keep the originals in the run folders without using extra disk space;
            source and target must be on the same filesystem.
  4) YAML:  samples = final *.fq.gz names in --target minus _1/_2.fq.gz
            (e.g. LA001_S1_L004_JLLW_Nova1479P_Williams_L2_1.fq.gz -> LA001_S1_L004_JLLW_Nova1479P_Williams_L2),
            written to --yaml-dir (default: parent of --target, i.e. the project directory).

Re-running is safe: files already linked into --target are detected and skipped.
"""

import argparse
import glob
import os
import re
import shutil
import sys
from collections import defaultdict
from typing import Dict, Iterable, Optional, Set

import yaml

PAT_ORIG = re.compile(r'^(?P<prefix>.+)_(?P<read>R[12])_001\.fastq\.gz$')  # Illumina / renamed by step 2
PAT_DONE = re.compile(r'^(?P<prefix>.+)_(?P<read>[12])\.fastq\.gz$')      # already *_1/_2.fastq.gz
PAT_FQ = re.compile(r'^(?P<base>.+)_(?P<read>[12])\.fq\.gz$')             # final form in target


def final_name(fname: str) -> Optional[str]:
    """Return {prefix}_{1|2}.fq.gz for a recognised source filename, else None."""
    m = PAT_ORIG.match(fname)
    if m:
        return f"{m.group('prefix')}_{m.group('read')[1]}.fq.gz"
    m = PAT_DONE.match(fname)
    if m:
        return f"{m.group('prefix')}_{m.group('read')}.fq.gz"
    return None


def reads_by_sample(names: Iterable[str]) -> Dict[str, Set[str]]:
    """Map sample name (filename minus _1/_2.fq.gz) -> read numbers present."""
    reads = defaultdict(set)
    for n in names:
        m = PAT_FQ.match(n)
        if m:
            reads[m.group('base')].add(m.group('read'))
    return reads


def place(src: str, dest: str, mode: str) -> None:
    if mode == "hardlink":
        os.link(src, dest)
    elif mode == "symlink":
        os.symlink(src, dest)  # src is absolute
    else:
        shutil.move(src, dest)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Finalize FASTQ names (_R1/_R2_001.fastq.gz -> _1/_2.fq.gz), hard-link them into "
                    "--target, and write task_samples.yaml (sample names taken from final filenames)."
    )
    ap.add_argument("directories", nargs="+", help="Folders that directly contain the FASTQs (not searched recursively)")
    ap.add_argument("--target", required=True, help="Destination directory, normally <project>/01_raw_sequences")
    ap.add_argument("--apply", action="store_true", help="Actually create links/moves and write YAML (default: dry run)")
    ap.add_argument("--mode", choices=["hardlink", "symlink", "move"], default="hardlink",
                    help="How to place files in --target (default: hardlink)")
    ap.add_argument("--overwrite", action="store_true", help="Replace a different file already at the destination")
    ap.add_argument("--include-undetermined", action="store_true", help="Also place Undetermined_* reads (default: skip)")
    ap.add_argument("--genome", default="hg38", help="Genome written to the YAML (default: hg38)")
    ap.add_argument("--yaml-out", default="task_samples.yaml", help="YAML filename (default: task_samples.yaml)")
    ap.add_argument("--yaml-dir", default=None,
                    help="Where to write the YAML (default: parent of --target, i.e. the project directory)")
    args = ap.parse_args()

    target = os.path.abspath(args.target)
    if not os.path.isdir(target):
        sys.exit(f"❌ Target dir does not exist: {target}")
    yaml_dir = os.path.abspath(args.yaml_dir) if args.yaml_dir else os.path.dirname(target)
    yaml_path = os.path.join(yaml_dir, args.yaml_out)

    # ---- Step 1: plan -------------------------------------------------------------------
    planned = []  # (src_abs, dest_abs)
    unmatched, undetermined = [], []
    for d in args.directories:
        if not os.path.isdir(d):
            print(f"⚠️  Skipping {d} (not a directory)")
            continue
        for src in sorted(glob.glob(os.path.join(d, "*.fastq.gz"))):
            fname = os.path.basename(src)
            if fname.startswith("Undetermined") and not args.include_undetermined:
                undetermined.append(src)
                continue
            new = final_name(fname)
            if new is None:
                unmatched.append(src)
                continue
            planned.append((os.path.abspath(src), os.path.join(target, new)))
    planned = list(dict.fromkeys(planned))  # drop repeats if a folder was passed twice

    if not planned:
        sys.exit("❌ No FASTQs matched. Check that you passed the folders that directly contain *.fastq.gz.")

    # ---- Step 2: uniqueness check -------------------------------------------------------
    by_dest = defaultdict(list)
    for src, dest in planned:
        by_dest[dest].append(src)
    clashes = {dest: srcs for dest, srcs in by_dest.items() if len(srcs) > 1}
    if clashes:
        print("❌ Non-unique filenames: these source files would get the same final name in the target.")
        print("   Run 01_create_transfers_list_csv.py and 02_data_renaming.py first so every lane/run")
        print("   has a unique name. Nothing was changed.")
        for dest, srcs in sorted(clashes.items()):
            print(f"   {os.path.basename(dest)}")
            for s in srcs:
                print(f"      <- {s}")
        sys.exit(1)

    # ---- Step 3: place files ------------------------------------------------------------
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

        verb = f"{'overwrite with ' if replacing else ''}{args.mode}"
        if not args.apply:
            print(f"{pre}would {verb}: {src} -> {dest}")
            n["would"] += 1
            continue
        try:
            if replacing:
                os.remove(dest)
            place(src, dest, args.mode)
            print(f"✅ {verb}: {src} -> {dest}")
            n["done"] += 1
        except OSError as e:
            hint = ""
            if args.mode == "hardlink":
                hint = "\n   (hard links need source and target on the same filesystem; otherwise use --mode symlink)"
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
                  "(files already placed are skipped).")
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
    print(f"Source files matched:         {len(planned)}")
    print(f"Unmatched names (skipped):    {len(unmatched)}")
    print(f"Undetermined (skipped):       {len(undetermined)}")
    if args.apply:
        print(f"Placed:                       {n['done']}")
    else:
        print(f"Would place:                  {n['would']}")
    print(f"Already in place:             {n['already']}")
    print(f"Destination conflicts:        {n['exists']}")
    print(f"Sources missing:              {n['missing']}")
    print(f"Errors:                       {n['errors']}")
    print(f"YAML samples:                 {len(samples)}")
    if unmatched:
        print("\nUnmatched files:")
        for u in unmatched:
            print(f"   - {u}")


if __name__ == "__main__":
    main()