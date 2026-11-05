#!/usr/bin/env python3
"""
Step 4: list finished DMRichR comparisons and write the inputs for the overlap/UpSet script
(05_dmrichr_overlap.R, with helper_overlap.R next to it).

Writes, in --out_dir:
  manifest.tsv  one row per finished comparison: set_id, contrast (= run), sex (= comparison), path
  groups.tsv    which sets are overlapped together (group_id, set_id)
  run_overlap.sh the Rscript command, ready to sbatch or run

In the R script's terms, "outcome" = run (e.g. DMR_AZ_04) and "timepoint" = comparison, so its
built-in groupings mean:
  within_outcome_across_timepoints  -> all comparisons within one run
  within_timepoint_across_outcomes  -> the same comparison across runs (e.g. AZ_04 vs AZ_05)
groups.tsv adds the LACOFD-specific groupings below (per run):
  exposure_by_career    recruit_post_vs_pre + instructor_post_vs_pre
  career_by_timepoint   pre_recruit_vs_instructor + post_recruit_vs_instructor
  baseline_effect       pre_recruit_vs_instructor + pre_recruit_vs_instructor_withbaseline

Usage (from the repository root):
  python3 scripts/dmrichr/04_make_overlap_sets.py --runs DMR_AZ_04,DMR_AZ_05
  # with the regioneR permutation test of locus overlaps (100 permutations, 8 cores):
  python3 scripts/dmrichr/04_make_overlap_sets.py --runs DMR_AZ_04 --perm_test
  bash results/dmrichr/overlap/DMR_AZ_04__DMR_AZ_05/run_overlap.sh      # or sbatch it

A comparison is included only if <comparison>/<dmr_file> exists (finished runs only).
"""

import argparse
import os
import re
import sys
from pathlib import Path


def find_repo_root():
    p = Path(__file__).resolve().parent
    for q in [p, *p.parents]:
        if (q / '.git').exists():
            return q
    return Path(__file__).resolve().parents[2]


REPO_ROOT = find_repo_root()
GROUPS = {
    'exposure_by_career': ['recruit_post_vs_pre', 'instructor_post_vs_pre'],
    'career_by_timepoint': ['pre_recruit_vs_instructor', 'post_recruit_vs_instructor'],
    'baseline_effect': ['pre_recruit_vs_instructor', 'pre_recruit_vs_instructor_withbaseline'],
}


def short(run):
    """DMR_AZ_04 -> AZ04 (keeps set_ids short in plots)."""
    return re.sub(r'^DMR_', '', run).replace('_', '')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', required=True, help='comma list of run IDs, e.g. DMR_AZ_04,DMR_AZ_05')
    ap.add_argument('--run_root', default='results/dmrichr')
    ap.add_argument('--comparisons', default=None, help='comma list to restrict comparisons (default: all finished)')
    ap.add_argument('--dmr_file', default='DMRs/DMRs_annotated.xlsx')
    ap.add_argument('--out_dir', default=None, help='default <run_root>/overlap/<runs joined by __>')
    ap.add_argument('--r_script', default='scripts/dmrichr/05_dmrichr_overlap.R')
    ap.add_argument('--overlap_mode', default='both', help='coord | gene | both')
    ap.add_argument('--perm_test', action='store_true',
                    help='add the regioneR permutation test of locus (coordinate) overlaps')
    ap.add_argument('--perm_n', type=int, default=100, help='permutations per test [default 100]')
    ap.add_argument('--perm_cores', type=int, default=8, help='cores for regioneR [default 8]')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()]
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    root = REPO_ROOT / args.run_root
    out = Path(args.out_dir) if args.out_dir else root / 'overlap' / '__'.join(runs)
    if not out.is_absolute():
        out = REPO_ROOT / out

    rows, missing = [], []
    for run in runs:
        rd = root / run
        if not rd.is_dir():
            sys.exit(f"Run folder not found: {rd}")
        for comp in sorted(p.name for p in rd.iterdir() if p.is_dir()):
            if keep and comp not in keep:
                continue
            if not (rd / comp / 'sample_info.xlsx').exists():
                continue                                    # not a comparison folder
            if (rd / comp / args.dmr_file).exists():
                rows.append((f"{short(run)}_{comp}", run, comp, str(rd / comp)))
            else:
                missing.append(f"{run}/{comp}")
    if not rows:
        sys.exit("No finished comparisons found (no " + args.dmr_file + ")")

    out.mkdir(parents=True, exist_ok=True)
    with open(out / 'manifest.tsv', 'w') as f:
        f.write('set_id\tcontrast\tsex\tpath\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')

    have = {(r[1], r[2]): r[0] for r in rows}
    grp_lines = []
    for run in runs:
        for gid, comps in GROUPS.items():
            sids = [have[(run, c)] for c in comps if (run, c) in have]
            if len(sids) >= 2:
                grp_lines += [f"{short(run)}_{gid}\t{s}" for s in sids]
    with open(out / 'groups.tsv', 'w') as f:
        f.write('group_id\tset_id\n' + ''.join(l + '\n' for l in grp_lines))

    rs = REPO_ROOT / args.r_script
    cmd = (f"Rscript {rs} \\\n  --manifest {out / 'manifest.tsv'} \\\n  --groups_file {out / 'groups.tsv'} \\\n"
           f"  --out_dir {out} \\\n  --overlap_mode {args.overlap_mode} \\\n  --upset_genes 5")
    if args.perm_test:
        if args.overlap_mode == 'gene':
            sys.exit("--perm_test needs --overlap_mode coord or both")
        cmd += f" \\\n  --perm_test --perm_n {args.perm_n} --perm_cores {args.perm_cores}"
    (out / 'run_overlap.sh').write_text(
        "#!/bin/bash\n#SBATCH --job-name=dmr_overlap\n#SBATCH --account=lasallegrp\n#SBATCH --partition=high\n"
        f"#SBATCH --cpus-per-task={args.perm_cores if args.perm_test else 2}\n#SBATCH --mem=32G\n"
        f"#SBATCH --time={'24:00:00' if args.perm_test else '04:00:00'}\n"
        f"#SBATCH --output={out}/overlap_%j.log\n\nset -euo pipefail\n"
        "# activate the same pixi environment as DMRichR (needs GenomicRanges, ggplot2, UpSet packages)\n"
        "# eval \"$(pixi shell-hook --manifest-path <pixi.toml>)\"\n\n" + cmd + "\n")

    print(f"{len(rows)} finished comparisons:")
    for r in rows:
        print(f"  {r[0]}")
    if missing:
        print("Not finished (no DMR file), skipped: " + ', '.join(missing))
    print(f"\nGroups: {len(set(l.split(chr(9))[0] for l in grp_lines))} custom "
          f"(+ all comparisons per run, + each comparison across runs if >1 run)")
    print(f"Wrote {os.path.relpath(out, REPO_ROOT)}/manifest.tsv, groups.tsv, run_overlap.sh\n\n{cmd}")


if __name__ == '__main__':
    main()