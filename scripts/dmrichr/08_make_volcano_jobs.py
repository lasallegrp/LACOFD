#!/usr/bin/env python3
"""
Step 8: volcano plots (08_volcano_plot.R) for every finished comparison in one or more runs.

Lists comparisons that have DMRs/background_annotated.xlsx and writes, in --out_dir:
  volcano_manifest.tsv   one row per comparison: task, run, comparison, background, outdir, prefix
  run_volcano.sbatch     SLURM array job, one task per comparison (reads the manifest)

Each plot goes to <run>/<comparison>/DMRs/volcano/<comparison>_volcano_full.pdf/.png.
With --collect, the job also copies every PNG into <out_dir>/png/ so you can flip through them.

Usage (from the repository root):
  python3 scripts/dmrichr/08_make_volcano_jobs.py --runs DMR_AZ_04
  python3 scripts/dmrichr/08_make_volcano_jobs.py --runs DMR_AZ_04,DMR_AZ_05 --collect
  python3 scripts/dmrichr/08_make_volcano_jobs.py --runs DMR_AZ_04 --sig_by auto --n_labels 30
  sbatch results/dmrichr/volcano/DMR_AZ_04/run_volcano.sbatch

--sig_by pvalue (default) colours every region with p < 0.05. If DMRichR called this
comparison's DMRs by q (>= 100 regions with q < 0.05), the coloured count will be larger
than the DMRs in DMRs.bed; --sig_by auto colours exactly DMRichR's DMRs instead.
"""

import argparse
import os
import shlex
import sys
from pathlib import Path


def find_repo_root():
    p = Path(__file__).resolve().parent
    for q in [p, *p.parents]:
        if (q / '.git').exists():
            return q
    return Path(__file__).resolve().parents[2]


REPO_ROOT = find_repo_root()
DEFAULT_PIXI = 'reproducibility/dmrichr/pixi/pixi.toml'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', required=True, help='comma list of run IDs, e.g. DMR_AZ_04,DMR_AZ_05')
    ap.add_argument('--run_root', default='results/dmrichr')
    ap.add_argument('--comparisons', default=None, help='comma list to restrict comparisons (default: all finished)')
    ap.add_argument('--out_dir', default=None, help='default <run_root>/volcano/<runs joined by __>')
    ap.add_argument('--r_script', default='scripts/dmrichr/08_volcano_plot.R')
    ap.add_argument('--pixi_manifest', default=DEFAULT_PIXI)
    ap.add_argument('--sig_by', default='pvalue', choices=['pvalue', 'qvalue', 'auto'],
                    help='pvalue (default) = colour regions with p < --cutoff; qvalue = q < --cutoff; '
                         'auto = DMRichR rule (q if >= 100 regions pass q, else p)')
    ap.add_argument('--cutoff', type=float, default=0.05)
    ap.add_argument('--n_labels', type=int, default=20)
    ap.add_argument('--collect', action='store_true', help='also copy all PNGs into <out_dir>/png/')
    ap.add_argument('--force', action='store_true', help='redo comparisons that already have a volcano plot')
    ap.add_argument('--mem', default='16G')
    ap.add_argument('--time', default='02:00:00')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()]
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    root = REPO_ROOT / args.run_root
    out = Path(args.out_dir) if args.out_dir else root / 'volcano' / '__'.join(runs)
    if not out.is_absolute():
        out = REPO_ROOT / out
    rs = REPO_ROOT / args.r_script
    pixi = Path(args.pixi_manifest)
    if not pixi.is_absolute():
        pixi = REPO_ROOT / pixi
    for f, label in [(rs, 'R script'), (pixi, 'pixi manifest')]:
        if not f.exists():
            print(f"WARNING: {label} not found here: {f} (fine if you run this off Hive)")

    rows, missing, done = [], [], []
    for run in runs:
        rd = root / run
        if not rd.is_dir():
            sys.exit(f"Run folder not found: {rd}")
        for comp in sorted(p.name for p in rd.iterdir() if p.is_dir()):
            if keep and comp not in keep:
                continue
            cd = rd / comp
            if not (cd / 'sample_info.xlsx').exists():
                continue                                        # not a comparison folder
            bg = cd / 'DMRs' / 'background_annotated.xlsx'
            if not bg.exists():
                missing.append(f"{run}/{comp}")
                continue
            od = cd / 'DMRs' / 'volcano'
            if (od / f"{comp}_volcano_full.pdf").exists() and not args.force:
                done.append(f"{run}/{comp}")
                continue
            rows.append((run, comp, str(bg), str(od), comp))
    if not rows:
        msg = "Nothing to do."
        if done:
            msg += f" Already plotted (use --force to redo): {', '.join(done)}"
        if missing:
            msg += f" Not finished: {', '.join(missing)}"
        sys.exit(msg)

    out.mkdir(parents=True, exist_ok=True)
    (out / 'logs').mkdir(exist_ok=True)
    man = out / 'volcano_manifest.tsv'
    with open(man, 'w') as f:
        f.write('task\trun\tcomparison\tbackground\toutdir\tprefix\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(r) + '\n')

    collect = (f'\nmkdir -p {shlex.quote(str(out / "png"))}\n'
               f'cp "$OUTDIR/${{PREFIX}}_volcano_full.png" '
               f'{shlex.quote(str(out / "png"))}/"${{RUN}}__${{PREFIX}}_volcano_full.png"\n') if args.collect else ''
    sb = f"""#!/bin/bash
#SBATCH --job-name=dmr_volcano
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={out}/logs/volcano_%A_%a.log
# Written by 08_make_volcano_jobs.py: one array task per comparison in volcano_manifest.tsv
set -euo pipefail

MANIFEST={shlex.quote(str(man))}
LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' "$MANIFEST")
[ -n "$LINE" ] || {{ echo "No manifest row for task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK RUN COMP BACKGROUND OUTDIR PREFIX <<< "$LINE"
echo "Task $TASK: $RUN / $COMP  ($(date))"

eval "$(pixi shell-hook --manifest-path {shlex.quote(str(pixi))})"
for p in optparse readxl ggplot2 ggrepel dplyr; do
  Rscript -e "q(status = !requireNamespace('$p', quietly = TRUE))" >/dev/null 2>&1 \\
    || {{ echo "R package $p is missing from the pixi env"; exit 1; }}
done

mkdir -p "$OUTDIR"
Rscript {shlex.quote(str(rs))} \\
  --background "$BACKGROUND" \\
  --outdir "$OUTDIR" \\
  --prefix "$PREFIX" \\
  --sig_by {args.sig_by} \\
  --cutoff {args.cutoff} \\
  --n_labels {args.n_labels}
{collect}
echo "Done: $OUTDIR/${{PREFIX}}_volcano_full.pdf"
"""
    (out / 'run_volcano.sbatch').write_text(sb)

    print(f"{len(rows)} comparisons to plot:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}/{r[1]}")
    if done:
        print("Already plotted, skipped (--force to redo): " + ', '.join(done))
    if missing:
        print("Not finished (no background_annotated.xlsx), skipped: " + ', '.join(missing))
    rel = os.path.relpath(out, REPO_ROOT)
    print(f"\nWrote {rel}/volcano_manifest.tsv and run_volcano.sbatch\nSubmit with:\n  sbatch {rel}/run_volcano.sbatch")


if __name__ == '__main__':
    main()