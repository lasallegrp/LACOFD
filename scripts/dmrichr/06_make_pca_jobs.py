#!/usr/bin/env python3
"""
Step 6 job builder: PCA plots (06_dmr_pca.R) for every finished comparison in one or more runs.

Writes, in --out_dir (default results/dmrichr/pca/<runs>/):
  pca_manifest.tsv   one row per comparison x region set
  run_pca.sbatch     SLURM array job, one task per row

Plots go where 06_dmr_pca.R puts them: <run>/<comparison>/DMRs/PCA/<comparison>_<DMR|background>_PCA.pdf/.xlsx
With --collect, the job also copies every PDF into <out_dir>/pdf/ to flip through.

Usage (from the repository root):
  python3 scripts/dmrichr/06_make_pca_jobs.py --runs DMR_AZ_04
  python3 scripts/dmrichr/06_make_pca_jobs.py --runs DMR_AZ_04,DMR_AZ_05 --regions both --pcs 3 --collect
  python3 scripts/dmrichr/06_make_pca_jobs.py --runs DMR_AZ_04 --comparisons pre_recruit_vs_instructor --label
  sbatch results/dmrichr/pca/DMR_AZ_04/run_pca.sbatch
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
METH = {'dmr': ('DMR', 'DMRs/DMR_individual_smoothed_methylation.txt'),
        'background': ('background', 'DMRs/background_region_individual_smoothed_methylation.txt')}


def absol(p):
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', required=True, help='comma list of run IDs, e.g. DMR_AZ_04,DMR_AZ_05')
    ap.add_argument('--run_root', default='results/dmrichr')
    ap.add_argument('--comparisons', default=None, help='comma list to restrict comparisons (default: all finished)')
    ap.add_argument('--out_dir', default=None, help='default <run_root>/pca/<runs joined by __>')
    ap.add_argument('--regions', default='dmr', choices=['dmr', 'background', 'both'],
                    help='dmr = DMRs only; background = all tested regions; both [default dmr]')
    ap.add_argument('--shape', default=None, help='sample_info column for point shape, e.g. Sex')
    ap.add_argument('--colour', default=None, help='sample_info column to colour by (default: test covariate)')
    ap.add_argument('--label', action='store_true', help='label points with sample names')
    ap.add_argument('--show_p', action='store_true', help='show PC group-test p-values in the subtitle')
    ap.add_argument('--pcs', type=int, default=2, choices=[2, 3])
    ap.add_argument('--no_scale', action='store_true', help='do not scale regions to unit variance')
    ap.add_argument('--collect', action='store_true', help='also copy all PDFs into <out_dir>/pdf/')
    ap.add_argument('--force', action='store_true', help='redo comparisons that already have a PCA plot')
    ap.add_argument('--pixi_manifest', default='reproducibility/dmrichr/pixi/pixi.toml')
    ap.add_argument('--mem', default=None, help='default 16G (dmr) or 64G (background/both)')
    ap.add_argument('--time', default='04:00:00')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()]
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    regions = ['dmr', 'background'] if args.regions == 'both' else [args.regions]
    root = absol(args.run_root)
    out = absol(args.out_dir) if args.out_dir else root / 'pca' / '__'.join(runs)
    rs = REPO_ROOT / 'scripts/dmrichr/06_dmr_pca.R'
    pixi = absol(args.pixi_manifest)
    for f in [rs, pixi]:
        if not f.exists():
            print(f"WARNING: not found here: {f} (fine if you run this off Hive)")

    rows, skipped = [], []
    for run in runs:
        rd = root / run
        if not rd.is_dir():
            sys.exit(f"Run folder not found: {rd}")
        for comp in sorted(p.name for p in rd.iterdir() if p.is_dir()):
            cd = rd / comp
            if (keep and comp not in keep) or not (cd / 'sample_info.xlsx').exists():
                continue
            for reg in regions:
                tag, mf = METH[reg]
                if not (cd / mf).exists():
                    skipped.append(f"{run}/{comp} [{reg}] (no {Path(mf).name})")
                elif (cd / 'DMRs/PCA' / f"{comp}_{tag}_PCA.pdf").exists() and not args.force:
                    skipped.append(f"{run}/{comp} [{reg}] (done; --force to redo)")
                else:
                    rows.append((run, comp, reg, str(cd)))
    if not rows:
        sys.exit("Nothing to do.\n  " + '\n  '.join(skipped))

    (out / 'logs').mkdir(parents=True, exist_ok=True)
    man = out / 'pca_manifest.tsv'
    with open(man, 'w') as f:
        f.write('task\trun\tcomparison\tregions\tcomparison_dir\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(r) + '\n')

    q = lambda p: shlex.quote(str(p))
    opts = f"--pcs {args.pcs}"
    opts += f" --shape {q(args.shape)}" if args.shape else ''
    opts += f" --colour {q(args.colour)}" if args.colour else ''
    opts += ' --label' if args.label else ''
    opts += ' --show_p' if args.show_p else ''
    opts += ' --scale FALSE' if args.no_scale else ''
    mem = args.mem or ('16G' if args.regions == 'dmr' else '64G')
    collect = (f'mkdir -p {q(out / "pdf")}\n'
               f'TAG=$([ "$REG" = dmr ] && echo DMR || echo background)\n'
               f'cp "$CDIR/DMRs/PCA/${{COMP}}_${{TAG}}_PCA.pdf" {q(out / "pdf")}/"${{RUN}}__${{COMP}}_${{TAG}}_PCA.pdf"\n'
               ) if args.collect else ''
    (out / 'run_pca.sbatch').write_text(f"""#!/bin/bash
#SBATCH --job-name=dmr_pca
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={out}/logs/pca_%A_%a.log
# Written by 06_make_pca_jobs.py: one array task per row of pca_manifest.tsv
set -euo pipefail

LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' {q(man)})
[ -n "$LINE" ] || {{ echo "No task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK RUN COMP REG CDIR <<< "$LINE"
echo "Task $TASK: $RUN / $COMP [$REG]  ($(date))"

eval "$(pixi shell-hook --manifest-path {q(pixi)})"
Rscript {q(rs)} --comparison_dir "$CDIR" --regions "$REG" {opts}
{collect}echo "Done ($(date))"
""")

    print(f"{len(rows)} PCA tasks:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}/{r[1]}  [{r[2]}]")
    if skipped:
        print("Skipped:\n  " + '\n  '.join(skipped))
    rel = os.path.relpath(out, REPO_ROOT)
    print(f"\nWrote {rel}/pca_manifest.tsv and run_pca.sbatch\nSubmit with:\n  sbatch {rel}/run_pca.sbatch")


if __name__ == '__main__':
    main()