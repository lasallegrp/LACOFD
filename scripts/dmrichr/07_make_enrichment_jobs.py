#!/usr/bin/env python3
"""
Step 7: enrichment (07_enrichment_plot.R + helper.R) for every finished comparison in one or more runs.

Writes, in --out_dir (default results/dmrichr/enrichment/<runs>/):
  enrichment_manifest.tsv   one row per comparison to run
  run_enrichment.sbatch     SLURM array job, one task per comparison

Each comparison: --dmrs DMRs/DMRs_annotated.xlsx --background_regions DMRs/backgroundRegions.bed
-> <run>/<comparison>/DMRs/enrichment/ (GREAT + clusterProfiler + Enrichr; log in enrichment.log)

Usage (from the repository root):
  python3 scripts/dmrichr/07_make_enrichment_jobs.py --runs DMR_AZ_04
  python3 scripts/dmrichr/07_make_enrichment_jobs.py --runs DMR_AZ_04,DMR_AZ_05 --min_dmrs 100
  python3 scripts/dmrichr/07_make_enrichment_jobs.py --runs DMR_AZ_04 --comparisons pre_recruit_vs_instructor
  sbatch results/dmrichr/enrichment/DMR_AZ_04/run_enrichment.sbatch
python3 scripts/dmrichr/07_make_enrichment_jobs.py --runs DMR_AZ_04 --force \
  --great_gene_table /quobyte/lasallegrp/projects/LACOFD/reference/GREATv4.genes.hg38.tsv.gz
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
GENE_SETS = '/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reference/gene_sets_hg38.rds'


def absol(p):
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def count_bed(f):
    with open(f) as fh:
        return sum(1 for l in fh if l.strip() and not l.startswith(('#', 'track', 'browser')))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', required=True, help='comma list of run IDs, e.g. DMR_AZ_04,DMR_AZ_05')
    ap.add_argument('--run_root', default='results/dmrichr')
    ap.add_argument('--comparisons', default=None, help='comma list to restrict comparisons (default: all finished)')
    ap.add_argument('--out_dir', default=None, help='default <run_root>/enrichment/<runs joined by __>')
    ap.add_argument('--min_dmrs', type=int, default=50, help='skip comparisons with fewer DMRs [default 50]')
    ap.add_argument('--gene_set_cache', default=GENE_SETS)
    ap.add_argument('--pixi_manifest', default='reproducibility/dmrichr/pixi/pixi.toml')
    ap.add_argument('--great_gene_table', default='', help='local GREATv4.genes.hg38.tsv.gz, if rGREAT cannot download it')
    ap.add_argument('--enrich_args', default='', help='extra options for 07_enrichment_plot.R, quoted, e.g. "--mode presentation"')
    ap.add_argument('--force', action='store_true', help='redo comparisons that already have DMRs/enrichment/')
    ap.add_argument('--mem', default='48G')
    ap.add_argument('--cpus', type=int, default=4)
    ap.add_argument('--time', default='24:00:00')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()]
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    root = absol(args.run_root)
    out = absol(args.out_dir) if args.out_dir else root / 'enrichment' / '__'.join(runs)
    rs = REPO_ROOT / 'scripts/dmrichr/07_enrichment_plot.R'
    helper = REPO_ROOT / 'scripts/dmrichr/helper.R'
    pixi = absol(args.pixi_manifest)
    for f in [rs, helper, pixi, Path(args.gene_set_cache)]:
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
            dmrs, bg, bed = cd / 'DMRs/DMRs_annotated.xlsx', cd / 'DMRs/backgroundRegions.bed', cd / 'DMRs/DMRs.bed'
            if not (dmrs.exists() and bg.exists() and bed.exists()):
                skipped.append(f"{run}/{comp} (not finished)")
                continue
            n = count_bed(bed)
            if n < args.min_dmrs:
                skipped.append(f"{run}/{comp} ({n} DMRs < {args.min_dmrs})")
                continue
            od = cd / 'DMRs/enrichment'
            if od.exists() and any(od.iterdir()) and not args.force:
                skipped.append(f"{run}/{comp} (done; --force to redo)")
                continue
            rows.append((run, comp, n, str(dmrs), str(bg), str(od)))
    if not rows:
        sys.exit("Nothing to do.\n  " + '\n  '.join(skipped))

    (out / 'logs').mkdir(parents=True, exist_ok=True)
    man = out / 'enrichment_manifest.tsv'
    with open(man, 'w') as f:
        f.write('task\trun\tcomparison\tn_dmrs\tdmrs\tbackground\toutdir\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(map(str, r)) + '\n')

    q = lambda p: shlex.quote(str(p))
    (out / 'run_enrichment.sbatch').write_text(f"""#!/bin/bash
#SBATCH --job-name=dmr_enrich
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task={args.cpus}
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={out}/logs/enrich_%A_%a.log
# Written by 07_make_enrichment_jobs.py: one array task per comparison in enrichment_manifest.tsv
set -euo pipefail

LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' {q(man)})
[ -n "$LINE" ] || {{ echo "No task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK RUN COMP N DMRS BG OUT <<< "$LINE"
echo "Task $TASK: $RUN / $COMP, $N DMRs  ($(date))"

eval "$(pixi shell-hook --manifest-path {q(pixi)})"
export GREAT_GENE_TABLE={q(args.great_gene_table)}
mkdir -p "$OUT"
Rscript {q(rs)} \\
  --dmrs "$DMRS" \\
  --background_regions "$BG" \\
  --gene_set_cache {q(args.gene_set_cache)} \\
  --helper_file {q(helper)} \\
  --outdir "$OUT" \\
  --prefix "$COMP" {args.enrich_args} 2>&1 | tee "$OUT/enrichment.log"
echo "Done: $OUT ($(date))"
""")

    print(f"{len(rows)} comparisons:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}/{r[1]}  ({r[2]:,} DMRs)")
    if skipped:
        print("Skipped:\n  " + '\n  '.join(skipped))
    rel = os.path.relpath(out, REPO_ROOT)
    print(f"\nWrote {rel}/enrichment_manifest.tsv and run_enrichment.sbatch\nSubmit with:\n  sbatch {rel}/run_enrichment.sbatch")


if __name__ == '__main__':
    main()
