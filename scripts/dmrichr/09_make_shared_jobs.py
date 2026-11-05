#!/usr/bin/env python3
"""
Step 9: enrichment of what is SHARED between comparisons (after 05_dmrichr_overlap.R).

--mode gene  (from gene/<group>_overlap_gene.xlsx)
  09_shared_genes.R: genes with a DMR in all comparisons of the group; background = genes
  testable in ALL of them. 07_enrichment_plot.R --genes ... --gene_assignment annotation
  -> clusterProfiler + Enrichr (no GREAT: it needs regions).
  Output: <overlap_dir>/<family>/<group>/gene/shared_gene_enrichment/

--mode coord (from coord/<group>_overlap_coord.xlsx)
  09_shared_loci.R: loci where DMRs of all comparisons physically overlap; background =
  territory tested in ALL of them (intersection of backgroundRegions.bed).
  07_enrichment_plot.R --dmrs ... --background_regions ... --gene_assignment great
  -> GREAT, plus clusterProfiler + Enrichr on GREAT-assigned genes.
  Output: <overlap_dir>/<family>/<group>/coord/shared_loci_enrichment/

--mode both (default) = one array task per group per mode.

Usage (from the repository root):
  python3 scripts/dmrichr/09_make_shared_jobs.py --overlap_dir results/dmrichr/overlap/DMR_AZ_04
  python3 scripts/dmrichr/09_make_shared_jobs.py --overlap_dir results/dmrichr/overlap/DMR_AZ_04 \\
      --groups custom_AZ04_career_by_timepoint --mode coord --concordant --split_direction
  sbatch results/dmrichr/overlap/DMR_AZ_04/shared_enrichment/run_shared_enrichment.sbatch
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--overlap_dir', required=True, help='folder from 04/05, holding manifest.tsv')
    ap.add_argument('--groups', default=None, help='comma list of group folder names (default: all)')
    ap.add_argument('--mode', default='both', choices=['gene', 'coord', 'both'])
    ap.add_argument('--shared', default='all', choices=['all', 'any2'],
                    help='all = gene in every comparison of the group; any2 = in at least two [default all]')
    ap.add_argument('--concordant', action='store_true', help='same direction in every comparison only')
    ap.add_argument('--split_direction', action='store_true', help='also run shared-hyper and shared-hypo lists')
    ap.add_argument('--min_items', '--min_genes', type=int, default=20,
                    help='lists with fewer genes/loci are not tested [default 20]')
    ap.add_argument('--gene_set_cache', default=GENE_SETS)
    ap.add_argument('--great_gene_table', default='', help='local GREATv4.genes.hg38.tsv(.gz), needed for GREAT on compute nodes')
    ap.add_argument('--pixi_manifest', default='reproducibility/dmrichr/pixi/pixi.toml')
    ap.add_argument('--enrich_args', default='', help='extra options for 07_enrichment_plot.R, quoted')
    ap.add_argument('--force', action='store_true', help='redo groups that already have results')
    ap.add_argument('--mem', default='32G')
    ap.add_argument('--time', default='12:00:00')
    args = ap.parse_args()

    od = absol(args.overlap_dir)
    man = od / 'manifest.tsv'
    if not man.exists():
        sys.exit(f"No manifest.tsv in {od} (run 04_make_overlap_sets.py and 05 first)")
    keep = set(g.strip() for g in args.groups.split(',')) if args.groups else None
    prep = {'gene': REPO_ROOT / 'scripts/dmrichr/09_shared_genes.R',
            'coord': REPO_ROOT / 'scripts/dmrichr/09_shared_loci.R'}
    enrich = REPO_ROOT / 'scripts/dmrichr/07_enrichment_plot.R'
    helper = REPO_ROOT / 'scripts/dmrichr/helper.R'
    pixi = absol(args.pixi_manifest)
    for f in [*prep.values(), enrich, helper, pixi, Path(args.gene_set_cache)]:
        if not f.exists():
            print(f"WARNING: not found here: {f} (fine if you run this off Hive)")

    rows, done = [], []
    modes = ['gene', 'coord'] if args.mode == 'both' else [args.mode]
    sub = {'gene': 'shared_gene_enrichment', 'coord': 'shared_loci_enrichment'}
    for mode in modes:
        for x in sorted(od.glob(f'*/*/{mode}/*_overlap_{mode}.xlsx')):
            group = x.parent.parent.name
            if keep and group not in keep:
                continue
            out = x.parent / sub[mode]
            if (out / 'lists.tsv').exists() and not args.force:
                done.append(f"{group} ({mode})")
                continue
            rows.append((group, mode, str(x), str(out)))
    if not rows:
        sys.exit("Nothing to do." + (f" Already done (--force to redo): {', '.join(done)}" if done else
                                     f" No */*/<gene|coord>/*_overlap_*.xlsx under {od} for --mode {args.mode}"))

    jd = od / 'shared_enrichment'
    (jd / 'logs').mkdir(parents=True, exist_ok=True)
    tasks = jd / 'tasks.tsv'
    with open(tasks, 'w') as f:
        f.write('task\tgroup\tmode\toverlap_xlsx\tout_dir\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(r) + '\n')

    prep_opts = f"--shared {args.shared}"
    prep_opts += ' --concordant' if args.concordant else ''
    prep_opts += ' --split_direction' if args.split_direction else ''
    q = lambda p: shlex.quote(str(p))
    sb = f"""#!/bin/bash
#SBATCH --job-name=shared_enrich
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=2
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={jd}/logs/shared_%A_%a.log
# Written by 09_make_shared_jobs.py: one array task per overlap group in tasks.tsv
set -uo pipefail

LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' {q(tasks)})
[ -n "$LINE" ] || {{ echo "No task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK GROUP MODE XLSX OUT <<< "$LINE"
echo "Task $TASK: $GROUP, $MODE overlap  ($(date))"
eval "$(pixi shell-hook --manifest-path {q(pixi)})"
export GREAT_GENE_TABLE={q(args.great_gene_table) if args.great_gene_table else '"${GREAT_GENE_TABLE:-}"'}

# 1. shared genes/loci + background testable in ALL comparisons of the group
if [ "$MODE" = gene ]; then
  Rscript {q(prep['gene'])} --overlap_xlsx "$XLSX" --manifest {q(man)} --out_dir "$OUT" {prep_opts} --min_genes {args.min_items} || exit 1
else
  Rscript {q(prep['coord'])} --overlap_xlsx "$XLSX" --manifest {q(man)} --out_dir "$OUT" {prep_opts} --min_items {args.min_items} || exit 1
fi

# 2. enrichment of each list that is big enough
status=0
while IFS=$'\\t' read -r LIST N RUN; do
  [ "$LIST" = "list" ] && continue
  if [ "$RUN" != "TRUE" ]; then echo "$LIST: $N (< {args.min_items}) -> skipped"; continue; fi
  echo "$LIST: $N -> $OUT/$LIST"
  mkdir -p "$OUT/$LIST"
  if [ "$MODE" = gene ]; then
    INPUT=(--genes "$OUT/${{LIST}}_genes.txt" --background_genes "$OUT/background_genes.txt" --gene_assignment annotation)
  else
    INPUT=(--dmrs "$OUT/${{LIST}}_loci.bed" --background_regions "$OUT/background_regions.bed" --gene_assignment great)
  fi
  Rscript {q(enrich)} "${{INPUT[@]}}" \\
    --gene_set_cache {q(args.gene_set_cache)} \\
    --helper_file {q(helper)} \\
    --outdir "$OUT/$LIST" \\
    --prefix "${{GROUP}}_${{MODE}}_${{LIST}}" {args.enrich_args} > "$OUT/$LIST/enrichment.log" 2>&1 \\
    || {{ echo "  FAILED, see $OUT/$LIST/enrichment.log"; tail -n 5 "$OUT/$LIST/enrichment.log"; status=1; }}
done < "$OUT/lists.tsv"
echo "Done $GROUP ($(date))"
exit $status
"""
    (jd / 'run_shared_enrichment.sbatch').write_text(sb)

    print(f"{len(rows)} tasks:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}  [{r[1]}]")
    if done:
        print("Already done, skipped (--force to redo): " + ', '.join(done))
    rel = os.path.relpath(jd, REPO_ROOT)
    print(f"\nSelection: genes/loci in {'every comparison' if args.shared == 'all' else '>= 2 comparisons'}"
          f"{', same direction' if args.concordant else ''}{'; + hyper/hypo lists' if args.split_direction else ''}"
          f"; background = genes/territory testable in all comparisons of the group")
    print(f"Wrote {rel}/tasks.tsv and run_shared_enrichment.sbatch\nSubmit with:\n  sbatch {rel}/run_shared_enrichment.sbatch")


if __name__ == '__main__':
    main()