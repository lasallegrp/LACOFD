#!/usr/bin/env python3
"""
Step 3 (optional, recommended before the first DMRichR run): check that the two samples
of each pre/post pair come from a person of the same genetic sex.

The sample sheet can only show that labels are internally consistent. This check uses the
sequencing data itself: for every sample it sums read coverage on chrX, chrY and the
autosomes from the CpG cytosine report, infers sex from the chrY/chrX coverage ratio, and
flags any pair whose two samples disagree (a sign of a sample swap or a wrong pair number).

A sex match does not prove identity (about half of random pairs also match); a mismatch
proves a problem. For a full identity check, compare SNP genotypes between pair members.

Usage (run from the repository root, ideally inside a SLURM job; reads every report once):

  python3 scripts/dmrichr/03_check_pair_sex.py \
      --config scripts/dmrichr/config/DMR_AZ_01.yaml --threads 16

  sbatch --account=lasallegrp --partition=high --cpus-per-task=16 --mem=8G --time=06:00:00 \
      --wrap "python3 scripts/dmrichr/03_check_pair_sex.py \
              --config scripts/dmrichr/config/DMR_AZ_01.yaml --threads 16"

Output: <sample_table folder>/<RUN_ID>_pair_sex_check.tsv, a summary on screen, and genetic sex per
sample (genetic_sex_file in the config, default data/processed/wgbs/genetic_sex.tsv), which the run
configs can use as the Sex covariate. Run with --all so every sequenced sample gets a value.
Exit status is 1 if any pair disagrees.
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import yaml

def find_repo_root():
    p = Path(__file__).resolve().parent
    for q in [p, *p.parents]:
        if (q / '.git').exists():
            return q
    return Path(__file__).resolve().parents[2]


REPO_ROOT = find_repo_root()
AWK = r'''{ c=$4+$5; k=$1; sub(/^chr/,"",k);
  if (k=="X") {x+=c; nx++} else if (k=="Y") {y+=c; ny++} else if (k ~ /^[0-9]+$/) {a+=c; na++} }
  END { printf "%d\t%d\t%d\t%d\t%d\t%d\n", x, nx, y, ny, a, na }'''


def rpath(p):
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def coverage(path):
    out = subprocess.run(f"zcat '{path}' | awk -F'\\t' '{AWK}'", shell=True, capture_output=True, text=True, check=True)
    x, nx, y, ny, a, na = map(int, out.stdout.split())
    return {'covX': x, 'cpgX': nx, 'covY': y, 'cpgY': ny, 'covAuto': a, 'cpgAuto': na}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--all', action='store_true', help='every sequenced sample (also writes genetic sex for all of them), not only kept paired ones')
    args = ap.parse_args()

    cfg = yaml.safe_load(open(rpath(args.config)))
    t = pd.read_csv(rpath(cfg['sample_table']), sep='\t', dtype={'Pair': str})
    if not args.all:
        t = t[t.keep & t.Pair.notna() & t.timepoint_group.isin(['pre', 'post'])]

    with ThreadPoolExecutor(args.threads) as ex:
        res = list(ex.map(lambda p: coverage(rpath(p)), t.cytosine_report))
    c = pd.concat([t[['Name', 'Career', 'timepoint_group', 'Pair']].reset_index(drop=True), pd.DataFrame(res)], axis=1)
    # Mean coverage per CpG on each chromosome class, then Y relative to X
    c['meanX'] = c.covX / c.cpgX.clip(lower=1)
    c['meanY'] = c.covY / c.cpgY.clip(lower=1)
    c['meanAuto'] = c.covAuto / c.cpgAuto.clip(lower=1)
    c['Y_over_X'] = (c.meanY / c.meanX.where(c.meanX > 0)).round(4)
    c['X_over_auto'] = (c.meanX / c.meanAuto.where(c.meanAuto > 0)).round(3)

    # Split at the largest gap in the sorted Y/X ratios (two clusters: female low, male high)
    r = c.Y_over_X.dropna().sort_values()
    gaps = r.diff()
    if len(r) > 3 and gaps.max() > 0 and r.iloc[-1] > 2 * r.iloc[0] + 1e-9:
        cut = (r.loc[gaps.idxmax()] + r.iloc[r.index.get_loc(gaps.idxmax()) - 1]) / 2
        c['inferred_sex'] = c.Y_over_X.map(lambda v: 'M' if v > cut else 'F')
        print(f"Y/X coverage ratio split at {cut:.4f} (largest gap {gaps.max():.4f})")
    else:
        cut = None
        c['inferred_sex'] = 'unclear'
        print("No clear separation in Y/X ratios: all samples may be the same sex. Pair check is uninformative.")

    sex_out = rpath(cfg.get('genetic_sex_file', 'data/processed/wgbs/genetic_sex.tsv'))
    sex_out.parent.mkdir(parents=True, exist_ok=True)
    c.assign(genetic_sex=c.inferred_sex.map({'M': 'Male', 'F': 'Female'}))[
        ['Name', 'genetic_sex', 'Y_over_X', 'X_over_auto']].to_csv(sex_out, sep='\t', index=False)
    print(f"Wrote genetic sex for {len(c)} samples to {os.path.relpath(sex_out, REPO_ROOT)}")
    c = c[c.Pair.notna() & c.timepoint_group.isin(['pre', 'post'])]
    per = c.groupby('Pair').agg(Career=('Career', 'first'), samples=('Name', lambda x: ','.join(x)),
                                sexes=('inferred_sex', lambda x: ','.join(x)), n=('Name', 'size'))
    per['status'] = per.apply(lambda r: 'single sample' if r.n < 2 else
                              ('MISMATCH' if len(set(r.sexes.split(','))) > 1 else 'match'), axis=1)
    c = c.merge(per[['status']], left_on='Pair', right_index=True)

    out = rpath(cfg['sample_table']).parent / f"{cfg['run_id']}_pair_sex_check.tsv"
    c.drop(columns=['meanX', 'meanY', 'meanAuto']).to_csv(out, sep='\t', index=False)
    print(f"\nInferred sex: {c.inferred_sex.value_counts().to_dict()}")
    print(per.status.value_counts().to_string())
    bad = per[per.status == 'MISMATCH']
    if len(bad):
        print('\nPairs whose samples disagree on genetic sex (check for swaps):')
        print(bad.to_string())
    print(f"\nWrote {os.path.relpath(out, REPO_ROOT)}")
    sys.exit(1 if len(bad) else 0)


if __name__ == '__main__':
    main()
