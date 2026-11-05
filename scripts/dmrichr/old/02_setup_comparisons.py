#!/usr/bin/env python3
"""
Step 2 of 3: set up one folder per DMRichR comparison.

For each comparison in the run config, this script:
  1) selects samples from the sample table made by 01_build_sample_table.py
     (kept samples only; for paired comparisons, complete pre/post pairs only);
  2) hard-links each sample's CpG cytosine report into the comparison folder
     (no extra disk space; the epigenerator file is never moved or renamed);
  3) writes sample_info.xlsx in the form DMRichR expects (first column "Name"),
     with group labels ordered so DMRichR reports "test vs reference";
  4) writes run_dmrichr.sbatch, a sample list and a snapshot of the settings.

Layout created (run ID and comparison names come from the config):

  <run_root>/<RUN_ID>/            (run_root from the config, e.g. results/dmrichr)
  ├── config_snapshot.yaml
  ├── scripts_version.txt
  ├── comparisons_summary.tsv
  └── <comparison>/
      ├── <Name>_*.CpG_report.txt.gz   hard links, one per sample
      ├── sample_info.xlsx             DMRichR design matrix
      ├── sample_list.tsv              same samples with QC and labels, for the record
      ├── run_dmrichr.sbatch
      └── logs/

Usage (run from the repository root):

  # Dry run (default): prints the plan, changes nothing
  python3 scripts/dmrichr/02_setup_comparisons.py \
      --config scripts/dmrichr/config/DMR_AZ_01.yaml

  # Create folders, links and files
  python3 scripts/dmrichr/02_setup_comparisons.py \
      --config scripts/dmrichr/config/DMR_AZ_01.yaml --apply

  # One comparison only, then submit it
  python3 scripts/dmrichr/02_setup_comparisons.py \
      --config scripts/dmrichr/config/DMR_AZ_01.yaml \
      --only recruit_post_vs_pre --apply --submit

Safety:
  - Re-running is safe: links that already point at the right file are skipped.
  - A comparison folder that already holds DMRichR output (RData/ or DMRs/) is not changed
    unless --force is given. Change the run ID in the config to start a new run instead.
  - Cytosine reports in a folder that are not in the current selection are reported;
    --prune removes those links (the epigenerator originals are untouched).
"""

import argparse
import os
import shutil
import subprocess
import sys
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
SCRIPTS_DIR = Path(__file__).resolve().parent
GROUP_LABEL = {'pre': 'Pre', 'post': 'Post', 'baseline_unpaired': 'Baseline'}


def rpath(p):
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def ordered_labels(ref, test):
    """DMRichR uses the alphabetically first level as the reference. Prefix if needed."""
    if ref.lower() < test.lower():
        return ref, test
    return f'1_{ref}', f'2_{test}'


def select(t, comp):
    s = t[t.keep & t.Career.isin(comp['career']) & t.timepoint_group.isin(comp['timepoint_group'])].copy()
    if comp.get('paired'):
        s = s[s.pair_complete_kept & s.timepoint_group.isin(['pre', 'post'])]
    return s


def design(s, comp):
    cov = comp['test_covariate']
    if cov == 'Timepoint':
        raw = s.timepoint_group.map(lambda g: 'pre' if g in ('pre', 'baseline_unpaired') else g)
    elif cov == 'Career':
        raw = s.Career
    else:
        sys.exit(f"test_covariate must be Timepoint or Career, got {cov}")
    ref, test = str(comp['reference_level']), str(comp['test_level'])
    if set(raw) - {ref, test}:
        sys.exit(f"{cov} has values {sorted(set(raw))}, expected only {ref} and {test}")
    ref_lab, test_lab = ordered_labels(GROUP_LABEL.get(ref, ref), GROUP_LABEL.get(test, test))
    d = pd.DataFrame({'Name': s.Name.values, cov: raw.map({ref: ref_lab, test: test_lab}).values})
    for a in comp.get('adjust_covariates') or []:
        if a == 'Pair':
            d['Pair'] = s.Pair.values
        elif a in s.columns:
            d[a] = s[a].values
        else:
            sys.exit(f"adjust covariate {a} is not a column of the sample table")
    return d, ref_lab, test_lab


def check_design(name, d, cov, ref_lab, test_lab, paired):
    problems = []
    n_ref, n_test = (d[cov] == ref_lab).sum(), (d[cov] == test_lab).sum()
    if min(n_ref, n_test) < 2:
        problems.append(f"{name}: fewer than 2 samples in a group ({ref_lab} {n_ref}, {test_lab} {n_test})")
    if paired:
        per = d.groupby('Pair')[cov].agg(lambda x: tuple(sorted(x)))
        badp = per[per != tuple(sorted((ref_lab, test_lab)))]
        if len(badp):
            problems.append(f"{name}: pairs without exactly one {ref_lab} and one {test_lab}: {list(badp.index)}")
    return problems, n_ref, n_test


def sbatch_text(cfg, comp_name, comp_dir, cov, adjust):
    dm, sl = cfg['dmrichr'], cfg['slurm']
    dm_script = str(dm.get('dm_r_script', 'auto'))
    if dm_script == 'auto':
        find = ('DMR_SCRIPT=$(Rscript -e \'cat(system.file("exec", "DM.R", package = "DMRichR"))\')\n'
                'if [ -z "$DMR_SCRIPT" ]; then echo "DM.R not found: is DMRichR installed in this environment?" >&2; exit 1; fi')
    else:
        find = f'DMR_SCRIPT={rpath(dm_script)}'
    adj = f"\\\n  --adjustCovariate '{';'.join(adjust)}' " if adjust else ''
    return f"""#!/bin/bash
#SBATCH --job-name={cfg['run_id']}_{comp_name}
#SBATCH --account={sl['account']}
#SBATCH --partition={sl['partition']}
#SBATCH --time={sl['time']}
#SBATCH --mem={sl['mem']}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={dm['cores']}
#SBATCH --chdir={comp_dir}
#SBATCH --output={comp_dir}/logs/%x_%j.out
#SBATCH --error={comp_dir}/logs/%x_%j.err

# Generated by scripts/dmrichr/02_setup_comparisons.py. Edit the run config
# and re-run the script rather than editing this file.

set -euo pipefail
{sl['env_setup'].rstrip()}

cd {comp_dir}
{find}
{{ echo "DMRichR $(Rscript -e 'cat(as.character(packageVersion("DMRichR")))')"; echo "DM.R $DMR_SCRIPT"; \\
  echo "R $(Rscript -e 'cat(R.version.string)')"; }} > dmrichr_version.txt

call="Rscript --vanilla $DMR_SCRIPT \\
  --genome {dm['genome']} \\
  --coverage {dm['coverage']} \\
  --perGroup '{dm['perGroup']}' \\
  --minCpGs {dm['minCpGs']} \\
  --maxPerms {dm['maxPerms']} \\
  --maxBlockPerms {dm['maxBlockPerms']} \\
  --cutoff '{dm['cutoff']}' \\
  --testCovariate {cov} {adj}\\
  --sexCheck {dm['sexCheck']} \\
  --GOfuncR {dm['GOfuncR']} \\
  --EnsDb {dm['EnsDb']} \\
  --cores {dm['cores']}"

echo "$call"
eval "$call"
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    ap.add_argument('--only', nargs='+', help='comparison names to set up (default: all)')
    ap.add_argument('--apply', action='store_true', help='create folders, links and files')
    ap.add_argument('--submit', action='store_true', help='sbatch each comparison after --apply')
    ap.add_argument('--prune', action='store_true', help='remove links to samples no longer selected')
    ap.add_argument('--force', action='store_true', help='allow changes to folders that already hold DMRichR output')
    args = ap.parse_args()

    cfg = yaml.safe_load(open(rpath(args.config)))
    t = pd.read_csv(rpath(cfg['sample_table']), sep='\t', dtype={'Pair': str})
    run_dir = rpath(cfg.get('run_root', 'results/dmrichr')) / cfg['run_id']
    comps = cfg['comparisons']
    if args.only:
        missing = set(args.only) - set(comps)
        if missing:
            sys.exit(f"Unknown comparison(s): {sorted(missing)}. Known: {list(comps)}")
        comps = {k: comps[k] for k in args.only}

    dm_script = str(cfg['dmrichr'].get('dm_r_script', 'auto'))
    if dm_script != 'auto' and not rpath(dm_script).exists():
        print(f"WARNING: DM.R not found at {rpath(dm_script)}")

    summary, all_problems = [], []
    print(f"{'DRY RUN' if not args.apply else 'APPLY'}: run {cfg['run_id']} -> {os.path.relpath(run_dir, REPO_ROOT)}\n")
    for name, comp in comps.items():
        s = select(t, comp)
        d, ref_lab, test_lab = design(s, comp)
        cov = comp['test_covariate']
        problems, n_ref, n_test = check_design(name, d, cov, ref_lab, test_lab, comp.get('paired'))
        all_problems += problems
        comp_dir = run_dir / name
        print(f"== {name}: {comp['description']}")
        print(f"   {cov}: {test_lab} n={n_test} vs {ref_lab} n={n_ref} (reference)"
              f"{'; pairs=' + str(d.Pair.nunique()) if comp.get('paired') else ''}"
              f"; adjust={comp.get('adjust_covariates') or 'none'}")
        for p in problems:
            print(f"   ERROR: {p}")

        existing = {p.name for p in comp_dir.glob('*CpG_report.txt.gz')} if comp_dir.exists() else set()
        wanted = {os.path.basename(x): rpath(x) for x in s.cytosine_report}
        stale = sorted(existing - set(wanted))
        to_link = sorted(set(wanted) - existing)
        done_output = comp_dir.exists() and any((comp_dir / x).exists() for x in ('RData', 'DMRs'))
        print(f"   links: {len(wanted)} wanted, {len(existing & set(wanted))} present, {len(to_link)} to add, {len(stale)} stale")
        if stale:
            print(f"   stale (not in selection): {', '.join(stale)}{' -> removing' if args.prune else ' (use --prune to remove)'}")
        if done_output and not args.force:
            print("   SKIP: folder already holds DMRichR output (RData/ or DMRs/). Use a new run ID, or --force.")
            all_problems.append(f"{name}: has DMRichR output; not modified")
            continue
        summary.append({'comparison': name, 'description': comp['description'], 'test_covariate': cov,
                        'test_label': test_lab, 'n_test': n_test, 'reference_label': ref_lab, 'n_reference': n_ref,
                        'pairs': d.Pair.nunique() if comp.get('paired') else '',
                        'adjust': ';'.join(comp.get('adjust_covariates') or [])})
        if not args.apply or problems:
            print()
            continue

        (comp_dir / 'logs').mkdir(parents=True, exist_ok=True)
        for fname, src in wanted.items():
            dst = comp_dir / fname
            if dst.exists():
                if not os.path.samefile(src, dst):
                    sys.exit(f"{dst} exists but is not the same file as {src}; remove it and re-run")
                continue
            os.link(src, dst)
        if args.prune:
            for fname in stale:
                (comp_dir / fname).unlink()
        d.to_excel(comp_dir / 'sample_info.xlsx', index=False)
        s.to_csv(comp_dir / 'sample_list.tsv', sep='\t', index=False)
        (comp_dir / 'run_dmrichr.sbatch').write_text(sbatch_text(cfg, name, comp_dir, cov, comp.get('adjust_covariates') or []))
        (comp_dir / 'comparison.yaml').write_text(yaml.safe_dump({'run_id': cfg['run_id'], 'comparison': name, **comp,
                                                                  'labels': {'reference': ref_lab, 'test': test_lab}},
                                                                 sort_keys=False))
        print(f"   wrote {os.path.relpath(comp_dir, REPO_ROOT)}/ (links, sample_info.xlsx, sample_list.tsv, run_dmrichr.sbatch)")
        if args.submit:
            r = subprocess.run(['sbatch', str(comp_dir / 'run_dmrichr.sbatch')], capture_output=True, text=True)
            print(f"   sbatch: {r.stdout.strip() or r.stderr.strip()}")
        print()

    if args.apply:
        run_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(rpath(args.config), run_dir / 'config_snapshot.yaml')
        shutil.copy(rpath(cfg['sample_table']), run_dir / 'sample_table_snapshot.tsv')
        ver = subprocess.run(['git', '-C', str(SCRIPTS_DIR), 'log', '-1', '--format=%H %ad', '--', '.'],
                             capture_output=True, text=True).stdout.strip()
        (run_dir / 'scripts_version.txt').write_text((ver or 'uncommitted scripts') + '\n')
        lock = cfg.get('slurm', {}).get('pixi_lock')
        if lock and Path(lock).exists():
            shutil.copy(lock, run_dir / 'pixi.lock')
        elif lock:
            print(f"WARNING: pixi_lock not found: {lock}")
        pd.DataFrame(summary).to_csv(run_dir / 'comparisons_summary.tsv', sep='\t', index=False)

    if summary:
        print(pd.DataFrame(summary)[['comparison', 'test_label', 'n_test', 'reference_label', 'n_reference', 'pairs']]
              .to_string(index=False))
    if all_problems:
        print('\nProblems:\n  ' + '\n  '.join(all_problems))
        sys.exit(1)


if __name__ == '__main__':
    main()
