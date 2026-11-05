#!/usr/bin/env python3
"""
Step 1 of 3: build and validate the sample table for a DMRichR run.

Joins the sample sheet (career status, timepoint, pair) with QC from each epigenerator
run's MultiQC data folder, flags samples, checks that the pre/post labels and pairs are
internally consistent, and decides which samples the run keeps.

Nothing is linked or submitted here. The outputs are a sample table (TSV) and a
validation report (Markdown) at the paths named in the run config.

Usage (run from the repository root):

  python3 scripts/dmrichr/01_build_sample_table.py \
      --config scripts/dmrichr/config/DMR_AZ_01.yaml

Exit status is 1 if any check fails (an "ERROR" line in the report). Warnings are listed
but do not stop the script.

Sample sheet columns (matched case-insensitively, by the start of the header):
  Sequencing Submission Name -> Name        (must match the cytosine report prefix, e.g. LA001)
  Career Status              -> Career      (Recruit / Instructor)
  Timepoint                  -> Timepoint   (pre-exposure / post-exposure / baseline /
                                             not paired - baseline)
  Pair #                     -> Pair        (number shared by a person's pre and post samples)
  Specimen.ID                -> specimen_id (optional; used for the specimen-number check)
  Size-Selected (Y/N)        -> size_selected (optional)
  "...actual volume corrected (ng)" -> input_ng (optional)

Optional inputs (config "inputs"), used only for checks; they never leave the machine:
  voc_results  VOC/PAH results (Pair #, plasma_pre_date, plasma_post_date, VOC_specimenid, ...)
               -> checks each pair's plasma pre date is before its post date; adds plasma_date.
  id_map       FFCCS_ID_Map.csv (record_id, VOC_specimenid; urine specimens)
               -> with voc_results: finds the person behind each pair through its urine specimens;
                  checks a pair is one person and a person is in only one pair; flags unpaired
                  baselines who are also in a pair. Adds record_id.

MultiQC data files read from each run's multiqc_data folder:
  multiqc_bismark_dedup.txt, multiqc_bismark_methextract.txt,
  multiqc_bismark_alignment.txt, multiqc_fastqc.txt
"""

import argparse
import glob
import os
import re
import sys
from datetime import datetime
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
NAME_RE = re.compile(r'^[A-Za-z0-9]+$')          # DMRichR cuts sample names at the first "_"
SAMPLE_PREFIX = re.compile(r'^([A-Za-z0-9]+)')    # MultiQC sample -> Name (text before first "_")
LANE_RE = re.compile(r'_(L\d{3})_')

TIMEPOINT_MAP = {
    'pre-exposure': 'pre',
    'post-exposure': 'post',
    'baseline': 'baseline',                        # resolved below: paired -> pre, else unpaired
    'not paired - baseline': 'baseline_unpaired',
}


def rpath(p):
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def pick(df, startswith, required=True):
    for c in df.columns:
        if str(c).strip().lower().startswith(startswith.lower()):
            return c
    if required:
        sys.exit(f"Sample sheet: no column starting with '{startswith}'. Columns: {list(df.columns)}")
    return None


def pick_contains(df, text):
    for c in df.columns:
        if text.lower() in str(c).lower():
            return c
    return None


def read_mqc(folder, fname, needed):
    f = Path(folder) / fname
    if not f.exists():
        sys.exit(f"Missing MultiQC table: {f}")
    df = pd.read_csv(f, sep='\t')
    miss = [c for c in needed if c not in df.columns]
    if miss:
        sys.exit(f"{f}: missing columns {miss}. Columns present: {list(df.columns)}")
    return df


def find_report(folder, name):
    hits = sorted(set(glob.glob(os.path.join(folder, f'{name}_*CpG_report.txt.gz')) +
                      glob.glob(os.path.join(folder, f'{name}.CpG_report.txt.gz'))))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(rpath(args.config)))
    qc = cfg['qc']
    errors, warnings, notes = [], [], []

    # ------------------------------------------------------------------ sample sheet
    si_path = rpath(cfg['inputs']['sample_info'])
    raw = pd.read_excel(si_path) if si_path.suffix in ('.xlsx', '.xls') else pd.read_csv(si_path, sep=None, engine='python')
    raw = raw.dropna(how='all')
    c_name, c_car, c_tp, c_pair = (pick(raw, 'Sequencing Submission Name'), pick(raw, 'Career'),
                                   pick(raw, 'Timepoint'), pick(raw, 'Pair'))
    c_spec = pick(raw, 'Specimen', required=False)
    c_size = pick(raw, 'Size-Selected', required=False)
    c_ng = pick_contains(raw, 'actual volume corrected')

    si = pd.DataFrame({
        'Name': raw[c_name].astype(str).str.strip(),
        'Career': raw[c_car].astype(str).str.strip(),
        'Timepoint_label': raw[c_tp].astype(str).str.strip(),
        'Pair': raw[c_pair].map(lambda x: None if pd.isna(x) else f'P{int(x)}'),
        'specimen_id': raw[c_spec].astype(str).str.strip() if c_spec else None,
        'size_selected': raw[c_size] if c_size else None,
        'input_ng': raw[c_ng].round(3) if c_ng else None,
    })

    # Names
    dup = si[si.Name.duplicated(keep=False)].Name.unique()
    if len(dup):
        errors.append(f"Duplicated sample names in sample sheet: {', '.join(dup)}")
    bad = si[~si.Name.str.match(NAME_RE)].Name.tolist()
    if bad:
        errors.append(f"Names with characters other than letters/digits (DMRichR cuts at '_'): {bad}")
    names = sorted(si.Name)
    for a in names:
        for b in names:
            if a != b and b.startswith(a):
                errors.append(f"Name {a} is a prefix of {b}; DMRichR's partial matching would confuse them")

    # Career / timepoint values
    badc = sorted(set(si.Career) - {'Recruit', 'Instructor'})
    if badc:
        errors.append(f"Unexpected Career values: {badc}")
    tpl = si.Timepoint_label.str.lower()
    badt = sorted(set(tpl) - set(TIMEPOINT_MAP))
    if badt:
        errors.append(f"Unexpected Timepoint values: {badt}")
    si['timepoint_group'] = tpl.map(TIMEPOINT_MAP)

    # Baselines: paired -> pre ; unpaired label with a pair number -> error
    pb = si[(si.timepoint_group == 'baseline') & si.Pair.notna()]
    si.loc[pb.index, 'timepoint_group'] = 'pre'
    if len(pb):
        warnings.append(f"{len(pb)} samples labelled 'baseline' have a pair number and are used as the pre-exposure "
                        f"sample of their pair: {', '.join(pb.Name)}. Confirm they were collected at the pre-exposure point.")
    ub = si[(si.timepoint_group == 'baseline') & si.Pair.isna()]
    si.loc[ub.index, 'timepoint_group'] = 'baseline_unpaired'
    wrong = si[(si.timepoint_group == 'baseline_unpaired') & si.Pair.notna()]
    if len(wrong):
        errors.append(f"Labelled 'not paired - baseline' but has a pair number: {', '.join(wrong.Name)}")
    nopair = si[si.timepoint_group.isin(['pre', 'post']) & si.Pair.isna() & ~si.index.isin(pb.index)]
    if len(nopair):
        warnings.append(f"Pre/post samples without a pair number: {', '.join(nopair.Name)}")

    # Pair structure
    for pid, g in si[si.Pair.notna()].groupby('Pair'):
        if g.Career.nunique() > 1:
            errors.append(f"Pair {pid} mixes careers: {dict(zip(g.Name, g.Career))}")
        npre, npost = (g.timepoint_group == 'pre').sum(), (g.timepoint_group == 'post').sum()
        if npre > 1 or npost > 1:
            errors.append(f"Pair {pid} has {npre} pre and {npost} post samples: {dict(zip(g.Name, g.Timepoint_label))}")
        # Specimen numbers: within a series (same prefix), pre should precede post
        if c_spec and npre == 1 and npost == 1:
            pre, post = g[g.timepoint_group == 'pre'].iloc[0], g[g.timepoint_group == 'post'].iloc[0]
            mp, mq = re.match(r'^([A-Z]+)-(\d+)$', str(pre.specimen_id)), re.match(r'^([A-Z]+)-(\d+)$', str(post.specimen_id))
            if mp and mq and mp.group(1) == mq.group(1) and int(mp.group(2)) > int(mq.group(2)):
                errors.append(f"Pair {pid}: pre specimen {pre.specimen_id} is numbered after post specimen {post.specimen_id}; "
                              f"check whether the pre/post labels are swapped")

    # ------------------------------------------------------------------ pair dates and identity (optional)
    # VOC/PAH results: one row per urine specimen, with the pair's plasma pre/post dates.
    # ID map: urine specimen -> record_id (person). Together they give the person behind each pair.
    norm = lambda x: str(x).strip().upper().replace(' ', '')
    si['record_id'] = None
    si['plasma_date'] = None
    vcfg, idm = cfg['inputs'].get('voc_results'), cfg['inputs'].get('id_map')
    voc = None
    if vcfg:
        vp = rpath(vcfg['file'] if isinstance(vcfg, dict) else vcfg)
        if not vp.exists():
            errors.append(f"VOC/PAH results file not found: {vp}")
        else:
            voc = pd.read_excel(vp) if vp.suffix in ('.xlsx', '.xls') else pd.read_csv(vp, sep=None, engine='python')
            voc.columns = [str(c).strip() for c in voc.columns]
            need = ['Pair #', 'plasma_pre_date', 'plasma_post_date', 'VOC_specimenid']
            miss = [c for c in need if c not in voc.columns]
            if miss:
                sys.exit(f"{vp}: missing columns {miss}. Columns: {list(voc.columns)}")
            voc = voc.dropna(subset=['Pair #'])
            voc['Pair'] = voc['Pair #'].map(lambda x: f'P{int(float(x))}')
            for c in ('plasma_pre_date', 'plasma_post_date'):
                voc[c] = pd.to_datetime(voc[c], errors='coerce')
            # One pre and one post date per pair
            dates = voc.groupby('Pair').agg(pre=('plasma_pre_date', lambda x: x.dropna().unique()),
                                            post=('plasma_post_date', lambda x: x.dropna().unique()))
            for pid, r in dates.iterrows():
                if len(r.pre) > 1 or len(r.post) > 1:
                    errors.append(f"Pair {pid}: more than one plasma date in the VOC file (pre {list(r.pre)}, post {list(r.post)})")
                elif len(r.pre) == 1 and len(r.post) == 1 and not r.pre[0] < r.post[0]:
                    errors.append(f"Pair {pid}: plasma pre date {pd.Timestamp(r.pre[0]).date()} is not before "
                                  f"post date {pd.Timestamp(r.post[0]).date()}")
            pre_d = {k: pd.Timestamp(v[0]) for k, v in dates.pre.items() if len(v) == 1}
            post_d = {k: pd.Timestamp(v[0]) for k, v in dates.post.items() if len(v) == 1}
            si['plasma_date'] = [(pre_d if g == 'pre' else post_d).get(pr) if pr else None
                                 for pr, g in zip(si.Pair, si.timepoint_group)]
            si['plasma_date'] = si['plasma_date'].map(lambda d: d.date().isoformat() if d is not None and pd.notna(d) else None)
            paired_sheet = set(si[si.Pair.notna()].Pair)
            nodate = sorted(paired_sheet - set(pre_d) - set(post_d), key=lambda x: int(x[1:]))
            if nodate:
                warnings.append(f"Pairs with no plasma dates in the VOC file (order unchecked): {', '.join(nodate)}")
            notes.append(f"Plasma dates found for {len(paired_sheet & (set(pre_d) | set(post_d)))} of {len(paired_sheet)} pairs")

    if idm:
        idp = rpath(idm['file'] if isinstance(idm, dict) else idm)
        id_col = (idm.get('id_column') if isinstance(idm, dict) else None) or 'record_id'
        sp_col = (idm.get('specimen_column') if isinstance(idm, dict) else None) or 'VOC_specimenid'
        if not idp.exists():
            errors.append(f"ID map not found: {idp}")
        else:
            m = pd.read_csv(idp, dtype=str)
            if id_col not in m.columns or sp_col not in m.columns:
                sys.exit(f"{idp}: expected columns {id_col} and {sp_col}, found {list(m.columns)}")
            m = m.dropna(subset=[id_col, sp_col])
            m['spec'] = m[sp_col].map(norm)
            m['rid'] = m[id_col].str.strip()
            multi = m.groupby('spec').rid.nunique()
            for sp in multi[multi > 1].index:
                errors.append(f"ID map: specimen {sp} is assigned to more than one record_id: "
                              f"{sorted(m[m.spec == sp].rid.unique())}")
            lut = m.drop_duplicates('spec').set_index('spec').rid

            # Person behind each pair, through the pair's urine specimens
            pair_rid = {}
            if voc is not None:
                voc['rid'] = voc.VOC_specimenid.map(lambda x: lut.get(norm(x)))
                for pid, g in voc.groupby('Pair'):
                    rids = sorted(set(g.rid.dropna()))
                    if len(rids) > 1:
                        errors.append(f"Pair {pid}: its urine specimens belong to different people: "
                                      f"{', '.join(f'{a}={b}' for a, b in zip(g.VOC_specimenid, g.rid))}")
                    elif len(rids) == 1:
                        pair_rid[pid] = rids[0]
                unmapped = voc[voc.rid.isna()].VOC_specimenid.tolist()
                if unmapped:
                    warnings.append(f"{len(unmapped)} urine specimens in the VOC file are not in the ID map: {', '.join(map(str, unmapped))}")
                by_person = {}
                for pid, rid in pair_rid.items():
                    by_person.setdefault(rid, []).append(pid)
                for rid, pl in by_person.items():
                    if len(pl) > 1:
                        errors.append(f"record_id {rid} is the person behind more than one pair: {sorted(pl)}")

            # Per sample: the pair's person; otherwise a direct match of the plasma specimen ID
            direct = si.specimen_id.map(lambda x: lut.get(norm(x))) if c_spec else pd.Series(None, index=si.index)
            si['record_id'] = [pair_rid.get(pr) if pr in pair_rid else d for pr, d in zip(si.Pair, direct)]
            for r in si[si.Pair.isin(pair_rid.keys())].itertuples():
                d = direct.loc[r.Index]
                if isinstance(d, str) and d != pair_rid[r.Pair]:
                    errors.append(f"{r.Name}: plasma specimen {r.specimen_id} maps to record {d}, "
                                  f"but pair {r.Pair} belongs to record {pair_rid[r.Pair]}")
            again = si[(si.timepoint_group == 'baseline_unpaired') & si.record_id.isin(set(pair_rid.values()))]
            if len(again):
                warnings.append("Unpaired baseline samples from people who are also in a pair (not independent "
                                "in the pre-exposure comparison with baselines): "
                                + ', '.join(f'{r.Name} (record {r.record_id})' for r in again.itertuples()))
            paired = si[si.Pair.notna()]
            notes.append(f"Person identified for {paired.record_id.notna().sum()} of {len(paired)} paired samples and "
                         f"{si[si.Pair.isna()].record_id.notna().sum()} of {int(si.Pair.isna().sum())} unpaired samples")

    # ------------------------------------------------------------------ QC from MultiQC
    qc_rows, lane_issues, report_paths = [], {}, {}
    for epi, paths in cfg['inputs']['epigenerator_runs'].items():
        mq = rpath(paths['multiqc_data'])
        dd = read_mqc(mq, 'multiqc_bismark_dedup.txt', ['Sample', 'dedup_reads', 'dup_reads_percent'])
        me = read_mqc(mq, 'multiqc_bismark_methextract.txt', ['Sample', 'total_c'])
        # Older MultiQC versions give counts only; Bismark's percentage is meth / (meth + unmeth) * 100
        for ctx in ('cpg', 'chh'):
            pc = f'percent_{ctx}_meth'
            if pc not in me.columns:
                mc, uc = f'meth_{ctx}', f'unmeth_{ctx}'
                if mc not in me.columns or uc not in me.columns:
                    sys.exit(f"{mq}/multiqc_bismark_methextract.txt: need {pc}, or {mc} and {uc}. Columns: {list(me.columns)}")
                me[pc] = (100 * me[mc] / (me[mc] + me[uc])).round(2)
        al = read_mqc(mq, 'multiqc_bismark_alignment.txt', ['Sample', 'total_reads', 'aligned_reads'])
        fq = read_mqc(mq, 'multiqc_fastqc.txt', ['Sample', 'Total Sequences'])
        for d in (dd, me, al, fq):
            d['Name'] = d.Sample.astype(str).str.extract(SAMPLE_PREFIX, expand=False)
        s = dd.set_index('Name')[['dedup_reads', 'dup_reads_percent']].join(
            me.set_index('Name')[['percent_cpg_meth', 'percent_chh_meth', 'total_c']])
        s['epi_run'] = epi
        qc_rows.append(s)

        # Lane completeness: Bismark pairs analysed vs FastQC pairs (R1) per lane
        r1 = fq[fq.Sample.str.endswith('_val_1')].copy()
        r1['key'] = r1.Sample.str.replace(r'_val_1$', '', regex=True)
        al['key'] = al.Sample.str.replace(r'_bismark.*$', '', regex=True).str.replace(r'_val_1$', '', regex=True)
        lanes = r1.merge(al[['key', 'total_reads', 'aligned_reads']], on='key', how='left')
        for _, r in lanes.iterrows():
            lane = (LANE_RE.search(r.Sample) or [None, r.key])[1]
            prob = None
            if pd.isna(r.total_reads):
                prob = f'{lane}: no Bismark report'
            elif r.total_reads / r['Total Sequences'] < qc['lane_min_analysed_frac']:
                prob = f"{lane}: Bismark analysed {r.total_reads / r['Total Sequences']:.0%} of pairs"
            elif 100 * r.aligned_reads / r.total_reads < qc['lane_min_unique_pct']:
                prob = f"{lane}: {100 * r.aligned_reads / r.total_reads:.0f}% unique alignment"
            if prob:
                lane_issues.setdefault(r.Name, []).append(prob)

        cr = rpath(paths['cytosine_reports'])
        if not cr.is_dir():
            errors.append(f"{epi}: cytosine report folder not found: {cr}")
        for n in s.index:
            hits = find_report(str(cr), n)
            if len(hits) == 1:
                report_paths[n] = os.path.relpath(hits[0], REPO_ROOT)
            elif len(hits) > 1:
                errors.append(f"{epi}: {n} matches several cytosine reports: {[os.path.basename(h) for h in hits]}")
            elif cr.is_dir():
                errors.append(f"{epi}: no CpG cytosine report for {n} in {cr}")

    q = pd.concat(qc_rows)
    if q.index.duplicated().any():
        errors.append(f"Samples present in more than one epigenerator run: {sorted(set(q.index[q.index.duplicated()]))}")
    q = q[~q.index.duplicated()]
    q['dedup_pairs_M'] = (q.dedup_reads / 1e6).round(2)
    q['total_c_M'] = (q.total_c / 1e6).round(1)

    # Metadata vs sequencing
    seq_only = sorted(set(q.index) - set(si.Name))
    if seq_only:
        errors.append(f"Sequenced samples missing from the sample sheet: {seq_only}")
    meta_only = sorted(set(si.Name) - set(q.index))
    if meta_only:
        warnings.append(f"In the sample sheet but not sequenced (dropped): {meta_only}")

    t = si.set_index('Name').join(q[['epi_run', 'dedup_pairs_M', 'dup_reads_percent', 'percent_cpg_meth',
                                       'percent_chh_meth', 'total_c_M']], how='inner')
    t['cytosine_report'] = t.index.map(report_paths)
    t['lane_problems'] = t.index.map(lambda n: '; '.join(lane_issues.get(n, [])))

    # ------------------------------------------------------------------ flags
    def flags(r):
        f = []
        chh = r.percent_chh_meth
        if chh >= qc['chh_fail_pct']:
            f.append('failed_conversion')
        elif chh >= qc['chh_poor_pct']:
            f.append('poor_conversion')
        elif chh >= qc['chh_marginal_pct']:
            f.append('marginal_conversion')
        if r.dedup_pairs_M < qc['min_dedup_pairs_M']:
            f.append('low_depth')
        if r.lane_problems:
            f.append('lane_problem')
        if r.name in (qc.get('manual_exclude') or {}):
            f.append('manual_exclude')
        return ' '.join(f)

    t['flags'] = t.apply(flags, axis=1)
    allowed = set(cfg.get('allow_flags') or [])
    t['keep'] = t['flags'].map(lambda s: all(f in allowed for f in s.split()))
    for n, why in (qc.get('manual_exclude') or {}).items():
        if n in t.index:
            notes.append(f"Manual exclusion {n}: {why}")

    # Complete pairs among kept samples
    kept = t[t.keep & t.Pair.notna() & t.timepoint_group.isin(['pre', 'post'])]
    complete = kept.groupby('Pair').timepoint_group.nunique()
    complete = set(complete[complete == 2].index)
    t['pair_complete_kept'] = t.Pair.isin(complete) & t.keep

    # ------------------------------------------------------------------ write
    out = t.reset_index().rename(columns={'index': 'Name'})
    cols = ['Name', 'Career', 'Timepoint_label', 'timepoint_group', 'Pair', 'pair_complete_kept', 'keep', 'flags',
            'lane_problems', 'dedup_pairs_M', 'dup_reads_percent', 'percent_cpg_meth', 'percent_chh_meth', 'total_c_M',
            'input_ng', 'size_selected', 'specimen_id', 'record_id', 'plasma_date', 'epi_run', 'cytosine_report']
    out = out[cols].sort_values('Name')
    st = rpath(cfg['sample_table'])
    st.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(st, sep='\t', index=False)

    # Report
    lines = [f"# Sample table validation: {cfg['run_id']}", '',
             f"Generated {datetime.now():%Y-%m-%d %H:%M} from `{os.path.relpath(si_path, REPO_ROOT)}` and "
             f"{', '.join(cfg['inputs']['epigenerator_runs'])}.", '',
             f"Allowed flags: {', '.join(sorted(allowed)) or 'none'}", '',
             '## Checks', '']
    lines += [f"- ERROR: {e}" for e in errors] or ['- No errors.']
    lines += [f"- WARNING: {w}" for w in warnings]
    lines += [f"- NOTE: {n}" for n in notes]
    lines += ['', '## Samples kept', '', '| Career | Timepoint group | Sequenced | Kept |', '|---|---|---|---|']
    for (c, g), d in out.groupby(['Career', 'timepoint_group']):
        lines.append(f"| {c} | {g} | {len(d)} | {int(d.keep.sum())} |")
    lines += ['', f"Complete pre/post pairs kept: "
              f"{', '.join(f'{c} {len(set(out[(out.Career == c) & out.pair_complete_kept].Pair))}' for c in ['Recruit', 'Instructor'])}",
              '', '## Excluded samples', '', '| Name | Career | Timepoint | Flags | Lane problems |', '|---|---|---|---|---|']
    for _, r in out[~out.keep].iterrows():
        lines.append(f"| {r.Name} | {r.Career} | {r.Timepoint_label} | {r['flags']} | {r.lane_problems} |")
    vr = rpath(cfg['validation_report'])
    vr.write_text('\n'.join(lines) + '\n')

    print('\n'.join(lines))
    print(f"\nWrote {os.path.relpath(st, REPO_ROOT)}\nWrote {os.path.relpath(vr, REPO_ROOT)}")
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()