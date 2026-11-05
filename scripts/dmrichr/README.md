# DMRichR runs

Differentially methylated region (DMR) analysis with [DMRichR](https://github.com/ben-laufer/DMRichR), set up so that:

* sample selection and pre/post labels are checked before anything runs;
* every comparison has its own folder, named by one rule;
* cytosine reports are **hard-linked**, never copied, so no comparison folder uses extra disk space and the epigenerator outputs are never moved or renamed;
* every run records which samples, settings and script version it used.

Code and configs are version-controlled in `scripts/dmrichr/`; each run is a folder under `results/dmrichr/` that git ignores.

## Layout

Code and run configs are tracked in `scripts/dmrichr/`. Runs are written to `results/dmrichr/`, which git ignores. DMRichR itself comes from the existing pixi environment, so nothing is cloned or copied per run.

```text
scripts/dmrichr/                   # tracked
├── README.md
├── config/                        # one YAML per run
│   └── DMR_AZ_01.yaml
├── 01_build_sample_table.py       # join sample sheet + QC, validate labels, decide who is kept
├── 02_setup_comparisons.py        # one folder per comparison: hard links, sample_info.xlsx, sbatch
└── 03_check_pair_sex.py           # genetic check that each pre/post pair is the same sex

results/dmrichr/                   # not tracked
└── DMR_AZ_01/                     # one run directory per config
    ├── config_snapshot.yaml
    ├── sample_table_snapshot.tsv
    ├── scripts_version.txt
    ├── pixi.lock                  # copy of the environment lock file used
    ├── comparisons_summary.tsv
    ├── recruit_post_vs_pre/
    │   ├── LA003_merged_name_sorted.CpG_report.txt.gz   # hard link
    │   ├── ...
    │   ├── sample_info.xlsx        # DMRichR design matrix
    │   ├── sample_list.tsv         # same samples with QC and labels
    │   ├── comparison.yaml
    │   ├── run_dmrichr.sbatch
    │   ├── dmrichr_version.txt     # DMRichR, R and DM.R path used
    │   ├── logs/
    │   └── DMRs/ Blocks/ Global/ RData/ ...   # written by DMRichR
    ├── instructor_post_vs_pre/
    ├── pre_recruit_vs_instructor/
    ├── pre_recruit_vs_instructor_withbaseline/
    └── post_recruit_vs_instructor/
```

The validated sample table and its report are written to `data/processed/wgbs/sample_tables/`, alongside other analysis-ready data.

## Naming

**Run ID:** `DMR_<GROUP>_<RUN>`, following the repository's other identifiers (`EPI_AZ_01`, `DNAEXT_LA_01`). The config file is named after the run (`scripts/dmrichr/config/DMR_AZ_01.yaml`) and the run directory takes the same name (`results/dmrichr/DMR_AZ_01/`). A new sample selection or new DMRichR settings means a new run ID, so earlier results are never overwritten.

**Comparison folders:** `<subset>_<test group>_vs_<reference group>`, lower case, words joined by `_`. The reference group is the one DMRichR treats as control, so results read "test vs reference".

| Folder                                     | Samples                                      | Test vs reference     | Design            |
| ------------------------------------------ | -------------------------------------------- | --------------------- | ----------------- |
| `recruit_post_vs_pre`                    | Recruit complete pairs                       | Post vs Pre           | adjusted for Pair |
| `instructor_post_vs_pre`                 | Paired-instructor complete pairs             | Post vs Pre           | adjusted for Pair |
| `pre_recruit_vs_instructor`              | Pre-exposure samples                         | Recruit vs Instructor | unpaired          |
| `pre_recruit_vs_instructor_withbaseline` | Pre-exposure + unpaired baseline instructors | Recruit vs Instructor | unpaired          |
| `post_recruit_vs_instructor`             | Post-exposure samples                        | Recruit vs Instructor | unpaired          |

Comparisons are defined in the config, not in code. To add one, add a block under `comparisons:` and re-run step 2.

## How the labels are handled

* **Pre / post.** "pre-exposure" and paired "baseline" (a baseline sample that has a pair number) are the *pre* sample of a pair; "post-exposure" is *post*. "not paired - baseline" samples are unpaired baseline instructors.
* **Pairs.** `Pair` is written as `P<number>`. For within-person comparisons only complete pairs are kept (both samples pass QC), and `Pair` is passed to DMRichR as an adjustment covariate.
* **Reference level.** DMRichR treats the alphabetically first label as the control. "Post" sorts before "Pre", so the script writes `1_Pre` and `2_Post`. "Instructor" already sorts before "Recruit", so those labels are kept as they are.
* **Sample names.** DMRichR takes the sample name from the cytosine report file name up to the first `_` (`LA003_merged_name_sorted.CpG_report.txt.gz` → `LA003`), so names in the sample sheet must contain only letters and digits.

## Checks run before anything is linked (step 1)

Errors stop the script; warnings are listed in the report.

| Check                                                                                                     | Level   |
| --------------------------------------------------------------------------------------------------------- | ------- |
| Sample names unique, letters and digits only, none a prefix of another                                    | error   |
| Career is Recruit or Instructor; Timepoint is one of the four known labels                                | error   |
| Each pair has one career, at most one pre and at most one post                                            | error   |
| "not paired - baseline" samples have no pair number                                                       | error   |
| Within a specimen series (e.g. two `REP-` IDs), the pre specimen is numbered before the post            | error   |
| Every sequenced sample is in the sample sheet; every sample has exactly one CpG report                    | error   |
| With the VOC/PAH file: each pair's plasma pre date is before its post date                                | error   |
| With the VOC/PAH file + ID map: a pair's urine specimens belong to one person; no person behind two pairs | error   |
| Pairs with no dates; urine specimens missing from the ID map; unpaired baselines who are also in a pair   | warning |
| Paired "baseline" samples used as pre                                                                     | warning |
| Samples in the sheet that were not sequenced (e.g. LA074)                                                 | warning |

Step 3 adds a check that uses the data itself: within each pair, both samples must have the same genetic sex (from chrY vs chrX coverage). A mismatch means a sample swap or a wrong pair number. A match is necessary but not sufficient; a SNP-genotype comparison between pair members is the full identity check.

## QC flags and which samples are kept

Set in the config under `qc:` and `allow_flags:`. Each sample gets zero or more flags:

| Flag                    | Rule (defaults)                                                                   |
| ----------------------- | --------------------------------------------------------------------------------- |
| `failed_conversion`   | CHH methylation ≥ 20%                                                            |
| `poor_conversion`     | CHH methylation 5–20%                                                            |
| `marginal_conversion` | CHH methylation 2–5%                                                             |
| `low_depth`           | fewer than 10M deduplicated read pairs                                            |
| `lane_problem`        | a lane with no Bismark report, < 98% of pairs analysed, or < 20% unique alignment |
| `manual_exclude`      | listed under `qc: manual_exclude` with a reason                                 |

CHH methylation is taken from Bismark's splitting reports (`percent_chh_meth` in `multiqc_bismark_methextract.txt`). A sample is kept only if every flag it carries is in `allow_flags`. `DMR_AZ_01` allows `low_depth` and `marginal_conversion`.

## Environment

DMRichR is used from the lab's existing pixi environment. The generated sbatch script activates it with the `env_setup` line in the config (`pixi shell-hook --manifest-path …`) and finds `DM.R` inside the installed DMRichR package (`system.file("exec", "DM.R", package = "DMRichR")`), so no separate DMRichR clone is needed. Each run copies the environment's `pixi.lock` and writes `dmrichr_version.txt`, which record exactly which package versions produced the results.

The setup scripts themselves need Python 3 with `pandas`, `pyyaml` and `openpyxl`.

## Running

From the repository root:

```bash
# 1) Build and validate the sample table
python3 scripts/dmrichr/01_build_sample_table.py --config scripts/dmrichr/config/DMR_AZ_01.yaml
#    read data/processed/wgbs/sample_tables/DMR_AZ_01_validation_report.md

# 2) (recommended) Pair sex check, as a SLURM job
sbatch --account=lasallegrp --partition=high --cpus-per-task=16 --mem=8G --time=06:00:00 \
    --wrap "python3 scripts/dmrichr/03_check_pair_sex.py \
            --config scripts/dmrichr/config/DMR_AZ_01.yaml --threads 16 --all"

# 3) Preview the comparison folders (dry run)
python3 scripts/dmrichr/02_setup_comparisons.py --config scripts/dmrichr/config/DMR_AZ_01.yaml

# 4) Create them, then submit one or all
python3 scripts/dmrichr/02_setup_comparisons.py --config scripts/dmrichr/config/DMR_AZ_01.yaml --apply
sbatch results/dmrichr/DMR_AZ_01/recruit_post_vs_pre/run_dmrichr.sbatch
#    or: add --submit to step 4 to submit every comparison
```

Commit the scripts and config before step 4 so `scripts_version.txt` records the exact version used.

Hard links need the cytosine reports and the run directory on the same filesystem (both under `/quobyte/lasallegrp/projects/LACOFD`). Deleting a run directory removes only the links; the epigenerator files remain.

## Results

DMRichR writes its output inside each comparison folder under `results/dmrichr/<RUN_ID>/<comparison>/`. Figures or tables chosen for a report or manuscript can be referenced from there directly; there is no second copy to keep in sync.
