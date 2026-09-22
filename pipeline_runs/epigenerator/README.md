# Epigenerator: A Practical Workflow for Whole Genome Bisulfite Sequencing (WGBS) Analysis

## Table of Contents

* [Updates](#updates)
* [Project Set-Up](#project-set-up)
	* [Setting Up Your Project Directory](#setting-up-your-project-directory)
	* [LaSalle Lab Layout: Shared Scripts and README](#lasalle-lab-layout-shared-scripts-and-readme)
	* [Installation](#installation)
	* [Genome Preparation](#genome-preparation)
* [1. Preparing Raw Sequences](#1-preparing-raw-sequences)
	* [Unique File Names Are Required](#unique-file-names-are-required)
	* [Option A: Renaming Scripts (Recommended)](#option-a-renaming-scripts-recommended)
		* [Step 0: Soft-Link the Raw Data](#step-0-soft-link-the-raw-data)
		* [Step 1: Build the Renaming Table](#step-1-build-the-renaming-table)
		* [Step 2: Rename the Files in Place](#step-2-rename-the-files-in-place)
		* [Step 3: Place Files and Write `task_samples.yaml`](#step-3-place-files-and-write-task_samplesyaml)
		* [Check the Result](#check-the-result)
	* [Option B: `FASTQ_Me2` (Legacy / SLIMS Downloads)](#option-b-fastq_me2-legacy--slims-downloads)
	* [`task_samples.yaml` Format](#task_samplesyaml-format)
* [2. `CpG_Me2`](#2-cpg_me2)
	* [Running `CpG_Me2` Locally](#running-cpg_me2-Locally)
	* [Running `CpG_Me2` on SLURM (Recommended)](#Running-CpG_Me2-on-SLURM-Recommended)
* [3. `DMRichR`](#3-dmrichr)
* [4. `Comethyl`](#4-comethyl)
* [Interpretting Outputs](#Interpretting-Outputs)
* [Acknowledgements](#acknowledgements)
* [Citation](#citation)

## Updates

The original [CpG_Me](https://github.com/ben-laufer/CpG_Me), written by Dr. Ben Laufer, was intended to merge WGBS sequencing lanes and generate cytosine reports to be used in further analysis, such as in DMRichR and Comethyl. Upon trying the pipeline for the first time, I ran into some complications specific to my data (which was not downloaded from SLIMS), and I had the idea to make it more modular and a little easier to use for first-time users. Here are some of the updates I made:

* `FASTQ_Me2` has been rewritten as a python script
* Locally downloaded sequencing data can now be used in `FASTQ_Me2`
* Users will only interact with command line prompts instead of adapting any of the scripts
* `CpG_Me2` can be run locally (i.e. on screen or at your terminal) as a script instead of only through SLURM 
* `CpG_Me2` has been written as a snakemake file, yielding the following advantages:
    * Multi-threading of jobs can be handled locally or on SLURM
    * Jobs remove corrupted intermediate files if they fail
    * Snakemake's "memory" prevents re-running of samples and files that have already been generated
* Three renaming scripts in the shared `scripts/` folder (`pipeline_runs/epigenerator/scripts/` for LaSalle Lab) are now the recommended way to prepare raw sequences, especially for data spread across multiple lanes, runs, or flow cells (see [Option A](#option-a-renaming-scripts-recommended)). They guarantee that every FASTQ has a unique name before `CpG_Me2` is run. `FASTQ_Me2` is kept for downloading data from SLIMS.

## Project Set-Up

### Setting Up Your Project Directory

Clone the repository using your project name in the directory you plan to host the project:

```
git clone https://github.com/vhaghani26/epigenerator {project_name}
```

Enter the directory

```
cd {project_name}
```

### LaSalle Lab Layout: Shared Scripts and README

In LaSalle Lab projects, each pipeline run is its own epigenerator clone inside `pipeline_runs/epigenerator/`. The renaming scripts and this README live **once**, one level above the runs, and are shared by every run:

```
{lab_project}/pipeline_runs/epigenerator/
├── README.md                  # this file, shared by all runs
├── scripts/                   # renaming scripts, shared by all runs
│   ├── 01_create_transfers_list_csv.py
│   ├── 02_data_renaming.py
│   └── 03_removeR_change_to_fq_gen_yaml.py
├── EPI_AZ_01/                 # one epigenerator clone per run ("run directory")
│   ├── 01_raw_sequences/
│   │   ├── 2026_LACOFD_WGBS_cellfree_Logan -> /quobyte/lasallegrp/data/...   # soft link (Step 0)
│   │   └── *_1.fq.gz, *_2.fq.gz                                             # hard links (Step 3)
│   └── task_samples.yaml
├── EPI_AZ_PILOT_01/
├── EPI_LA_01/
└── EPI_LA_02/
```

Clone a new run next to the existing ones (`git clone https://github.com/vhaghani26/epigenerator EPI_XX_01`) and run everything from inside that run directory. The shared scripts are therefore called as `../scripts/...`.

Why keep one shared copy instead of a copy in each run:

* **Consistency.** Every run uses the same code, so runs are processed identically and a bug fix applies to all future runs at once. Separate copies drift apart over time, and it becomes unclear which run used which version.
* **Lab code stays separate from upstream code.** Each run directory is a clone of the upstream epigenerator repository. Lab-specific scripts placed inside a clone get mixed in with upstream files and are lost if the clone is deleted or re-cloned.
* **Version control in one place.** The shared `scripts/` and `README.md` are tracked once in the lab project's git repository, so every change is recorded with its commit history.
* **Traceability.** Because shared scripts can change later, record which version a run used by running this from the run directory before Step 1 (requires the scripts to be committed):

  ```
  git -C ../scripts log -1 --format='%H %ad' -- . > scripts_version.txt
  ```

### Installation

Run the following command in your project directory. It will clone the conda environment with all dependencies needed in order to run the workflow outlined here. This creates an environment called `epigenerator`. If you would like to change the name, feel free to do so where the command says `epigenerator`. Please note that this may take quite a few minutes to run.

```
conda env create -f 00_software/environment.yml --name epigenerator
```

Activate your environment using

```
conda activate epigenerator
```

Run everything downstream of this point in this conda environment. Note that you must activate this environment every time you restart your terminal.

**For LaSalle Lab**

The environment is already installed and ready in a shared space, so all you will need to do is run:

```
conda activate /quobyte/lasallegrp/programs/.conda/epigenerator
```

### Genome Preparation

**1. `01_genomes/` Setup**

Several steps require genomes for alignment. One of the steps, FastQ-Screen, allows you to align your reads to multiple genomes to check for sources of contamination. In your project directory, create a subdirectory called `01_genomes/`:

```
mkdir 01_genomes/
```

Within `01_genomes/`, create subdirectories corresponding to each genome of interest. The directory structure should look something like the following:

```
project_directory/
	01_genomes/
		hg38/
		Lamba/
		mm10/
		PhiX/
		rheMac10/
		rn6/
```

Each subdirectory containing your genome of interest should contain the appropriate genome files as described by the [FastQ-Screen documentation](https://stevenwingett.github.io/FastQ-Screen/). Please do not worry about downloading or installing FastQ-Screen. The environment you created has all the software you will need. Activate the environment you created in the previous section, download the reference genome, and index the reference genome using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome). If you have already processed and indexed a genome(s), then you can instead alias those files/directories into `01_genomes` as shown in the "For LaSalle Lab" section below.

**2. `fastq_screen.conf` Setup**

In the directory `00_software`, locate and open the file `fastq_screen.conf`. Change lines 34+ to reflect the genomes and location of the genomes that you are aligning to. Provided paths can be absolute or relative paths.

**For LaSalle Lab**

Create `01_genomes/`:

```
mkdir 01_genomes
```

Change into the directory:

```
cd 01_genomes
```

Run the following commands to link our standard genome directories into this directory:

```
ln -s /quobyte/lasallegrp/genomes/hg19 .
ln -s /quobyte/lasallegrp/genomes/hg38 .
ln -s /quobyte/lasallegrp/genomes/Lambda .
ln -s /quobyte/lasallegrp/genomes/mm10 .
ln -s /quobyte/lasallegrp/genomes/PhiX .
ln -s /quobyte/lasallegrp/genomes/rheMac10 .
ln -s /quobyte/lasallegrp/genomes/rn6 .
```

You do not need to update `fastq_screen.conf`. 

## 1. Preparing Raw Sequences

`CpG_Me2` aligns each lane of each sample separately and then merges the aligned lanes, which is more efficient and needs smaller resource requests than merging raw FASTQs first. Before running `CpG_Me2` you need two things:

1. All raw FASTQs in one directory, `01_raw_sequences/`, named `{sample}_1.fq.gz` (forward) and `{sample}_2.fq.gz` (reverse), where `{sample}` identifies **one lane of one sample**.
2. A configuration file, `task_samples.yaml`, in the project directory that lists every `{sample}`. `CpG_Me2` will not run without it.

### Unique File Names Are Required

> [!IMPORTANT]
> **Every FASTQ must have a unique file name across all lanes, runs, and flow cells before it is placed in `01_raw_sequences/`.**

Illumina file names (`{sample}_S{n}_L{lane}_R{1|2}_001.fastq.gz`) only encode the sample ID, sample index, and lane number. They do **not** encode the run or flow cell. When the same samples are sequenced on more than one run, different runs routinely produce identical file names. For example, each of these folders can contain its own `LA001_S1_L004_R1_001.fastq.gz`:

```
01_raw_sequences/2026_LACOFD_WGBS_cellfree_Logan/     (soft link to the lab data folder, see Step 0)
├── 260714_DTSA1275_1284_NovaX25B/
│   └── JLLW_Nova1479P_Williams/            LA001_S1_L004_R1_001.fastq.gz ...
├── 260819_DTSA1302_1303_1304_1305_NovaX25B/
│   └── JLLW_Nova1479P_Williams_L2/         LA001_S1_L004_R1_001.fastq.gz ...
└── 260909_DTSA1315_1316_1317_1318_NovaX25B/
    ├── JLLW_Nova1479P_Williams_L3/         LA001_S1_L004_R1_001.fastq.gz ...
    ├── JLLW_Nova1479P_Williams_L4/
    ├── JLLW_Nova1479P_Williams_L5/
    └── JLLW_Nova1479P_Williams_L6/
```

If these files are collected into `01_raw_sequences/` without renaming, same-named files overwrite each other or are skipped. The lanes also collapse into a single entry in `task_samples.yaml`. Either way, sequencing data is **silently lost** from the analysis. The renaming scripts below prevent this by adding the name of the lane/run folder to every file name.

### Option A: Renaming Scripts (Recommended)

Use this for any locally stored data and **always** when data comes from more than one lane, run, or flow cell. The three scripts are in the shared `scripts/` folder one level above the run directories (see [LaSalle Lab Layout](#lasalle-lab-layout-shared-scripts-and-readme)). Run them in order **from your run directory** (e.g. `pipeline_runs/epigenerator/EPI_AZ_01/`) with the environment activated. Replace `{batch}` with the folder that holds your sequencing runs (e.g. `2026_LACOFD_WGBS_cellfree_Logan`).

#### Step 0: Soft-Link the Raw Data

For LaSalle Lab, raw sequencing data is delivered to and kept in the lab data folder (`/quobyte/lasallegrp/data/`). Do **not** copy or move it into the run. Instead, soft-link (symlink) the batch folder into `01_raw_sequences/`:

```
mkdir -p 01_raw_sequences
ln -s /quobyte/lasallegrp/data/{batch} 01_raw_sequences/{batch}
ls -l 01_raw_sequences/          # should show: {batch} -> /quobyte/lasallegrp/data/{batch}
```

Use the **absolute** path to the data, as shown. A relative link breaks if the run directory is moved.

Why soft-link instead of copying or moving:

* **One source of truth.** Raw data stays in the lab's data folder, where lab storage, permissions, and backups are managed, and where other lab members expect to find it. Moving it into a run directory scatters raw data across analysis folders.
* **No duplicated storage.** WGBS FASTQs take hundreds of gigabytes to terabytes per batch. A copy doubles that for no benefit, and a soft link takes essentially no space.
* **Several runs can share one dataset.** Pilot runs, re-runs, and runs with different settings (e.g. `EPI_AZ_PILOT_01`, `EPI_AZ_01`) can all link the same data without extra copies.
* **The raw data is protected from run cleanup.** Deleting a run directory, or the intermediate files `CpG_Me2` produces, removes only the link, not the data.
* **Provenance is visible.** `ls -l 01_raw_sequences/` shows exactly which data folder a run used.

> [!WARNING]
> **A soft link points to the real files, so Step 2 renames the files in the lab data folder itself.** This is a permanent change to the shared raw data, which is why the renaming CSV must be kept (see [Step 2](#step-2-rename-the-files-in-place)). Confirm that renaming is acceptable for that data folder before running Step 2. It also means that if the data were already renamed for an earlier run, Steps 1–2 are already done: skip straight to Step 3.

> [!CAUTION]
> To remove the link, use `unlink 01_raw_sequences/{batch}` (or `rm` with no trailing slash). Never run `rm -r` on the link, especially with a trailing slash (`rm -r 01_raw_sequences/{batch}/`): that can follow the link and delete the raw data.

You can link the whole batch folder (as above) or link run folders individually inside `01_raw_sequences/{batch}/`; the scripts follow soft links either way. Each lane/run folder must have a **distinct name**, because that name becomes part of every file name. Folder names that are already distinct (such as `JLLW_Nova1479P_Williams` and `JLLW_Nova1479P_Williams_L2`) can be left as delivered.

#### Step 1: Build the Renaming Table

Script: `01_create_transfers_list_csv.py`. This searches `{batch}` recursively and writes a CSV mapping each original file name to a new name with its folder name inserted before the read direction. It does not rename anything.

```
python3 ../scripts/01_create_transfers_list_csv.py \
  --indir 01_raw_sequences/{batch} \
  --output {batch}_transfers.csv
```

Example: `LA001_S1_L004_R1_001.fastq.gz` in `JLLW_Nova1479P_Williams_L2/` becomes `LA001_S1_L004_JLLW_Nova1479P_Williams_L2_R1_001.fastq.gz`.

Open the CSV and check the `New_File_Name` column. The script prints a warning if two files would get the same new name (i.e. two folders share a name). It also lists any files that do not follow the Illumina naming pattern, which includes files that were already renamed. Keep this CSV: it is the record of original → new names.

#### Step 2: Rename the Files in Place

Script: `02_data_renaming.py`. Preview first, then apply:

```
python3 ../scripts/02_data_renaming.py --file {batch}_transfers.csv --dry-run
python3 ../scripts/02_data_renaming.py --file {batch}_transfers.csv
```

The files are renamed inside their original lane/run folders; with a soft-linked batch (Step 0), that means in the lab data folder. Existing files are never overwritten unless you add `--overwrite`. Because the renaming changes the shared raw data, also keep a copy of the CSV next to the data so anyone using the data folder can trace the original names:

```
cp {batch}_transfers.csv /quobyte/lasallegrp/data/{batch}/
```

#### Step 3: Place Files and Write `task_samples.yaml`

Script: `03_removeR_change_to_fq_gen_yaml.py`. This converts `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` to `*_1.fq.gz` / `*_2.fq.gz` and **hard-links** each file into `01_raw_sequences/` under its final name. It then writes `task_samples.yaml` to the project directory. Hard links point directly at the files in the lab data folder, even when the data were reached through a soft link. They use no extra disk space, keep working if the soft link is removed, and deleting the hard links in `01_raw_sequences/` does not delete the raw data. The script only looks at files **directly** inside the folders you pass it, so pass the lane/run folders themselves. With the layout above, the shell glob `*/*/` does this.

Dry run (default; prints the plan and a preview of the YAML, changes nothing):

```
python3 ../scripts/03_removeR_change_to_fq_gen_yaml.py \
  01_raw_sequences/{batch}/*/*/ \
  --target 01_raw_sequences
```

Apply:

```
python3 ../scripts/03_removeR_change_to_fq_gen_yaml.py \
  01_raw_sequences/{batch}/*/*/ \
  --target 01_raw_sequences --apply
```

Notes on step 3:

* If any two files would end up with the same final name, the script stops before changing anything and lists the conflicting files. Go back to Steps 1–2.
* `Undetermined_*` reads are skipped by default (`--include-undetermined` to keep them).
* The sample list is built from **every** `*.fq.gz` in `01_raw_sequences/`, including files from earlier batches; the script reports any that are not from the current batch.
* The genome defaults to `hg38` (`--genome mm10`, etc. to change it).
* Hard links require the data folder and the run directory to be on the same filesystem (both under `/quobyte/lasallegrp/`). If linking fails, use `--mode symlink`. Do **not** use `--mode move` with soft-linked data: it would move the raw data out of the lab data folder.
* Re-running is safe: files already linked are detected and skipped, and the YAML is only written if every file was placed successfully.

#### Check the Result

```
ls 01_raw_sequences/*.fq.gz | wc -l          # should equal 2 x the number of samples in the YAML
stat -c '%h %n' 01_raw_sequences/*.fq.gz | head   # link count of 2 = hard link to the original
cat task_samples.yaml
```

Each YAML entry is one lane of one sample (e.g. `LA001_S1_L004_JLLW_Nova1479P_Williams_L2`); `CpG_Me2` merges lanes after alignment.

### Option B: `FASTQ_Me2` (Legacy / SLIMS Downloads)

`FASTQ_Me2` was previously the standard first step. It is still useful for **downloading data from SLIMS**, but for local data from multiple lanes or runs, use [Option A](#option-a-renaming-scripts-recommended). If you do use `FASTQ_Me2` on local data, first confirm that every file name is unique (see [Unique File Names Are Required](#unique-file-names-are-required)).

**SLIMS Data**

If your data is on SLIMS, you will be prompted for your SLIMS string and SLIMS directory. This triggers the data download in your project directory. It also computes the MD5 checksum and compares it to the expected checksum hashes to ensure that the files were not corrupted during download.

**Local Data**

If you've already downloaded your data, you will be prompted to give the absolute path to the raw data.

`FASTQ_Me2` does not merge lanes; it generates `task_samples.yaml`, asking for your confirmation along the way to make sure that your samples are being handled correctly. If you have local data instead of data from SLIMS, put all your raw sequence data in the same directory and name it `01_raw_sequences/`. Note that you can use symbolic links within this directory to link to the original files as well.

**Running `FASTQ_Me2`**

To run `FASTQ_Me2`, run the following in your project directory:

```
python3 01_FASTQ_Me2.py
```

### `task_samples.yaml` Format

If you would like to skip both options above, you can write the configuration file yourself. It must be named `task_samples.yaml`, be placed in the project directory, and be formatted as follows:

```
---
sequence_data: /quobyte/lasallegrp/Viki/project_name/01_raw_sequences
genome: mm10
samples:
  - FA114_FKDN220207669-1A_H5GY3DSX3_L3
  - FA114_FKDN220207669-1A_H5GY3DSX3_L2
  - FA114_FKDN220207669-1A_H5JF2DSX3_L1
```

The samples list can go on for as many as you need, but they should include the full file name EXCEPT the read orientation and file extension (i.e. `_1.fq.gz`, `_2.fq.gz`, `_1.fastq.gz`, `_2.fastq.gz`, `_R1.fq.gz`, `_R2.fq.gz`, `_R1.fastq.gz`, `_R2.fastq.gz`, `_R1_001.fq.gz`, `_R2_001.fq.gz`, `_R1_001.fastq.gz`, `_R2_001.fastq.gz`). The same uniqueness rule applies: every sample entry must correspond to exactly one pair of files.

## 2. `CpG_Me2`

`CpG_Me2` carries out a number of steps to process your sample, which can be seen in the figure below.

![Workflow](https://github.com/ben-laufer/CpG_Me/blob/master/Examples/CpG_Me_Flowchart.png)

Note that **02_CpG_Me2_PE is for use with paired-end data and 02_CpG_Me2_SE is for use with single-end data**. Please substitute the script names as needed in the below commands. I have gone with the deafult of paired-end data in this tutorial since paired-end data is more commonly used than single-end data in our lab.

There are two ways that you can run `CpG_Me2`. 

### Running `CpG_Me2` Locally

Since all the setup was done in [Section 1](#1-preparing-raw-sequences) (`task_samples.yaml` must be in your project directory), all you have to do now is run the following to run `CpG_Me2` from start to finish:

```
snakemake -j 1 -p -s 02_CpG_Me2_PE
```

The `-j` option for `jobs` means how many jobs are able to be run in parallel. The most resource-intensive step uses ~50-100 GB of RAM per sample, so be mindful of the resources you have available if you choose to run more than one job at a time. If, for some reason, you are disconnected from your terminal or your job fails, re-run the above command and it will pick up where it left off. It even deletes possibly corrupted files from where it was cut off to ensure ALL outputs are properly generated. 

Notice some issues here.

1. Because the alignment step uses so many resources, you are limited by the number of jobs you can run at the same time due to required resource allocations
2. If you get disconnected from your terminal, your job fails. This means that you need to be logged in and have an active terminal for days at a time. This can be circumvented by running `CpG_Me2` in `screen`, but you will still need the proper resource allocation for that length of time.

For the above reasons, it is HIGHLY recommended that you instead run it via SLURM. This is similar to what was initially written by Ben, where individual jobs get submitted to SLURM.

### Running `CpG_Me2` on SLURM (Recommended)

In the directory `00_slurm/`, there is a file named `config.yaml`. You will need to modify two things:

1. Update your SLURM partition for the **two** lines (line 18 and line 19) containing `--account` and `--partition` by inputting a string. This will look something like `--partition=production` 
2. Change the `conda_prefix` (line 33). It should look something like `/software/anaconda3/4.8.3/lssc0-linux/`, `/home/vhaghani/anaconda3/`, or `/quobyte/lasallegrp/programs/.conda/`

Once you have updated `config.yaml`, go back to your project directory. Snakemake manages the submission of jobs, so wherever you run it, it will need to stay open. As such, I recommend running it in [screen](https://linuxize.com/post/how-to-use-linux-screen/). Activate the conda environment (confirm you are in the environment if you are using screen).  When you are ready, run:

```
snakemake -s 02_CpG_Me2_PE --profile 00_slurm/
```

**For LaSalle Lab**

There are some bugs with the way your home directory is mounted on the Genome Center cluster. This will result in Snakemake erroring out when it tries to cache files. To fix it, we will need to make some changes. Briefly, you are creating a cache directory within your project directory. Since this is in a shared lab space, it is properly mounted for use. You are giving full permissions to the cache directory to ensure that Snakemake can store files in the directory. Then, you are setting the Snakemake cache variable to the proper directory so Snakemake knows to use that directory for cache. After these changes are made, you can submit to SLURM. Please run the following in your project directory (ideally in screen with your conda environment activated):

```
cd 00_slurm/
mkdir .cache
chmod -R 777 .cache
cd ..
```

Change the paths to whatever path you are using and put the following in your `.bashrc` or `.profile` (or whatever you're using):

```
export SNAKEMAKE_OUTPUT_CACHE=/quobyte/lasallegrp/{your_directory}/{your_project}/00_slurm/.cache
export XDG_CACHE_HOME=/quobyte/lasallegrp/{your_directory}/{your_project}/00_slurm/.cache
```

Note that these will have to be updated every time you run `CpG_Me2`, but it is currently the only workaround I'm aware of for dealing with the issues the Genome Center cluster has.

Then, run the following:

```
snakemake -s 02_CpG_Me2_PE --profile 00_slurm/ --cache 00_slurm/.cache 
```

## 3. `DMRichR`

I am hoping to one day better incorporate `DMRichR` here, but in the meantime. please refer to Dr. Ben Laufer's [`DMRichR` package](https://github.com/ben-laufer/DMRichR) for further analysis. This package will aid in the identification of differentially methylated regions within your dataset. 

If you are interested in running `DMRichR` on SLURM, there is a template SLURM script in `00_slurm/` titled `DMRichR_SLURM_template.sh`. If you are using it, please ensure that you change the appropriate variables, including the `#SBATCH` parameter section and the path to the `DM.R` script.

## 4. `Comethyl`

If you are interested in weighted region comethylation network analysis, please refer to Dr. Charles Mordaunt's [`Comethyl` package](https://github.com/cemordaunt/comethyl) for further analysis.

## Interpretting Outputs

The most important outputs you will find are `01_raw_sequences/` and `08_cytosine_reports/`. You can feel free to delete the intermediate files, but they are included for your reference in case something goes wrong when you try to run `CpG_Me2`. The `03_screened` and `09_multiqc` directories may also be good to keep so you can check on the qualities of your samples at various stages. Below is a description of each of the output directories and what is contained in each. Note that the numbers preceding the directory names are reflective of the order they were generated.

### `00_std_err_logs/`

This directory contains the standard error per job per sample. If, for some reason, you encounter an error for a sample, this allows you to view the full error message and debug.

### `00_time_logs/`

This directory documents the time it took each job to run. This could help inform future decisions about time allocations in SLURM should you need to adjust the time for a submitted job.

### `01_genomes/`

This directory contains the primary genome you want to align your data to as well as whatever genomes you want to screen your samples for using FastQ-Screen.

### `01_raw_sequences/`

This directory contains the raw sequences for each lane of each sample, with unique names, separated into forward (`{sample}_1.fq.gz`) and reverse (`{sample}_2.fq.gz`) reads. When using the renaming scripts, these are hard links to the raw files in the lab data folder, and `01_raw_sequences/{batch}` is the soft link to that folder.

### `02_trimmed/`

This contains the reads after they have been trimmed, meaning that the adapters and low quality base pairs at sequence ends have been removed.

### `03_screened/`

All samples are screened against some different background genomes to determine sample read origins. This could be helpful if you are concerned about possible contamination of your samples.

### `04_aligned/`

This directory contains the aligned sequences, meaning that samples are mapped against your reference genome.

### `05_deduplicated/`

This contains your reads with PCR duplicates removed.

### `06_sorted/`

This contains your aligned and deduplicated sequence reads sorted in chromosome coordinate order.

### `06_nt_coverage/`

Nucleotide coverage of reads are calculated for input BAM files and output here.

### `06_methylation/`

This contains the various outputs of Bismark associated with extracting methylation information from each sample.

### `07_size_metrics/`

According to the documentation, this "provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries." 

### `08_cytosine_reports/`

This directory contains the cytosine reports that are required for downstream analysis.

### `09_multiqc/`

This contains quality read outs for the eligible files maintained in the other directories. It can be helpful to view if any of your samples are being problematic.

### `logs/`

This directory is generated when you use SLURM to run `CpG_Me2`. It contains the `err` and `out` files typically generated by SLURM per job submission, which corresponds to individual rules run per sample. 

## Acknowledgements

This work was largely adapted from Dr. Ben Laufer's original [CpG_Me Program](https://github.com/ben-laufer/CpG_Me). These updates could not have been made without the help of Jules Mouat and Aron Mendiola.

## Citation

If you are using `epigenerator`, please cite the following (noting that it may be subject to change upon further publication): 

* Laufer BI*, Neier KE*, Valenzuela AE, Yasui DH, Lein PJ, LaSalle JM. Placenta and Fetal Brain Share a Neurodevelopmental Disorder DNA Methylation Profile in a Mouse Model of Prenatal PCB Exposure. Cell Reports, 2022. doi: [10.1016/j.celrep.2022.110442](https://www.sciencedirect.com/science/article/pii/S2211124722001693?via%3Dihub)

* Krueger F, Andrews SR. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics, 2011. doi: [10.1093/bioinformatics/btr167](https://doi.org/10.1093/bioinformatics/btr167)

* Matrin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 2011. doi: [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)

* Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016. doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

* Haghani V, Mouat JS. 2023. Epigenerator. GitHub. https://github.com/vhaghani26/epigenerator.

If you are using `Comethyl`, also make sure to cite the following:

* Mordaunt CE, Mouat JS, Schmidt RJ, and LaSalle JM. (2022) Comethyl: a network-based methylome approach to investigate the multivariate nature of health and disease. Briefings in Bioinformatics bbab554.