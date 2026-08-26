# LACOFD

## Examining longitudinal changes in DNA methylation in firefighters exposed to products of combustion

This repository contains data-processing workflows, quality-control procedures, and analysis code for the **Los Angeles County Fire Department (LACoFD) firefighter DNA methylation study**.

The study is part of the California Firefighter Cancer Prevention and Research Program and examines molecular changes associated with occupational exposure to products of combustion. Firefighters can be repeatedly exposed to potentially harmful chemicals during training and active duty. The study uses longitudinal samples collected before and after exposure to investigate whether these exposures are associated with changes in DNA methylation.

The project includes new firefighter recruits and experienced instructors, who differ in their histories and frequency of occupational exposure. Chemical exposure measurements currently include volatile organic compounds (VOCs) and polycyclic aromatic hydrocarbons (PAHs). Cell-free DNA (cfDNA) is isolated from plasma and characterized before whole-genome bisulfite sequencing (WGBS). DNA methylation profiles generated from WGBS will subsequently be examined in relation to exposure measures and participant characteristics.

The long-term goal of the project is to better understand molecular changes associated with repeated occupational exposures in firefighters and to identify potential blood-based biomarkers relevant to cancer risk.

## Study workflow

The primary workflow is:

```text
LACoFD firefighter recruitment
        │
        ├── New recruits
        └── Experienced instructors
                │
                ▼
       Pre-exposure plasma
                │
                ▼
       Training / occupational exposure
                │
                ▼
       Post-exposure plasma
                │
                ▼
        cfDNA isolation
                │
                ▼
      DNA quantity and quality QC
         ├── Qubit
         └── Bioanalyzer
                │
                ▼
              WGBS
                │
                ▼
     Epigenerator processing
                │
                ▼
       Cytosine reports
                │
                ▼
   DNA methylation analysis
```

Downstream DNA methylation analyses are under development. Approaches being considered include region-based differential methylation analyses and DNA methylation network analyses using tools such as **DMRichR** and **CoMethyL**.

## Repository organization

The repository separates original data, intermediate processing, analysis-ready datasets, computational pipeline runs, scripts, and generated results.

```text
LACOFD/
├── config/
│
├── data/
│   ├── raw/
│   ├── interim/
│   ├── processed/
│   └── metadata/
│
├── docs/
│
├── pipeline_runs/
│   └── epigenerator/
│
├── results/
│   ├── sequencing_qc/
│   └── wet_lab_qc/
│
├── scripts/
│   ├── repo/
│   ├── slurm/
│   └── wet_lab/
│
└── tests/
```

### `data/raw/`

Contains source files in their original form. These files should not be modified during analysis.

Current raw wet-lab data include:

* sample manifests;
* DNA isolation records;
* Bioanalyzer submission forms; and
* Bioanalyzer reports.

DNA isolation records are organized by processing batch.

```text
data/raw/wet_lab/dna_isolation/
├── arizona/
│   └── DNAEXT_AZ_01/
└── la_county/
    ├── DNAEXT_LA_01/
    └── DNAEXT_LA_02/
```

The `arizona` designation reflects sample-processing and shipping origin and does **not** represent a separate Arizona firefighter cohort. These samples originated from the LACoFD study but were handled through an Arizona-associated workflow. This distinction is retained in the repository to preserve sample-processing origin.

### `data/interim/`

Contains cleaned or standardized datasets derived from the original source files.

Examples include datasets with:

* standardized variable names;
* harmonized timepoint definitions;
* standardized sample identifiers;
* extraction batch identifiers; and
* calculated DNA yield variables.

Intermediate files may be regenerated from the corresponding raw files and scripts.

### `data/processed/`

Contains datasets prepared for statistical analysis or downstream computational workflows.

For example:

```text
data/processed/wet_lab/dna_isolation/
```

contains analysis-ready DNA isolation datasets.

As WGBS processing is completed, final cytosine reports selected for downstream DNA methylation analyses will also be stored within the processed-data structure. Full pipeline working directories remain separate under `pipeline_runs/`.

### `data/metadata/`

Contains study metadata, data dictionaries, variable descriptions, and related supporting information.

### `pipeline_runs/`

Contains working directories generated during computational processing.

WGBS preprocessing is performed using[**epigenerator**](https://github.com/vhaghani26/epigenerator), a reproducible workflow for processing bisulfite sequencing data.

Epigenerator project repository: [**vhaghani26/epigenerator**](https://github.com/vhaghani26/epigenerator)

Individual processing runs are maintained separately:

```text
pipeline_runs/epigenerator/
├── EPI_PILOT_01/
├── EPI_AZ_01/
├── EPI_LA_01/
└── EPI_LA_02/
```

These directories may contain raw-sequence links, trimmed reads, alignment outputs, methylation extraction files, cytosine reports, quality-control files, and computational logs.

`pipeline_runs/` should therefore be considered a computational workspace rather than the permanent location of analysis-ready data.

### `results/`

Contains generated figures, tables, and quality-control summaries.

Current result categories include:

```text
results/
├── sequencing_qc/
└── wet_lab_qc/
```

Wet-lab DNA isolation results are further organized by processing batch where appropriate.

For example:

```text
results/wet_lab_qc/dna_isolation/la_county/
└── DNAEXT_LA_01/
```

### `scripts/`

Contains reproducible scripts used for data cleaning, quality control, plotting, and computational processing.

Scripts are grouped according to their purpose rather than by individual dataset whenever possible.

For example:

```text
scripts/
├── repo/
│   └── reorganize_repo.sh
│
├── slurm/
│
└── wet_lab/
    └── dna_isolation/
        ├── 01_clean_dna_isolation.R
        ├── 02_plot_dna_isolation_pre_post.R
        └── submit_rscript.sbatch
```

The DNA isolation scripts accept command-line arguments so the same workflow can be applied across multiple extraction batches without modifying hard-coded file paths.

Detailed execution instructions should be documented within README files located in the relevant script directories.

### `docs/`

Contains supporting documentation including study methods, wet-lab protocols, research plans, and reproducibility documentation.

## Naming conventions

Processing steps are assigned stable identifiers so that sample origin can be followed across wet-lab and computational workflows.

### DNA extraction batches

DNA extraction batches use:

```text
DNAEXT_<GROUP>_<BATCH>
```

Examples:

```text
DNAEXT_AZ_01
DNAEXT_LA_01
DNAEXT_LA_02
```

where:

* `DNAEXT` = DNA extraction;
* `AZ` = Arizona-associated processing;
* `LA` = Los Angeles County processing; and
* `01`, `02`, etc. = sequential processing batch.

The batch identifier describes the **DNA extraction event** and should not be interpreted as a sequencing or computational batch.

### Bioanalyzer runs

Bioanalyzer processing uses:

```text
BIOA_<GROUP>_<RUN>
```

Examples:

```text
BIOA_LA_01
BIOA_LA_02
BIOA_LA_03
```

Each Bioanalyzer run can contain the associated submission information and instrument-generated reports.

### Sequencing batches

Sequencing batches, when assigned, use:

```text
SEQ_<GROUP>_<BATCH>
```

For example:

```text
SEQ_LA_01
```

Sequencing batches are tracked independently from DNA extraction batches because samples from multiple extraction batches may be sequenced together.

### Epigenerator runs

Epigenerator processing runs use:

```text
EPI_<GROUP>_<RUN>
```

Examples:

```text
EPI_PILOT_01
EPI_AZ_01
EPI_LA_01
EPI_LA_02
```

where `EPI` identifies an Epigenerator processing run.

DNA extraction, Bioanalyzer, sequencing, and Epigenerator identifiers therefore describe different stages of sample processing:

```text
DNAEXT_LA_01   DNA extraction batch
BIOA_LA_01     Bioanalyzer run
SEQ_LA_01      Sequencing batch
EPI_LA_01      Epigenerator processing run
```

These identifiers should be maintained as separate metadata variables rather than inferred from one another.

## DNA isolation data processing

Cell-free DNA isolation data are processed using reusable R scripts located in:

```text
scripts/wet_lab/dna_isolation/
```

The general data flow is:

```text
Original extraction workbook
        │
        ▼
data/raw/
        │
        ▼
Cleaning and standardization
        │
        ▼
data/interim/
        │
        ▼
Analysis-ready dataset
        │
        ▼
data/processed/
        │
        ▼
QC figures and summaries
        │
        ▼
results/wet_lab_qc/
```

Current DNA isolation QC includes evaluation of:

* starting plasma volume;
* DNA concentration;
* total DNA yield;
* DNA yield normalized to starting plasma volume; and
* paired pre- and post-exposure measurements where available.

Detailed instructions for executing these scripts are maintained separately within the DNA isolation scripts directory.

## WGBS processing

Whole-genome bisulfite sequencing data are processed using the [**vhaghani26/epigenerator**](https://github.com/vhaghani26/epigenerator) workflow.

The full Epigenerator run is retained under:

```text
pipeline_runs/epigenerator/<RUN_ID>/
```

A typical run produces sequential outputs including:

```text
raw sequence links
      ↓
read trimming
      ↓
sequence screening
      ↓
alignment
      ↓
deduplication
      ↓
methylation extraction
      ↓
cytosine reports
      ↓
MultiQC
```

Large pipeline intermediates remain within the pipeline working directory. Final files selected for downstream analyses are copied or linked into the appropriate processed-data directory.

This distinction allows the complete computational workflow to be retained while keeping analysis-ready datasets separate from transient pipeline outputs.

## Project team

The LACOFD study is a collaborative project between the Los Angeles County Fire Department, the University of California, Davis, and collaborating investigators.

Current team members include:

| Team member                 | Role                                                         |
| --------------------------- | ------------------------------------------------------------ |
| **Thomas Sullivan**         | Los Angeles County Fire Department Co-Principal Investigator |
| **Janine LaSalle**          | UC Davis Co-Principal Investigator                           |
| **Shehnaz Hussain**         | UC Davis Co-Investigator                                     |
| **George Kuodza**           | UC Davis Postdoctoral Scholar                                |
| **Logan Williams**          | Former UC Davis Graduate Student                             |
| **Additional team members** | To be added                                                  |

The team list will be updated as project roles and contributors are finalized.

## Funding

This work is supported by the **California Firefighter Cancer Prevention and Research Program**, University of California Office of the President.

**Project:** *Examining longitudinal changes in DNA methylation in firefighters exposed to products of combustion*

**Award:** F01FF8767

**Fire service Co-Principal Investigator:** Thomas Sullivan, Los Angeles County Fire Department

**UC Co-Principal Investigator:** Janine LaSalle, University of California, Davis


Wet-lab processing, sample quality control, WGBS processing, and development of downstream DNA methylation analyses are ongoing. Repository organization, analysis scripts, and documentation will continue to be updated as additional samples and processing batches become available.
