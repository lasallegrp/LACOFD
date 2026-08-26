#!/bin/bash
#SBATCH --job-name=LACoFD_DNA_Isolation
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G #700
#SBATCH --account=lasallegrp
#SBATCH --output=00_std_err_logs/dna_report.out
#SBATCH --error=00_std_err_logs/dna_report.err
#SBATCH --chdir=/quobyte/lasallegrp/projects/LACOFD

set -euo pipefail

#This is an old env where i have Rscript installed. I will need to update this to a newer env with Rscript and other dependencies installed.
conda activate Comethyl_v1.3.0

#La County DNA Isolation Report Batch 1
Rscript scripts/wet_lab/dna_isolation/01_clean_dna_isolation.R \
  --input data/raw/wet_lab/dna_isolation/la_county/DNAEXT_LA_01/LACoFD_LaSalle_Sullivan_DNA_Isolation_01_29June2026_GK.xlsx \
  --output-dir data/processed/wet_lab/dna_isolation/la_county/DNAEXT_LA_01 \
  --site la_county \
  --batch DNAEXT_LA_01 \
  --operator george_kuodza

Rscript scripts/wet_lab/dna_isolation/02_plot_dna_isolation_pre_post.R \
  --input data/processed/wet_lab/dna_isolation/la_county/DNAEXT_LA_01/DNAEXT_LA_01_dna_isolation_analysis_ready.tsv \
  --output-dir results/wet_lab_qc/dna_isolation/la_county/DNAEXT_LA_01 \
  --label "LA County — DNA extraction batch 1" \
  --prefix DNAEXT_LA_01
  
echo "Done $(Sys.Date())"
