#!/bin/bash
#SBATCH --job-name=PROJECT_NAME_smk
#SBATCH --partition=production
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Example usage:
# sbatch scripts/slurm/submit_snakemake.sbatch

snakemake --profile workflow/profiles/slurm -j 100
