#!/bin/bash
#SBATCH --job-name=PROJECT_NAME
#SBATCH --partition=production
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

echo "SLURM job started on Georges-MacBook-Pro-4.local at Thu Mar  5 00:16:56 PST 2026"
echo "Customize this to your project"

# Example:
# module load R
# Rscript analysis/run_analysis.R

echo "Done at Thu Mar  5 00:16:56 PST 2026"
