#!/bin/bash
#SBATCH --job-name=fitting_timing
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --array=1-117
#SBATCH --output=results/log/output_%A_%a.out
#SBATCH --error=results/log/error_%A_%a.err

module load R/4.4.1

Rscript 01_fit_chunk.R ${SLURM_ARRAY_TASK_ID} 117
