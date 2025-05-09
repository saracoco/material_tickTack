#!/bin/bash
#SBATCH --job-name=PCAWG_BRCA_timing
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --array=1-117
#SBATCH --output=log/output_%A_%a.out
#SBATCH --error=log/error_%A_%a.err

module load R/4.4.1

Rscript 01_PCAWG_fit_chunk.R ${SLURM_ARRAY_TASK_ID} 117
