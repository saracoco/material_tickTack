#!/bin/bash
#SBATCH --job-name=PCAWG_timing
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00
#SBATCH --partition=THIN
#SBATCH --mem=100gb
#SBATCH --output=PCAWG_%j.out


module load R/4.4.1

R CMD BATCH 00_smoothing.R

module purge
 
