#!/bin/bash
#SBATCH --job-name=smooth
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:20:00
#SBATCH --partition=THIN
#SBATCH --mem=100gb
#SBATCH --output=smooth_%j.out


module load R/4.4.1

R CMD BATCH 00_samples_statistics.R

module purge

