#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=clustering
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --time=24:00:00
#SBATCH --output=out/clustering
#SBATCH --error=err/clustering

module load R
# Rscript 02_parse_results.R
# Rscript 02a_cluster_other_methods.R
Rscript 02b_cluster_other_methods_alternative_algorithms.R
