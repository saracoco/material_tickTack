#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=clustering
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=out/clusterin
#SBATCH --error=err/clustering

module load R
Rscript 02_parse_results.R
Rscript 02a_cluster_other_methods.R
