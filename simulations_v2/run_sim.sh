#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=sim_tickTack
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --time=48:00:00
#SBATCH --output=out/%a
#SBATCH --error=err/%a
#SBATCH --array=1-1500

module load R
echo $SLURM_ARRAY_TASK_ID

n_clocks=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$1; exit }" params_config.txt)
n_events=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$2; exit }" params_config.txt)
purity=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$3; exit }" params_config.txt)
coverage=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$4; exit }" params_config.txt)
epsilon=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$5; exit }" params_config.txt)
tolerance=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$6; exit }" params_config.txt)
max_attempts=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$7; exit }" params_config.txt)
seed=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$8; exit }" params_config.txt)
n_mutations=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$9; exit }" params_config.txt)

awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print ; exit }" params_config_reviews_1.txt >> tickTack_sim_1_${n_clocks}_${n_events}_${purity}_${coverage}/config

echo $n_clocks
echo $n_events
echo $purity
echo $coverage
echo $epsilon
echo $tolerance
echo $max_attempts
echo $seed
echo $n_mutations

echo "NR==${SLURM_ARRAY_TASK_ID}"

Rscript 01_sim_and_fit.R ${n_clocks} ${n_events} ${purity} ${coverage} ${epsilon} ${tolerance} ${max_attempts} ${seed} ${n_mutations}
