#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=races_sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200gb
#SBATCH --time=48:00:00
#SBATCH --output=races_sim_%A_%a
#SBATCH --error=races_sim_err_%A_%a
#SBATCH --array=1-6

module load R/4.4.1


echo $SLURM_ARRAY_TASK_ID

purity=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$1; exit }" params_config.txt)
coverage=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$2; exit }" params_config.txt)
n_clocks=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$3; exit }" params_config.txt)
n_events=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$4; exit }" params_config.txt)
epsilon=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$5; exit }" params_config.txt)
seed=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$6; exit }" params_config.txt)

mkdir -p rRACES_sim_${purity}_${coverage}_${n_clocks}_${n_events}/${seed}
mkdir -p rRACES_sim_${purity}_${coverage}_${n_clocks}_${n_events}/${seed}/results
mkdir -p rRACES_sim_${purity}_${coverage}_${n_clocks}_${n_events}/${seed}/plots

awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print ; exit }" params_config.txt >> rRACES_sim_${purity}_${coverage}_${n_clocks}_${n_events}/config

echo $purity
echo $coverage
echo $seed
echo $n_clocks
echo $n_events

Rscript scripts/simulate_rRACES_wrapper.R ${purity} ${coverage} ${n_clocks} ${n_events} ${epsilon} ${seed}