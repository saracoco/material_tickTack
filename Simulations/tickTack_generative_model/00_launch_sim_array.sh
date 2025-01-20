#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=Generative_model_tickTack
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=tickTack_sim_%A_%a
#SBATCH --error=tickTack_sim_err_%A_%a
#SBATCH --array=1-120

module load R/4.4.1
echo $SLURM_ARRAY_TASK_ID

cd results

n_clocks=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$1; exit }" params_config.txt)
n_events=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$2; exit }" params_config.txt)
purity=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$3; exit }" params_config.txt)
coverage=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$4; exit }" params_config.txt)
epsilon=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$5; exit }" params_config.txt)
tolerance=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$6; exit }" params_config.txt)
max_attempts=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$7; exit }" params_config.txt)
seed=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$8; exit }" params_config.txt)


for (( i = 1 ; i <= 20 ; i += 1 )) ; do
  echo "i ha valore $i"
  seed=$i
  mkdir -p tickTack_sim_${n_clocks}_${n_events}_${purity}_${coverage}/${seed}
  mkdir -p tickTack_sim_${n_clocks}_${n_events}_${purity}_${coverage}/${seed}/results
  mkdir -p tickTack_sim_${n_clocks}_${n_events}_${purity}_${coverage}/${seed}/plots

  awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print ; exit }" params_config.txt >> tickTack_sim_${n_clocks}_${n_events}_${purity}_${coverage}/config
  
  echo $n_clocks
  echo $n_events
  echo $purity
  echo $coverage
  echo $epsilon
  echo $tolerance
  echo $max_attempts
  echo $seed
  
  Rscript ../scripts/01_simulate_wrapper.R ${n_clocks} ${n_events} ${purity} ${coverage} ${epsilon} ${tolerance} ${max_attempts} ${seed}
  
done

