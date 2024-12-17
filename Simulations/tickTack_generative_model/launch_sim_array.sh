#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=Generative_model_tickTack
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --output=tickTack_sim_%A_%a
#SBATCH --error=tickTack_sim_err_%A_%a
#SBATCH --array=1-3

module load R/4.4.1
echo $SLURM_ARRAY_TASK_ID

cd results

purity=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$1; exit }" params_config.txt)
coverage=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$2; exit }" params_config.txt)
n_clocks=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$3; exit }" params_config.txt)
n_events=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$4; exit }" params_config.txt)
epsilon=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$5; exit }" params_config.txt)
tolerance=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$6; exit }" params_config.txt)
max_attempts=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$7; exit }" params_config.txt)
seed=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print \$8; exit }" params_config.txt)


for (( i = 1 ; i <= 20 ; i += 1 )) ; do
  echo "i ha valore $i"
  seed=$i
  mkdir -p tickTack_sim_${purity}_${coverage}_${n_clocks}_${n_events}/${seed}
  mkdir -p tickTack_sim_${purity}_${coverage}_${n_clocks}_${n_events}/${seed}/results
  mkdir -p tickTack_sim_${purity}_${coverage}_${n_clocks}_${n_events}/${seed}/plots

  awk -F',' "NR==${SLURM_ARRAY_TASK_ID} { print ; exit }" params_config.txt >> tickTack_sim_${purity}_${coverage}_${n_clocks}_${n_events}/config
  
  echo $purity
  echo $coverage
  echo $seed
  echo $n_clocks
  echo $n_events
  
  Rscript ../scripts/simulate_wrapper.R ${purity} ${coverage} ${n_clocks} ${n_events} ${epsilon} ${tolerance} ${max_attempts} ${seed}
  
done
