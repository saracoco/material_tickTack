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
#SBATCH --array=1-500

module load R
echo $SLURM_ARRAY_TASK_ID

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1000))

n_clocks=$(awk -F',' "NR==${LINE_NUMBER} { print \$1; exit }" params_config_reviews_1.txt)
n_events=$(awk -F',' "NR==${LINE_NUMBER} { print \$2; exit }" params_config_reviews_1.txt)
purity=$(awk -F',' "NR==${LINE_NUMBER} { print \$3; exit }" params_config_reviews_1.txt)
coverage=$(awk -F',' "NR==${LINE_NUMBER} { print \$4; exit }" params_config_reviews_1.txt)
epsilon=$(awk -F',' "NR==${LINE_NUMBER} { print \$5; exit }" params_config_reviews_1.txt)
tolerance=$(awk -F',' "NR==${LINE_NUMBER} { print \$6; exit }" params_config_reviews_1.txt)
max_attempts=$(awk -F',' "NR==${LINE_NUMBER} { print \$7; exit }" params_config_reviews_1.txt)
seed=$(awk -F',' "NR==${LINE_NUMBER} { print \$8; exit }" params_config_reviews_1.txt)
n_mutations=$(awk -F',' "NR==${LINE_NUMBER} { print \$9; exit }" params_config_reviews_1.txt)

awk -F',' "NR==${LINE_NUMBER} { print ; exit }" params_config_reviews_1.txt >> tickTack_sim_1_${n_clocks}_${n_events}_${purity}_${coverage}/config

# n_clocks=$(awk -F',' "NR==${LINE_NUMBER} { print \$1; exit }" params_config_reviews_2.txt)
# n_events=$(awk -F',' "NR==${LINE_NUMBER} { print \$2; exit }" params_config_reviews_2.txt)
# purity=$(awk -F',' "NR==${LINE_NUMBER} { print \$3; exit }" params_config_reviews_2.txt)
# coverage=$(awk -F',' "NR==${LINE_NUMBER} { print \$4; exit }" params_config_reviews_2.txt)
# epsilon=$(awk -F',' "NR==${LINE_NUMBER} { print \$5; exit }" params_config_reviews_2.txt)
# tolerance=$(awk -F',' "NR==${LINE_NUMBER} { print \$6; exit }" params_config_reviews_2.txt)
# max_attempts=$(awk -F',' "NR==${LINE_NUMBER} { print \$7; exit }" params_config_reviews_2.txt)
# seed=$(awk -F',' "NR==${LINE_NUMBER} { print \$8; exit }" params_config_reviews_2.txt)
# n_mutations=$(awk -F',' "NR==${LINE_NUMBER} { print \$9; exit }" params_config_reviews_2.txt)
# 
# awk -F',' "NR==${LINE_NUMBER} { print ; exit }" params_config_reviews_2.txt >> tickTack_sim_2_${n_clocks}_${n_events}_${purity}_${coverage}/config

# n_clocks=$(awk -F',' "NR==${LINE_NUMBER} { print \$1; exit }" params_config_reviews_3.txt)
# n_events=$(awk -F',' "NR==${LINE_NUMBER} { print \$2; exit }" params_config_reviews_3.txt)
# purity=$(awk -F',' "NR==${LINE_NUMBER} { print \$3; exit }" params_config_reviews_3.txt)
# coverage=$(awk -F',' "NR==${LINE_NUMBER} { print \$4; exit }" params_config_reviews_3.txt)
# epsilon=$(awk -F',' "NR==${LINE_NUMBER} { print \$5; exit }" params_config_reviews_3.txt)
# tolerance=$(awk -F',' "NR==${LINE_NUMBER} { print \$6; exit }" params_config_reviews_3.txt)
# max_attempts=$(awk -F',' "NR==${LINE_NUMBER} { print \$7; exit }" params_config_reviews_3.txt)
# seed=$(awk -F',' "NR==${LINE_NUMBER} { print \$8; exit }" params_config_reviews_3.txt)
# n_mutations=$(awk -F',' "NR==${LINE_NUMBER} { print \$9; exit }" params_config_reviews_3.txt)
# 
# awk -F',' "NR==${LINE_NUMBER} { print ; exit }" params_config_reviews_3.txt >> tickTack_sim_3_${n_clocks}_${n_events}_${purity}_${coverage}/config



echo $n_clocks
echo $n_events
echo $purity
echo $coverage
echo $epsilon
echo $tolerance
echo $max_attempts
echo $seed
echo $n_mutations

echo "NR==${LINE_NUMBER}"

Rscript 01_sim_and_fit.R ${n_clocks} ${n_events} ${purity} ${coverage} ${epsilon} ${tolerance} ${max_attempts} ${seed} ${n_mutations}
