#!/bin/bash

# Define the parameter ranges
# 1
# purity_values=(0.1 0.2 0.3 0.6 0.9)
# coverage_values=(100 50 20 10 5)
# n_clocks_values=(1 2 3 4 5)
# n_events_values=(5 10 20)
# n_mutations_values=(10 20 50 100)
# tolerance=0.01
# max_attempts=2
# seed=123
# 2
# purity_values=(0.1 0.2 0.3 0.6 0.9)
# coverage_values=(100 50 20 10 5)
# n_clocks_values=(6 8 10 15)
# n_events_values=(15 20 30 40)
# n_mutations_values=(10 20 50 100)
# tolerance=0.01
# max_attempts=2
# seed=123
# 3
purity_values=(0.1 0.2 0.3 0.6 0.9)
coverage_values=(100 50 20 10 5)
n_clocks_values=(6 8 10 15 20)
n_events_values=(6 8 10 15 20)
n_mutations_values=(10 20 50 100)
tolerance=0.01
max_attempts=2
seed=123

# Create the file
# 1
# echo "n_clocks,n_events,purity,coverage,epsilon,tolerance,max_attempts,seed,n_mutations" > params_config_reviews_1.txt
# 2
# echo "n_clocks,n_events,purity,coverage,epsilon,tolerance,max_attempts,seed,n_mutations" > params_config_reviews_2.txt
# 3
# echo "n_clocks,n_events,purity,coverage,epsilon,tolerance,max_attempts,seed,n_mutations" > params_config_reviews_3.txt


# Generate all parameter combinations
# 1, 2
# for purity in "${purity_values[@]}"; do
#     for coverage in "${coverage_values[@]}"; do
#         for n_clocks in "${n_clocks_values[@]}"; do
#             for n_events in "${n_events_values[@]}"; do
#                 for n_mutations in "${n_mutations_values[@]}"; do
#                     # Set epsilon based on n_clocks
#                     if [[ "$n_clocks" -le 3 ]]; then
#                         epsilon=0.2
#                     elif [[ "$n_clocks" -eq 4 ]]; then
#                         epsilon=0.15
#                     elif [[ "$n_clocks" -eq 5 ]]; then
#                         epsilon=0.1
#                     else 
#                         epsilon=0.06
#                     fi
# 
#                     # Write to file
#                     echo "$n_clocks,$n_events,$purity,$coverage,$epsilon,$tolerance,$max_attempts,$seed,$n_mutations" >> params_config_reviews_2.txt
#                 done
#             done
#         done
#     done
# done
# 
# echo "Generated params_config.txt with all combinations."


# 3 
# Generate all parameter combinations
for purity in "${purity_values[@]}"; do
    for coverage in "${coverage_values[@]}"; do
        for n_clocks in "${n_clocks_values[@]}"; do
            for n_mutations in "${n_mutations_values[@]}"; do
                n_events=$n_clocks  # Ensure n_events == n_clocks
                epsilon=0.04
                # Write to file
                echo "$n_clocks,$n_events,$purity,$coverage,$epsilon,$tolerance,$max_attempts,$seed,$n_mutations" >> params_config_reviews_3.txt
            done
        done
    done
done

echo "Generated params_config.txt with all combinations."
