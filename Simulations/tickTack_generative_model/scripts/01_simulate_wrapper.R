#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
n_clocks <- as.double(args[1])
n_events <- as.double(args[2])
purity <- as.double(args[3])
coverage <- as.double(args[4])
epsilon <- as.double(args[5])
tolerance <- as.double(args[6])
max_attempts <- as.double(args[7])
seed <- as.double(args[8])

print("01_simulate_wrapper")
print (paste0("n_clocks: ",n_clocks)) 
print (paste0("n_events: ",n_events))
print(paste0("purity: ",purity)) 
print(paste0("coverage: ",coverage))
print(paste0("epsilon: ",epsilon))
print(paste0("tolerance: ",tolerance)) 
print(paste0("max_attempts: ",max_attempts)) 
print(paste0("seed: ",seed)) 


.libPaths(new="~/R/rstudio_v3/") 
library(MutationTimeR)
library(dplyr)
library(tickTack)
library(tibble)
source("../scripts/02_simulation_tickTack.R")

print("01_simulate_wrapper")

original_dir <- getwd()
cat(original_dir)
self_name = as.character(paste0("/tickTack_sim_",n_clocks,"_",n_events,"_",purity,"_",coverage,"/",seed))
new_dir = paste0(original_dir,self_name)
cat(new_dir)
setwd(new_dir)


# al variare del seed dovrei farlo per 20 simulazioni quindi 20 seed diversi
res <- simulation_tickTack(n_clocks=n_clocks, 
                                  n_events=n_events, 
                                  purity=purity, 
                                  coverage=coverage, 
                                  epsilon=epsilon, 
                                  seed = seed, 
                                  tolerance = tolerance,
                                  max_attempts = max_attempts,
                                  INIT = TRUE,
                                  min_mutations_number = 4)


model_selection_tibble <- res$res_tickTack$results_model_selection$model_selection_tibble
unique_K <- res$res_tickTack$results_model_selection$model_selection_tibble$K
df_summary_tickTack = data.frame()
for (i in unique_K){
  fit_single <- res$res_tickTack$results$draws_and_summary[[as.character(i)]]
  summary <- fit_single$summarized_results
  summary <- summary %>% mutate(n_components = i)
  df_summary_tickTack <- bind_rows(df_summary_tickTack, summary)
}

results_summary <- list(df_summary=df_summary_tickTack, model_selection_tibble = model_selection_tibble, compare_assignment = res$compare_assignment )
saveRDS( results_summary , paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/res_", n_clocks,"_",n_events,"_",purity,"_",coverage,"_",seed, ".rds"))


saveRDS(res$res_MutTime, "results/res_MutTime.rds")
saveRDS(res$res_tickTack, "results/res_tickTack.rds")
saveRDS(res$res_SingleTT, "results/res_SingleTT.rds")
saveRDS(res$compare_assignment, "results/compare_assignment.rds")

# ggsave("plots/plot_Muttime.png", plot= res$plot_MutTime)
ggsave("plots/plot_SingleTT.png", plot= res$plot_SingleTT, height=5, width=10)
ggsave("plots/plot_tickTack.png", plot= res$plot_tickTack, height=5, width=10)

  
