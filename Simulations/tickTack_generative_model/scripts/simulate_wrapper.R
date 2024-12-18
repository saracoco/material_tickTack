#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
purity <- as.double(args[1])
coverage <- as.double(args[2])
n_clocks <- as.double(args[3])
n_events <- as.double(args[4])
epsilon <- as.double(args[5])
tolerance <- as.double(args[6])
max_attempts <- as.double(args[7])
seed <- as.double(args[8])


print("simulate wrapper")
print(purity) 
print(coverage)
print (n_clocks) 
print (n_events)
print(epsilon)
print(tolerance) 
print(max_attempts) 
print(seed) 


.libPaths(new="~/R/rstudio_v3/") 
library(MutationTimeR)
library(dplyr)
library(tickTack)
library(tibble)

source("../scripts/simulation_tickTack.R")

original_dir <- getwd()
cat(original_dir)
self_name = as.character(paste0("/tickTack_sim_",purity,"_",coverage,"_",n_clocks,"_",n_events,"/",seed))
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
                                  min_mutations_number = 3)


saveRDS(res$res_MutTime, "results/res_MutTime.rds")
saveRDS(res$res_tickTack, "results/res_tickTack.rds")
saveRDS(res$res_SingleTT, "results/res_SingleTT.rds")
saveRDS(res$compare_assignment, "results/compare_assignment.rds")

# ggsave("plots/plot_Muttime.png", plot= res$plot_MutTime)
ggsave("plots/plot_SingleTT.png", plot= res$plot_SingleTT, height=5, width=10)
ggsave("plots/plot_tickTack.png", plot= res$plot_tickTack, height=5, width=10)

  
