#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
purity <- as.double(args[1])
coverage <- as.double(args[2])
n_clocks <- as.double(args[3])
n_events <- as.double(args[4])
epsilon <- as.double(args[5])
seed <- as.double(args[6])
# tolerance <- as.double(args[7])
# max_attempts <- as.double(args[8])

.libPaths(new="~/R/rstudio_v3/") 

source("./scripts/simulate_rRACES.R")
source("./scripts/simulate_tissue.R")
source("./scripts/place_mutations.R")
source("./scripts/get_sequencing.R")
source("./scripts/plot_races.R")
source("./scripts/results_races.R")

#source("./scripts/simulate_rRACES.R")

original_dir <- (getwd())
cat(original_dir)

self_name = as.character(paste0("/rRACES_sim_",purity,"_",coverage,"_",n_clocks,"_",n_events,"/",seed))
new_dir = paste0(original_dir,self_name)
cat(new_dir)
setwd(new_dir)

# seed = 123
# coverage = 10
# purity = 0.4
# ####
# n_clocks = 3
# n_events = 8
# simulation_time = 500
# epsilon = 0.20
# chromosomes = c("5","6","10","12","13")
# type_length = c("long")

# forest <- load_samples_forest("samples_forest.sff")
# phylo_forest <- load_phylogenetic_forest("results/phylo_forest.sff")
# seq_results <- readRDS(paste0("results/seq_final_", coverage, "X.RDS"))
# seq_results_old <- readRDS(paste0("results/seq_", coverage, "X.RDS"))
# 
# 


simulate_rRACES(purity = purity, coverage = coverage, seed = seed, K_clocks = n_clocks, n_events = n_events, epsilon = epsilon, type_length = c("long"), chromosomes = c("5","6","10","12","13")
)

# simulate_rRACES(purity = 0.90, coverage = 70, seed = 123, K_clocks = 3, n_events = 8)
