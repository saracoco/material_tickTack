#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
n_clocks <- as.double(args[1])
n_events <- as.double(args[2])
pi = purity <- as.double(args[3])
coverage <- as.double(args[4])
epsilon <- as.double(args[5])
tolerance <- as.double(args[6])
max_attempts <- as.double(args[7])
seed <- as.double(args[8])
n_mutations <- as.numeric(args[9])
set.seed(seed)

source("utils.R")

MIN_MUTATIONS = 4

print (paste0("n_clocks: ",n_clocks)) 
print (paste0("n_events: ",n_events))
print(paste0("purity: ",purity)) 
print(paste0("coverage: ",coverage))
print(paste0("epsilon: ",epsilon))
print(paste0("tolerance: ",tolerance)) 
print(paste0("max_attempts: ",max_attempts)) 
print(paste0("seed: ",seed))

main_dir = paste0("results/sim_", n_clocks, "_", n_events, "_", purity, "_", coverage, "_", n_mutations, "/")
if (!dir.exists(main_dir)) {
  dir.create(main_dir)  
}

for (i.iter in 1:10) {
  sub_dir = paste0(main_dir, i.iter)
  dir.create(sub_dir)
  
  sim = simulate_dataset(n_events, n_clocks, n_mutations, pi, coverage, sigma_tau = .01, min_dist = epsilon)
  
  res_AmpTimeR <- fit_AmpTimeR(sim)
  res_MutTime <- fit_MutTimeR(sim, pi)
  res_tickTack_single <- fit_tickTack_single(sim, pi, MIN_MUTATIONS)
  res_tickTack_h <- fit_tickTack_h(sim, pi, MIN_MUTATIONS, INIT = TRUE, tolerance = tolerance, max_attempts = max_attempts)
  
  # Merge results
  merged_res = dplyr::tibble(
    segment_idx = res_AmpTimeR$segment_idx,
    true_tau = sim$true_taus,
    tau_AmpTimeR = res_AmpTimeR$tau,
    tau_MutTimeR = res_MutTime$cn_timed$time,
    tau_tickTack = res_tickTack_single$summarized_results$tau_mean,
    tau_tickTack_h = res_tickTack_h$results_model_selection$best_fit$summarized_results$clock_mean
  )
  
  saveRDS(sim, paste0(sub_dir, "/sim.rds"))
  saveRDS(merged_res, paste0(sub_dir, "/merged_res.rds"))
  saveRDS(res_AmpTimeR, paste0(sub_dir, "/res_AmpTimeR.rds"))
  saveRDS(res_MutTime, paste0(sub_dir, "/res_MutTime.rds"))
  saveRDS(res_tickTack_single, paste0(sub_dir, "/res_tickTack_single.rds"))
  saveRDS(res_tickTack_h, paste0(sub_dir, "/res_tickTack_h.rds"))
}







