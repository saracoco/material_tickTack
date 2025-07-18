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
source("utils_GRITIC.R")

MIN_MUTATIONS = 10 # write that is the minimum requested by gritic (check if it is like this)

print (paste0("n_clocks: ",n_clocks)) 
print (paste0("n_events: ",n_events))
print(paste0("purity: ",purity)) 
print(paste0("coverage: ",coverage))
print(paste0("epsilon: ",epsilon))
print(paste0("tolerance: ",tolerance)) 
print(paste0("max_attempts: ",max_attempts)) 
print(paste0("seed: ",seed))

main_dir = paste0("results_gritic/sim_", n_clocks, "_", n_events, "_", purity, "_", coverage, "_", n_mutations, "/")
if (!dir.exists(main_dir)) {
  dir.create(main_dir)  
}

for (i.iter in 1:5) {
  sub_dir = paste0(main_dir, i.iter)
  dir.create(sub_dir)
  
  sim = simulate_dataset(n_events, n_clocks, n_mutations, pi, coverage, sigma_tau = .01, min_dist = epsilon)
  
  # Define the log file
  error_log <- file(paste0(sub_dir, "/error_log_1.txt"), open = "wt")
  
  # Function to safely run and catch errors
  safe_run <- function(expr, name) {
    tryCatch(
      expr,
      error = function(e) {
        msg <- paste(Sys.time(), "-", name, "failed with error:", e$message, "\n")
        writeLines(msg, error_log)
        return(NULL)
      }
    )
  }
  
  # Run fits and catch errors
  mult_path <- paste0(sub_dir,"/sample_1_dpclust_info")
  get_multiplicities(sim, purity, mult_path, sub_dir)
  res_AmpTimeR <- safe_run(fit_AmpTimeR(sim,mult_path), "fit_AmpTimeR")
  res_MutTime  <- safe_run(fit_MutTimeR(sim, pi), "fit_MutTimeR")
  res_tickTack_single <- safe_run(fit_tickTack_single(sim, pi, MIN_MUTATIONS), "fit_tickTack_single")
  res_tickTack_h <- safe_run(fit_tickTack_h(sim, pi, MIN_MUTATIONS, INIT = TRUE, tolerance = tolerance, max_attempts = max_attempts), "fit_tickTack_h")
  
  
  # Save simulation results
  saveRDS(sim, paste0(sub_dir, "/sim.rds"))
  if (!is.null(res_AmpTimeR)) saveRDS(res_AmpTimeR, paste0(sub_dir, "/res_AmpTimeR.rds"))
  if (!is.null(res_MutTime)) saveRDS(res_MutTime, paste0(sub_dir, "/res_MutTime.rds"))
  if (!is.null(res_tickTack_single)) saveRDS(res_tickTack_single, paste0(sub_dir, "/res_tickTack_single.rds"))
  if (!is.null(res_tickTack_h)) saveRDS(res_tickTack_h, paste0(sub_dir, "/res_tickTack_h.rds"))
  
  
  # Check if any result is NULL before merging
  if (!is.null(res_AmpTimeR) && !is.null(res_MutTime) && !is.null(res_tickTack_single) && !is.null(res_tickTack_h)) {
    merged_res <- dplyr::tibble(
      segment_idx = res_AmpTimeR$segment_idx,
      true_tau = sim$true_taus,
      true_tau_cluster = sim$taus_clust,
      tau_AmpTimeR = res_AmpTimeR$tau,
      tau_MutTimeR = res_MutTime$cn_timed$time,
      tau_tickTack = res_tickTack_single$summarized_results$tau_mean,
      tau_tickTack_h = res_tickTack_h$results_model_selection$best_fit$summarized_results$clock_mean)  
      # tau_gritic = res_gritic )  # extract the right results && !is.null(res_gritic)
    
    saveRDS(merged_res, paste0(sub_dir, "/merged_res.rds"))
  }
  
  close(error_log)
  
}