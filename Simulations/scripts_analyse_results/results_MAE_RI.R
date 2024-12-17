library(dplyr)
library(Metrics)
library(fossil)
library(future)
library(future.apply)

# Define the main directory containing all configurations
base_dir <- "."

# Function to analyze results for a single simulation folder
analyze_results <- function(sim_folder) {
  result_files <- list.files(sim_folder, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(result_files) > 0) {
    print(sim_folder)
    
    all_data <- lapply(result_files, readRDS)
    single_inference <- all_data[[3]]$summarized_results %>%
      mutate(chr = as.integer(gsub("chr", "", chr)))
    
    compare_tau <- single_inference %>%
      left_join(all_data[[1]], by = join_by(chr)) %>%
      mutate(
        tau_single = tau_mean,
        tau_real = taus,
        tau_hierarchical = clock_mean,
        tau_MutTime = time_MutTime
      )
    
    results <- tibble(
      MAE_single = mae(compare_tau$tau_real, compare_tau$tau_single),
      MAE_hierarchical = mae(compare_tau$tau_real, compare_tau$tau_hierarchical),
      MAE_MutTime = mae(compare_tau$tau_real, compare_tau$tau_MutTime),
      RI_real = rand.index(all_data[[1]]$clock_real_factor, all_data[[1]]$tickTack_estimate_factor),
      RI_real_adj = adj.rand.index(all_data[[1]]$clock_real_factor, all_data[[1]]$tickTack_estimate_factor)
    )
    
    return(results)
  } else {
    return(NULL)
  }
}

# Function to process each configuration directory
process_config <- function(config_dir) {
  index_dirs <- file.path(config_dir, 1:20)
  result_dirs <- file.path(index_dirs, "results")
  
  # Parallelize the result analysis
  res_iter <- future_lapply(result_dirs, analyze_results)
  res_iter <- bind_rows(res_iter, .id = "simulation")  # Efficient binding
  
  # Calculate column means
  column_means <- res_iter %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  return(column_means)
}

# Set up parallel backend
library(parallelly)
plan(multisession, workers = availableCores() - 1)

# Get all parameter configuration directories
parameter_dirs <- list.dirs(base_dir, recursive = FALSE)

# Apply the processing function in parallel
res_13 <- future_lapply(parameter_dirs[1:3], process_config)
res_46 <- future_lapply(parameter_dirs[4:6], process_config)
res_79 <- future_lapply(parameter_dirs[7:9], process_config)
res_10_12 <- future_lapply(parameter_dirs[10:12], process_config)
res_13_16 <- future_lapply(parameter_dirs[13:16], process_config)
res_17_20 <- future_lapply(parameter_dirs[17:20], process_config)


# Combine results into a single data frame
# final_result <- bind_rows(final_result, res_79, res_10_12, res_13_16, res_17_20)


saveRDS(final_result, "results_MAE_RI_0001tol.rds")














# Function to load and save .rds files for a single simulation folder
load_simulation_files <- function(sim_folder) {
  result_files <- list.files(sim_folder, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(result_files) > 0) {
    print(sim_folder)
    all_data <- lapply(result_files, readRDS)
    return(all_data)
  } else {
    return(NULL)
  }
}

# Function to process each configuration directory
process_config <- function(config_dir) {
  index_dirs <- file.path(config_dir, 1:20)
  result_dirs <- file.path(index_dirs, "results")
  
  # Load all .rds files for each simulation in parallel
  res_iter <- future_lapply(result_dirs, load_simulation_files)
  
  return(res_iter)
}

# Set up parallel backend
library(parallelly)
plan(multisession, workers = availableCores() - 1)

# Function to load all configurations and organize the results
load_simulation_results <- function(base_dir) {
  parameter_dirs <- list.dirs(base_dir, recursive = FALSE)
  
  # Apply the processing function in parallel
  results <- future_lapply(parameter_dirs, process_config)
  
  # Organize results in a structured list
  named_results <- setNames(results, basename(parameter_dirs))
  
  return(named_results)
}

# Load and organize all results
all_results <- load_simulation_results(base_dir)

# Example: Accessing .rds files for a specific configuration
config_1_files <- all_results[["config_1"]]
