.libPaths(new="~/R/rstudio_v3/") 

library(dplyr)
library(Metrics)
library(fossil)


# parallel
library(future)
library(future.apply)
library(parallel)



#check PCAWG tickTack results ############################################################à
results_time_1 <- list.files("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_tickTack_1")
results_time_parallel <- list.files("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_tickTack_parallel_2")

fittable_samples <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/fittable_samples.RDS")

count_files = function(a, counter) {
  n_files <- length(list.files(paste0("./",a,"/results")))
  if (n_files == 0){
    return (0)
  }
  else { return(1)}
}

get_files_success = function(a, counter) {
  n_files <- length(list.files(paste0("./",a,"/results")))
  if (n_files == 0){
    return(a)
  }
  else { return(NULL)}
}

result_files <- list.files()
# result_files_with_inference <- sum(lapply(result_files, count_files)%>%unlist())
# successes <- sum(lapply(result_files, count_files)%>%unlist())/length(list.files())
# successes_files <-(lapply(result_files, get_files_success)%>%unlist())

successes_parallel <- sum(lapply(results_time_parallel, count_files)%>%unlist()) #/length(list.files())
successes_files_parallel <-(lapply(results_time_parallel, get_files_success)%>%unlist())
# successes_1 <- sum(lapply(results_time_1, count_files)%>%unlist())/length(list.files())
# successes_files_1 <-(lapply(results_time_1, get_files_success)%>%unlist())

normalize_ids <- function(ids) {
  sub("data/clonal_analysis_PCAWG//", "", ids)
}
norm_fittable_samples <- normalize_ids(fittable_samples)
intersect(norm_fittable_samples, successes_files_1)



# list1 <- successes_files
# list2 <- successes_files_1
# successes_files1 <- list1[grepl("-0$", list1)]
# successes_files2 <- list2[grepl("-0$", list2)]
# normalize_ids <- function(ids) {
#   sub("-0$", "", ids)
# }
# norm_list1 <- normalize_ids(list1)
# norm_list2 <- normalize_ids(list2)
# common_ids <- intersect(norm_list1, norm_list2)
# 
# matches_with_suffix <- sapply(common_ids, function(id) {
#   in_list1_with_suffix <- paste0(id, "-0") %in% list1
#   in_list2_without_suffix <- id %in% list2
#   in_list2_with_suffix <- paste0(id, "-0") %in% list2
#   in_list1_without_suffix <- id %in% list1
#   if ((in_list1_with_suffix && in_list2_without_suffix) || 
#       (in_list2_with_suffix && in_list1_without_suffix)) {
#     return(id)
#   } else {
#     return(NULL)
#   }
# })
# matches_with_suffix <- matches_with_suffix[!sapply(matches_with_suffix, is.null)]
# matches_with_suffix
# exact_matches <- intersect(list1, list2)
# matches_with_suffix <- exact_matches[grepl("-0$", exact_matches)]


#######################################################################



# Define the main directory containing all configurations
base_dir <- "."



# Function to analyze results for a single simulation folder
analyze_results <- function(sim_folder) {
  result_files <- list.files(sim_folder, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(result_files) > 0) {
    print(sim_folder)
    
    all_data <- lapply(result_files, readRDS)
    # single_inference <- all_data[[3]]$summarized_results %>%
    #   mutate(chr = as.integer(gsub("chr", "", chr)))
    
    compare_tau <- all_data[[1]]%>%
      mutate(
        tau_single = singleTT,
        tau_real = real_clocks,
        tau_hierarchical = time_tickTack,
        tau_MutTime = time_MutTime
      )
    
    results <- tibble(
      MAE_single = mae(compare_tau$tau_real, compare_tau$tau_single),
      MAE_hierarchical = mae(compare_tau$tau_real, compare_tau$tau_hierarchical),
      MAE_MutTime = mae(compare_tau$tau_real, compare_tau$tau_MutTime),
      RI_real = rand.index(all_data[[1]]$clock_real_factor, all_data[[1]]$tickTack_estimate_factor),
      RI_real_adj = adj.rand.index(all_data[[1]]$clock_real_factor, all_data[[1]]$tickTack_estimate_factor),
    )
    
    
    MS <- res_tickTack$results_model_selection$model_selection_tibble
    #take best K according to scores
    # select results as res_tickTack$results$draws_and_summary[[best_k]]
    
    return(results)
  } else {
    return(NULL)
  }
}

# Function to process each configuration directory
process_config <- function(config_dir) {
  print(config_dir)
  index_dirs <- file.path(config_dir, 1:20)
  result_dirs <- file.path(index_dirs, "results")
  
  # Parallelize the result analysis
  res_iter <- future_lapply(result_dirs, analyze_results)
  res_iter <- bind_rows(res_iter, .id = "simulation")  # Efficient binding
  
  # Calculate column means
  column_means <- res_iter %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  column_means <- column_means %>% mutate(simulation = config_dir)
  
  return(column_means)
}

# # Set up parallel backend
# library(parallelly)
# plan(multisession, workers = availableCores() - 1)

# Get all parameter configuration directories
parameter_dirs <- list.dirs(base_dir, recursive = FALSE)

# # Apply the processing function in parallel
# res_13 <- future_lapply(parameter_dirs[1:3], process_config)
# res_46 <- future_lapply(parameter_dirs[4:6], process_config)
# res_79 <- future_lapply(parameter_dirs[7:9], process_config)
# res_10_12 <- future_lapply(parameter_dirs[10:12], process_config)
# res_13_16 <- future_lapply(parameter_dirs[13:16], process_config)
# res_17_20 <- future_lapply(parameter_dirs[17:20], process_config)


# Combine results into a single data frame
# final_result <- bind_rows(final_result, res_79, res_10_12, res_13_16, res_17_20)






final_result <- parallel::mclapply(parameter_dirs, process_config, mc.cores = (parallel::detectCores()-1), mc.allow.recursive = TRUE)




saveRDS(final_result, "results_MAE_RI_0001tol.rds")

library(dplyr)
# capisci problema parallel computation 
results_MAE_RI_0001tol[c(23, 24, 25, 28, 29, 30)] = NULL
final_result <- bind_rows(results_MAE_RI_0001tol)

library(stringr)


final_result <- final_result %>%
  mutate(
    purity = str_extract(simulation, "\\d+\\.\\d+"),
    coverage = str_match(simulation, "\\d+\\.\\d+_(\\d+)")[, 2],
    n_clocks = str_match(simulation, "\\d+\\.\\d+_\\d+_(\\d+)")[, 2],
    n_events = str_match(simulation, "\\d+\\.\\d+_\\d+_\\d+_(\\d+)")[, 2]
  )
# Converting extracted values to numeric
final_result <- final_result %>%
  mutate(across(c(purity, coverage, n_clocks, n_events), as.numeric))

saveRDS(final_result, "tibble_results_MAE_RI_0001tol.rds")







# visualize results 

library(ggplot2)
library(dplyr)

# Example: Grouping by purity and visualizing MAE scores
final_result %>%
  mutate(purity = as.factor(purity)) %>% # Convert purity to a factor for grouping
  ggplot(aes(x = purity, y = MAE_single)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "MAE_single Scores by Purity",
    x = "Purity",
    y = "MAE_single"
  ) +
  theme_minimal()

# Example: Facet by n_clocks and visualize MAE scores
final_result %>%
  mutate(n_clocks = as.factor(n_clocks)) %>% # Convert n_clocks to a factor for faceting
  ggplot(aes(x = as.factor(coverage), y = MAE_hierarchical, fill = n_clocks)) +
  geom_boxplot() +
  facet_wrap(~n_clocks) +
  labs(
    title = "MAE_hierarchical by Coverage (Faceted by n_clocks)",
    x = "Coverage",
    y = "MAE_hierarchical",
    fill = "n_clocks"
  ) +
  theme_minimal()


# plot all mae

library(tidyr)

# Reshape data into long format
mae_long <- final_result %>%
  pivot_longer(
    cols = starts_with("MAE"), # Select all columns starting with "MAE"
    names_to = "MAE_type",     # New column for MAE types
    values_to = "MAE_value"    # New column for MAE values
  )

# Plot all MAE scores
ggplot(mae_long, aes(x = as.factor(purity), y = MAE_value, fill = MAE_type)) +
  geom_boxplot() +
  facet_wrap(~MAE_type, scales = "free_y") + # Facet by MAE type
  labs(
    title = "Distribution of MAE Scores by Purity",
    x = "Purity",
    y = "MAE Value",
    fill = "MAE Type"
  ) +
  theme_minimal()



###########################################
# fix coverage and purity to compare them 


# Filter data for a specific coverage value
filtered_df <- final_result %>% filter(coverage == 100)  # Example: Fix coverage at 100
# Reshape the filtered dataset for plotting
mae_long_filtered <- filtered_df %>%
  pivot_longer(
    cols = starts_with("MAE"), # Select all MAE columns
    names_to = "MAE_type",     # Column for MAE type
    values_to = "MAE_value"    # Column for MAE values
  )

# Plot MAE scores grouped by purity for a fixed coverage
ggplot(mae_long_filtered, aes(x = as.factor(purity), y = MAE_value, fill = MAE_type)) +
  geom_boxplot() +
  facet_wrap(~MAE_type, scales = "free_y") +
  labs(
    title = "MAE Scores by Purity (Fixed Coverage = 100)",
    x = "Purity",
    y = "MAE Value",
    fill = "MAE Type"
  ) +
  theme_minimal()








filtered_df <- final_result %>% filter(coverage == 100)  # Fix coverage at 100

# Reshape data for plotting
mae_long_filtered <- filtered_df %>%
  pivot_longer(
    cols = starts_with("MAE"), # Select all MAE columns
    names_to = "MAE_type",     # Column for MAE type
    values_to = "MAE_value"    # Column for MAE values
  )



# Plot MAE scores by purity and facet by n_clocks
ggplot(mae_long_filtered, aes(x = as.factor(purity), y = MAE_value, fill = MAE_type)) +
  geom_boxplot() +
  facet_wrap(~n_clocks, scales = "free_y") +  # Facet by n_clocks
  labs(
    title = "MAE Scores by Purity (Fixed Coverage = 100, Grouped by Clocks)",
    x = "Purity",
    y = "MAE Value",
    fill = "MAE Type"
  ) +
  theme_minimal()


ggplot(mae_long_filtered, aes(x = as.factor(purity), y = MAE_value, fill = as.factor(n_clocks))) +
  geom_boxplot(position = "dodge") +
  facet_wrap(~MAE_type, scales = "free_y") +  # Facet by MAE type
  labs(
    title = "MAE Scores by Purity and Clocks (Fixed Coverage = 100)",
    x = "Purity",
    y = "MAE Value",
    fill = "Number of Clocks"
  ) +
  theme_minimal()


# Filter data for a specific purity value
filtered_df <- final_result %>% filter(purity == 0.1)  # Fix purity at 0.1

# Reshape and plot as before, grouping by coverage and faceting by n_clocks
ggplot(mae_long_filtered, aes(x = as.factor(coverage), y = MAE_value, fill = as.factor(n_clocks))) +
  geom_boxplot(position = "dodge") +
  facet_wrap(~MAE_type, scales = "free_y") +  # Facet by MAE type
  labs(
    title = "MAE Scores by Coverage and Clocks (Fixed Purity = 0.1)",
    x = "Coverage",
    y = "MAE Value",
    fill = "Number of Clocks"
  ) +
  theme_minimal()


###############################################








# compare model selection scores results 






















