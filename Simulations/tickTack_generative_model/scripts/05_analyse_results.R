.libPaths(new="~/R/rstudio_v3/") 

library(dplyr)
library(Metrics)
library(fossil)

# parallel
library(future)
library(future.apply)
library(parallel)

# print(parallel::detectCores()-1)


assign_clusters <- function(values) {
  as.integer(as.factor(values)) # Group values into unique clusters
}

parse_simulation_name <- function(name) {
  simulation_details <- sub("res_", "", name)
  simulation_details <- sub(".rds", "", simulation_details)
  details <- strsplit(simulation_details, "_")[[1]]
  list(
    n_clocks = as.integer(details[1]),
    n_events = as.integer(details[2]),
    purity = as.numeric(details[3]),
    coverage = as.integer(details[4]),
    seed = as.integer(details[5])
  )
}


my_function <- function (single_inference){
  
  results <- single_inference$df_summary %>%
    inner_join(single_inference$compare_assignment, by = c("segment_id" = "segment_original_indx")) %>%
    group_by(n_components) %>%
    summarise(
      MSE = mse(clock_mean, real_clocks),
      RI = rand.index(
        as.numeric(assign_clusters(real_clocks)),
        as.numeric(tickTack_estimate_factor)
      ),
      # ARI = adj.rand.index(
      #   as.numeric(assign_clusters(real_clocks)),
      #   as.numeric(tickTack_estimate_factor)
      # ),
      .groups = "drop"
    )
  
  return(results)
  
}



calculate_mae <- function(compare_assignment) {
  compare_assignment %>%
    summarise(
      MAE_tickTack = mae(time_tickTack, real_clocks),
      MAE_singleTT = mae(singleTT, real_clocks),
      MAE_MutTimeR = mae(time_MutTime, real_clocks)
    )
}



check_unique_values <- function(compare_assignment) {
  compare_assignment %>%
    summarise(
      n_clocks_real = n_distinct(real_clocks),
      n_clocks_tickTack = n_distinct(time_tickTack),
      difference = n_clocks_real - n_clocks_tickTack
    )
}



find_best_models <- function(model_selection_tibble) {
  model_selection_tibble %>%
    summarise(
      best_BIC = K[which.min(BIC)],
      best_Log_lik = K[which.max(Log_lik)],
      best_ICL = K[which.min(ICL)],
      best_AIC = K[which.min(AIC)],
      best_LOO = K[which.min(LOO)]
    )
}





setwd("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary")
files <- list.files()


# Get the information about the simulation for each simulation in the summary_results
variable_names <- files

simulation_info <- lapply(variable_names, parse_simulation_name)
names(simulation_info) <- variable_names

# simulation_info # list with keys as the name of the files with the information about the simulation 






my_function_2 <- function(a){
  single_inference <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/",a))
  simulation_info_single <- simulation_info[[a]]
  
  MS_scores_per_component <- my_function(single_inference)
  
  
  mae_ri_results <- calculate_mae(single_inference$compare_assignment)
  
  check_components_difference <- check_unique_values(single_inference$compare_assignment)
  
  best_models_MS <- find_best_models(single_inference$model_selection_tibble)
  
  returned_item = list(MS_scores_per_component = MS_scores_per_component, mae_ri_results=mae_ri_results, check_components_difference=check_components_difference, best_models_MS = best_models_MS )
}

files <- files
res_single_inference <- lapply(files, my_function_2) #parallel::mclapply  , mc.cores = (parallel::detectCores()-1), mc.allow.recursive = TRUE
names(res_single_inference) <- variable_names
# res_single_inference[[a]]


##################################################

aggregated_results <- list()
# indx <- 1
for (indx in seq_along(res_single_inference)) {
  
  sim_info <- simulation_info[[indx]]
  n_clocks <- sim_info$n_clocks
  n_events <- sim_info$n_events
  purity <- sim_info$purity
  coverage <- sim_info$coverage
  seed <- sim_info$seed
  
  ms_scores <- res_single_inference[[indx]]$MS_scores_per_component
  mae_ri_results <- res_single_inference[[indx]]$mae_ri_results
  check_components_difference <- res_single_inference[[indx]]$check_components_difference
  best_models_MS <- res_single_inference[[indx]]$best_models_MS
  
  ms_scores_augmented <- ms_scores %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed
    )
  
  mae_ri_augmented <- mae_ri_results %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed
    )
  
  check_diff_augmented <- check_components_difference %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed
    )
  
  best_models_augmented <- best_models_MS %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed
    )
  
  aggregated_results[[indx]] <- list(
    MS_scores = ms_scores_augmented,
    MAE_RI = mae_ri_augmented,
    Component_Check = check_diff_augmented,
    Best_Models = best_models_augmented
  )
}


final_results <- list(
  MS_scores = bind_rows(lapply(aggregated_results, `[[`, "MS_scores")),
  MAE_RI = bind_rows(lapply(aggregated_results, `[[`, "MAE_RI")),
  Component_Check = bind_rows(lapply(aggregated_results, `[[`, "Component_Check")),
  Best_Models = bind_rows(lapply(aggregated_results, `[[`, "Best_Models"))
)



final_results_combined <- final_results$MS_scores %>%
  left_join(final_results$MAE_RI, by = c("n_clocks", "n_events", "purity", "coverage", "seed")) %>%
  left_join(final_results$Component_Check, by = c("n_clocks", "n_events", "purity", "coverage", "seed")) %>%
  left_join(final_results$Best_Models, by = c("n_clocks", "n_events", "purity", "coverage", "seed"))

############################################à

saveRDS(final_results_combined, "final_results_combined.rds")





  

  
  
.libPaths(new="~/R/rstudio_v3/") 

library(dplyr)
library(Metrics)
library(fossil)
library(tidyr)
library(ggplot2)

# parallel
library(future)
library(future.apply)
library(parallel)

# print(parallel::detectCores()-1)

print(parallel::detectCores()-1)

assign_clusters <- function(values) {
  as.integer(as.factor(values)) # Group values into unique clusters
}

parse_simulation_name <- function(name) {
  simulation_details <- sub("res_", "", name)
  simulation_details <- sub(".rds", "", simulation_details)
  details <- strsplit(simulation_details, "_")[[1]]
  list(
    n_clocks = as.integer(details[1]),
    n_events = as.integer(details[2]),
    purity = as.numeric(details[3]),
    coverage = as.integer(details[4]),
    mut_per_seg = as.integer(details[5]),
    seed = as.integer(details[6])
  )
}


my_function <- function (single_inference){
  
  results <- single_inference$df_summary %>%
    inner_join(single_inference$compare_assignment, by = c("segment_id" = "segment_original_indx")) %>%
    group_by(n_components) %>%
    summarise(
      MSE = mse(clock_mean, real_clocks),  # where n_components == unique(time_tickTack)
      RI = rand.index(
        as.numeric(assign_clusters(real_clocks)),
        as.numeric(assign_clusters(clock_mean))
      ),
      # ARI = adj.rand.index(
      #   as.numeric(assign_clusters(real_clocks)),
      #   as.numeric(tickTack_estimate_factor)
      # ),
      .groups = "drop"
    )
  
  return(results)
  
}




calculate_mse <- function(compare_assignment) {
  compare_assignment %>%
    summarise(
      MSE_tickTack = mse(time_tickTack, real_clocks),
      MSE_singleTT = mse(singleTT, real_clocks),
      MSE_MutTimeR = mse(time_MutTime, real_clocks)
    )
}



check_unique_values <- function(compare_assignment) {
  compare_assignment %>%
    summarise(
      n_clocks_real = n_distinct(real_clocks),
      n_clocks_tickTack = n_distinct(time_tickTack),
      difference = n_clocks_real - n_clocks_tickTack
    )
}



find_best_models <- function(model_selection_tibble) {
  model_selection_tibble %>%
    summarise(
      best_BIC = K[which.min(BIC)],
      best_Log_lik = K[which.max(Log_lik)],
      best_ICL = K[which.min(ICL)],
      best_AIC = K[which.min(AIC)],
      best_LOO = K[which.min(LOO)]
    )
}





setwd("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary_new")
files <- list.files()


# Get the information about the simulation for each simulation in the summary_results
variable_names <- files

simulation_info <- lapply(variable_names, parse_simulation_name)
names(simulation_info) <- variable_names

# simulation_info # list with keys as the name of the files with the information about the simulation 






my_function_2 <- function(a){
  single_inference <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary_new/",a))
  simulation_info_single <- simulation_info[[a]]
  
  MS_scores_per_component <- my_function(single_inference)
  
  
  mse_ri_results <- calculate_mse(single_inference$compare_assignment)
  
  check_components_difference <- check_unique_values(single_inference$compare_assignment)
  
  best_models_MS <- find_best_models(single_inference$model_selection_tibble)
  
  returned_item = list(MS_scores_per_component = MS_scores_per_component, mse_ri_results=mse_ri_results, check_components_difference=check_components_difference, best_models_MS = best_models_MS )
}

files <- files
res_single_inference <- parallel::mclapply(files, my_function_2, mc.cores = (parallel::detectCores()-1), mc.allow.recursive = TRUE) #parallel::mclapply  , mc.cores = (parallel::detectCores()-1), mc.allow.recursive = TRUE
names(res_single_inference) <- variable_names
# res_single_inference[[a]]


saveRDS(res_single_inference, "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/data/res_single_inference.rds")


##################################################

aggregated_results <- list()
# indx <- 1
for (indx in seq_along(res_single_inference)) {
  
  sim_info <- simulation_info[[indx]]
  n_clocks <- sim_info$n_clocks
  n_events <- sim_info$n_events
  purity <- sim_info$purity
  coverage <- sim_info$coverage
  seed <- sim_info$seed
  mut_per_seg <- sim_info$mut_per_seg
  
  ms_scores <- res_single_inference[[indx]]$MS_scores_per_component
  mse_ri_results <- res_single_inference[[indx]]$mse_ri_results
  check_components_difference <- res_single_inference[[indx]]$check_components_difference
  best_models_MS <- res_single_inference[[indx]]$best_models_MS
  
  ms_scores_augmented <- ms_scores %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed,
      mut_per_seg = mut_per_seg
    )
  
  mse_ri_augmented <- mse_ri_results %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed,
      mut_per_seg = mut_per_seg
    )
  
  check_diff_augmented <- check_components_difference %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed,
      mut_per_seg = mut_per_seg
    )
  
  best_models_augmented <- best_models_MS %>%
    mutate(
      n_clocks = n_clocks,
      n_events = n_events,
      purity = purity,
      coverage = coverage,
      seed = seed,
      mut_per_seg = mut_per_seg
    )
  
  aggregated_results[[indx]] <- list(
    MS_scores = ms_scores_augmented,
    MSE_RI = mse_ri_augmented,
    Component_Check = check_diff_augmented,
    Best_Models = best_models_augmented
  )
}


final_results <- list(
  MS_scores = bind_rows(lapply(aggregated_results, `[[`, "MS_scores")),
  MSE_RI = bind_rows(lapply(aggregated_results, `[[`, "MSE_RI")),
  Component_Check = bind_rows(lapply(aggregated_results, `[[`, "Component_Check")),
  Best_Models = bind_rows(lapply(aggregated_results, `[[`, "Best_Models"))
)



final_results_combined <- final_results$MS_scores %>%
  left_join(final_results$MSE_RI, by = c("n_clocks", "n_events", "purity", "coverage", "seed","mut_per_seg")) %>%
  left_join(final_results$Component_Check, by = c("n_clocks", "n_events", "purity", "coverage", "seed","mut_per_seg")) %>%
  left_join(final_results$Best_Models, by = c("n_clocks", "n_events", "purity", "coverage", "seed","mut_per_seg"))


saveRDS(final_results_combined, "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/data/final_results_combined.rds")

############################################à



# extract the average results for each parameter 
final_results_avg_per_type_sim <- final_results_combined %>% 
  group_by(n_clocks,n_events, purity, coverage, mut_per_seg) %>% 
  summarise(
    n = n(),
    across(starts_with("MSE"), ~ mean(.x, na.rm = TRUE))
  )


final_results_avg_per_type_sim
mean(final_results_avg_per_type_sim$MSE_tickTack)
# 0.06345808
mean(final_results_avg_per_type_sim$MSE_MutTimeR)
# 0.02549125
mean(final_results_avg_per_type_sim$MSE_singleTT)
# 0.0258387

final_results_purity <- final_results_avg_per_type_sim %>%  
  group_by(n_clocks) %>% 
  summarise(
    n = n(),
    across(starts_with("MSE_"), ~ mean(.x, na.rm = TRUE))
  )

long_data <- final_results_purity %>%
  pivot_longer(cols = starts_with("MSE"), 
               names_to = "Method", 
               values_to = "MSE")

ggplot(long_data, aes(x = n_clocks, y = MSE, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "MSE Comparison Across Methods",
    x = "n_clocks",
    y = "Mean Squared Error (MSE)",
    color = "Method"
  ) +
  scale_color_manual(values = c("blue", "red", "green"))



#############################################################
# MODEL SELECTION 
data <- final_results_combined


unique_simulations <- data %>%
  group_by(purity, coverage, n_events, n_clocks, mut_per_seg, seed) %>%
  slice(1) %>%
  ungroup()

performance <- unique_simulations %>%
  mutate(
    diff_LogLik = abs(best_Log_lik - n_clocks_real),
    diff_BIC = abs(best_BIC - n_clocks_real),
    diff_ICL = abs(best_ICL - n_clocks_real),
    diff_AIC = abs(best_AIC - n_clocks_real),
    diff_LOO = abs(best_LOO - n_clocks_real),
    diff_TickTack = (n_clocks_tickTack - n_clocks_real),
    
    correct_LogLik = (best_Log_lik == n_clocks_real),
    correct_BIC = (best_BIC == n_clocks_real),
    correct_ICL = (best_ICL == n_clocks_real),
    correct_AIC = (best_AIC == n_clocks_real),
    correct_LOO = (best_LOO == n_clocks_real),
    correct_TickTack = (n_clocks_tickTack == n_clocks_real)
    
  )


correct_summary <- performance %>%
  group_by(purity,coverage,n_events,n_clocks,mut_per_seg) %>%
  summarise(
    LogLik = mean(correct_LogLik),
    BIC = mean(correct_BIC),
    ICL = mean(correct_ICL),
    AIC = mean(correct_AIC),
    LOO = mean(correct_LOO),
    TickTack = mean(correct_TickTack)
    
  ) 


correct_summary <- correct_summary %>%
  pivot_longer(cols = c(TickTack), names_to = "Score", values_to = "Proportion") #LogLik, BIC,ICL,AIC,LOO,



ggplot(correct_summary, aes(x = Score, y = Proportion, fill = Score)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ n_clocks, ncol = 5) +
  labs(
    title = "Proportion of Correct Predictions for Each Score",
    x = "Score",
    y = "Proportion of Correct Predictions"
  ) +
  theme_minimal(base_size = 14)





differences <- performance %>%
  group_by(purity,coverage,n_events,n_clocks,mut_per_seg) %>%
  summarise(
    LogLik = mean(diff_LogLik),
    BIC = mean(diff_BIC),
    ICL = mean(diff_ICL),
    AIC = mean(diff_AIC),
    LOO = mean(diff_LOO),
    TickTack = mean(diff_TickTack)
    
  ) 

differences <- differences %>%
  pivot_longer(cols = c(TickTack), names_to = "Score", values_to = "Difference") #(LogLik, BIC,ICL,AIC,LOO,

# Boxplot for differences
ggplot(differences, aes(x = Score, y = Difference, fill = Score)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ n_clocks, ncol = 5) +
  labs(
    title = "Difference Between Chosen and Real Clocks",
    x = "Score",
    y = "Difference"
  ) +
  theme_minimal(base_size = 14)


###########################################################

