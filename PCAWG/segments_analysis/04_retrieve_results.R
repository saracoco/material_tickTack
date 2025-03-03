.libPaths(new="~/R/rstudio_v3/") 
library(tickTack)
library(parallel)
library(dplyr)
library(ggplot2)
set.seed(seed=123)

args <- commandArgs(trailingOnly = TRUE)
chunk_id <- as.integer(args[1])
total_chunks <- as.integer(args[2])

print(parallel::detectCores()-1)

tolerance = 0.0001
initial_iter = 200
grad_samples = 10
elbo_samples = 100
min_mutations_number = 40

vector_names <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/data/fittable_paths_min_mut_10.rds")
num_files <- length(vector_names)
print(num_files)

chunk_size <- ceiling(num_files / total_chunks)
start_index <- (chunk_id - 1) * chunk_size + 1
end_index <- min(chunk_id * chunk_size, num_files)

chunk_files <- vector_names[start_index:end_index]

message(sprintf("Processing files %d to %d out of %d", start_index, end_index, num_files))

process_file <- function(s) {
  tryCatch({
    print(s)
    fit <- readRDS(paste0("",s,"/fit.rds"))
    original_dir <- getwd()
    
    print(s)
    
    
    name <- basename(s)   
    new_dir = paste0(original_dir, paste0("/results_whole/",name))

    x = list( mutations = fit$snvs, cna = fit$cna, metadata= tibble(purity=fit$purity))
    data <- x
    
    x <- readRDS(paste0(new_dir,"/results/x_after_inference.rds"))
    results <- x$results_timing

    results_model_selection <-     x <- readRDS(paste0(new_dir,"/results/results_model_selection.rds"))

    summarized_results <- results_model_selection$best_fit$summarized_results
    saveRDS(summarized_results, paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_ZENODO/",name,".rds"))
    
    best_K <- results_model_selection$best_K
    model_selection_tibble <- results_model_selection$model_selection_tibble

    K = nrow(results_model_selection$model_selection_tibble)
  
    print(results$draws_and_summary[[best_K]]$summarized_results, n = nrow(results$draws_and_summary[[best_K]]$summarized_results))
    print(results_model_selection$model_selection_tibble)
    

    
    results_single <- readRDS(paste0(new_dir, "/results/results_single.rds"))
    
    
    result_tibble_single <- results_single$summarized_results
    saveRDS(result_tibble_single, paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_ZENODO/",name,"_result_tibble_single.rds"))
    
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", s, e$message))
    NULL  # Return NULL if there is an error
  })
  
}

parallel::mclapply(chunk_files, process_file, mc.cores = (parallel::detectCores()-1) )
