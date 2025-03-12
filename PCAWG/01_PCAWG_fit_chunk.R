rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
# .libPaths("~/R/rstudio_v3/")
library(tickTack)
library(parallel)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(seed=123)
source("./plot_cnaqc_conflicts.R")

args <- commandArgs(trailingOnly = TRUE)
chunk_id <- as.integer(args[1])
total_chunks <- as.integer(args[2])
# print(parallel::detectCores()-1)

tolerance = 0.0001
initial_iter = 200
grad_samples = 10
elbo_samples = 100
min_mutations_number = 10
max_distance_CNAqc = 1e6

vector_names <- readRDS("./data/fittable_samples_1Mb_smoothed_10mm.RDS")
num_files <- length(vector_names)
chunk_size <- ceiling(num_files / total_chunks)
start_index <- (chunk_id - 1) * chunk_size + 1
end_index <- min(chunk_id * chunk_size, num_files)
chunk_files <- vector_names[start_index:end_index]

message(sprintf("Processing files %d to %d out of %d", start_index, end_index, num_files))

process_file <- function(s) {
  tryCatch({
    print(s)
    fit <- readRDS(paste0("../../data/clonal_analysis_PCAWG/",s,"/fit.rds"))
    original_dir <- getwd()

    name <- basename(s)   
    dir.create(file.path(original_dir, paste0("results_whole_smoothing/")), showWarnings = TRUE)
    new_dir = paste0(original_dir, paste0("/results_whole_smoothing/"))
    dir.create(file.path(new_dir, paste0("results")), showWarnings = TRUE)
    dir.create(file.path(new_dir, paste0("plots")), showWarnings = TRUE)
    dir.create(file.path(new_dir, paste0("summary")), showWarnings = TRUE)

    fit$mutations <- fit$snvs 
    cnaqc_x <- CNAqc::smooth_segments(fit, maximum_distance = 1e6) 
    
    x <- cnaqc_x
    x$metadata = tibble(purity=cnaqc_x$purity)
  
    x <- tickTack::fit_h(x,
                         max_attempts=2,
                         INIT=TRUE,
                         tolerance = tolerance,
                         initial_iter = initial_iter,
                         grad_samples=grad_samples,
                         elbo_samples=elbo_samples,
                         min_mutations_number = min_mutations_number,
                         n_components = 0) 

    saveRDS(x, paste0(new_dir,"/results/",s,".rds"))
    results <- x$results_timing

    p <- plot_cnaqc(x)
    ggsave(paste0(new_dir,"/plots/plot_h_cnaqc_",s,".png"), plot = p, width = 15, height = 7)
    
    results_model_selection <- tickTack::model_selection_h(results, n_components = 0)
    summarized_results <- results_model_selection$best_fit$summarized_results
    saveRDS(summarized_results, paste0(new_dir,"/summary/smmarized_results_",s,".rds"))

    best_K <- results_model_selection$best_K
  
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", s, e$message))
    
    # Append the failed sample name to a file
    write(s, file = "failed_samples.txt", append = TRUE)
    
    NULL  
  })
  
}

parallel::mclapply(chunk_files, process_file, mc.cores = (detectCores()-1) )
# lapply(chunk_files, process_file)
