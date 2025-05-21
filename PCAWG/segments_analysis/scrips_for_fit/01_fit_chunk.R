rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
# .libPaths("~/R/rstudio_v3/")
library(tickTack)
library(parallel)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(seed=123)
# source("./plot_cnaqc_conflicts.R")

args <- commandArgs(trailingOnly = TRUE)
chunk_id <- as.integer(args[1])
total_chunks <- as.integer(args[2])
# print(parallel::detectCores()-1)

tolerance = 0.0001
initial_iter = 200
grad_samples = 10
elbo_samples = 100
min_mutations_number = 15
max_distance_CNAqc = 1e6

vector_names <- readRDS("./data/fittable_paths.rds")
num_files <- length(vector_names)
chunk_size <- ceiling(num_files / total_chunks)
start_index <- (chunk_id - 1) * chunk_size + 1
end_index <- min(chunk_id * chunk_size, num_files)
chunk_files <- vector_names[start_index:end_index]


message(sprintf("Processing files %d to %d out of %d", start_index, end_index, num_files))

process_file <- function(s) {
  tryCatch({
    print(s)

    fit <- readRDS(paste0(s,"/fit.rds"))
    original_dir <- getwd()

    name <- basename(s)   
    dir.create(file.path(original_dir, paste0("results")), showWarnings = TRUE)
    new_dir = paste0(original_dir, paste0("/results"))
    dir.create(file.path(new_dir, paste0("results_extended")), showWarnings = TRUE)
    dir.create(file.path(new_dir, paste0("plots")), showWarnings = TRUE)
    dir.create(file.path(new_dir, paste0("results_summary")), showWarnings = TRUE)
 
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

    saveRDS(x, paste0(new_dir,"/results_extended/",s=basename(s),".rds"))
    results <- x$results_timing
    # p <- plot_cnaqc(x)
    # ggsave(paste0(new_dir,"/plots/plot_h_cnaqc_",s,".png"), plot = p, width = 15, height = 7)
  
    results_model_selection <- tickTack::model_selection_h(results, n_components = 0)
    best_K <- results_model_selection$best_K

    summarized_results <- results_model_selection$best_fit$summarized_results
    saveRDS(summarized_results,paste0(new_dir,"/results_summary/",basename(s),".rds"))
    
    # best_K <- results_model_selection$best_K
    # model_selection_tibble <- results_model_selection$model_selection_tibble
    # entropy <- results_model_selection$entropy_list
    # posterior_clocks <- tickTack::plot_posterior_clocks_h(results, 2)
    # posterior_weights <- tickTack::plot_posterior_weights_h(results, 2)
    
    K = nrow(results_model_selection$model_selection_tibble)
    # 
    # p_elbo <- list()
    # for (i in 1:K){
    #   p_elbo[[i]] <- tickTack::plot_elbo_h(results$elbo_iterations[[as.character(i)]]) + ggplot2::ggtitle(paste0("K = ", i))
    # }
    # p_elbo <- gridExtra::grid.arrange(grobs = p_elbo, ncol = 2)  #add global title
    # ggplot2::ggsave(paste0(new_dir,"/plots/plot_elbo.png"),plot = p_elbo, width = 30, height = 30)
    
    
    p <- tickTack::plot_timing_h(results, best_K)
    ggsave(paste0(new_dir,"/plots/plot_timing_h",basename(s),".png"),plot = p, width = 25, height = 5)
    
    # print(results$draws_and_summary[[best_K]]$summarized_results, n = nrow(results$draws_and_summary[[best_K]]$summarized_results))
    print(results_model_selection$model_selection_tibble)
    
    plot_model_selection_inference <- list()
    for (i in 1:K){
      plot_model_selection_inference[[i]] <- tickTack::plot_timing_h(results, i) + ggplot2::ggtitle(paste0("K = ", i))
    }
    plot_model_selection_inference <- gridExtra::grid.arrange(grobs = plot_model_selection_inference, nrow = K) #add global title
    ggsave(paste0(new_dir,"/plots/plot_timing_all_K_h",basename(s),".png"),plot = plot_model_selection_inference, width = 25, height = 30)
    
    
    
    
    # Single segment inference # Single segment inference # Single segment inference 
    segments <- x$cna
    mutations <- x$mutations
    purity <- x$metadata$purity
    
    # Run the fit function
    results_single <- fit(
      segments = segments,
      mutations = mutations,
      purity = purity,
      possible_k = c("2:1", "2:2", "2:0"),
      min_mutations_number = 15,
      beta_binomial = FALSE
    )
    
    saveRDS(results_single, paste0(new_dir, "/results_extended/",basename(s),"_results_single.rds"))
    
    result_tibble_single <- results_single$summarized_results 
    saveRDS(result_tibble_single, paste0(new_dir, "/results_summary/",basename(s),"_result_tibble_single.rds"))
    
    
    p <- tickTack::plot_timing(results_single, segments, colour_by = "karyotype")
    ggsave(paste0(new_dir,"/plots/",basename(s),"plot_timing_s.png"),plot = p, width = 25, height = 5)
    
  
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", s=basename(s), e$message))
    
    # Append the failed sample name to a file
    write(basename(s), file = "failed_samples.txt", append = TRUE)
    
    NULL  
  })
  
}

parallel::mclapply(chunk_files, process_file, mc.cores = (detectCores()-1) )
# lapply(chunk_files, process_file)
