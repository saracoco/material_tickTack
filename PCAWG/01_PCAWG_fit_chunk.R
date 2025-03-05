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
min_mutations_number = 10

vector_names <- readRDS("segments_analysis/data/failed_inference.rds")
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

    print(s)


    name <- basename(s)   
    dir.create(file.path(original_dir, paste0("results_whole_check_fails/",name)), showWarnings = TRUE)
    new_dir = paste0(original_dir, paste0("/results_whole_check_fails/",name))
    dir.create(file.path(new_dir, paste0("results")), showWarnings = TRUE)
    dir.create(file.path(new_dir, paste0("plots")), showWarnings = TRUE)



    x = list( mutations = fit$snvs, cna = fit$cna, metadata= tibble(purity=fit$purity))

    data <- x

    x <- tickTack::fit_h(x,
                         max_attempts=2,
                         INIT=TRUE,
                         tolerance = tolerance,
                         initial_iter = initial_iter,
                         grad_samples=grad_samples,
                         elbo_samples=elbo_samples,
                         min_mutations_number = min_mutations_number)



    results <- x$results_timing
    saveRDS(x, paste0(new_dir,"/results/x_after_inference.rds"))

    results_model_selection <- tickTack::model_selection_h(results, n_components = 0)
    saveRDS(results_model_selection, paste0(new_dir,"/results/results_model_selection.rds"))

    summarized_results <- results_model_selection$best_fit$summarized_results
    saveRDS(summarized_results, paste0("results_check_fails/",name,".rds"))

    best_K <- results_model_selection$best_K
    model_selection_tibble <- results_model_selection$model_selection_tibble
    entropy <- results_model_selection$entropy_list

    K = nrow(results_model_selection$model_selection_tibble)

    p_elbo <- list()
    for (i in 1:K){
      p_elbo[[i]] <- tickTack::plot_elbo_h(results$elbo_iterations[[as.character(i)]]) + ggplot2::ggtitle(paste0("K = ", i))
    }
    p_elbo <- gridExtra::grid.arrange(grobs = p_elbo, ncol = 2)  #add global title
    ggplot2::ggsave(paste0(new_dir,"/plots/plot_elbo.png"),plot = p_elbo, width = 30, height = 30)



    p <- tickTack::plot_timing_h(results, best_K)
    ggsave(paste0(new_dir,"/plots/plot_timing_h.png"),plot = p, width = 25, height = 5)

    print(results$draws_and_summary[[best_K]]$summarized_results, n = nrow(results$draws_and_summary[[best_K]]$summarized_results))
    print(results_model_selection$model_selection_tibble)


    plot_model_selection_inference <- list()
    for (i in 1:K){
      plot_model_selection_inference[[i]] <- tickTack::plot_timing_h(results, i) + ggplot2::ggtitle(paste0("K = ", i))
    }
    plot_model_selection_inference <- gridExtra::grid.arrange(grobs = plot_model_selection_inference, nrow = K) 
    ggsave(paste0(new_dir,"/plots/plot_timing_all_K_h.png"),plot = plot_model_selection_inference, width = 25, height = 30)



    # Single segment inference 
    segments <- data$cna
    mutations <- data$mutations
    purity <- data$metadata$purity

    results_single <- fit(
      segments = segments,
      mutations = mutations,
      purity = purity,
      possible_k = c("2:1", "2:2", "2:0"),
      min_mutations_number = min_mutations_number,
      beta_binomial = TRUE
    )

    saveRDS(results_single, paste0(new_dir, "/results/results_single.rds"))

    result_tibble_single <- results_single$summarized_results
    saveRDS(result_tibble_single, paste0("results_check_fails/",name,"_result_tibble_single.rds"))

    p <- tickTack::plot_timing(results_single, segments, colour_by = "karyotype")
    ggsave(paste0(new_dir,"/plots/plot_timing.png"),plot = p, , width = 25, height = 5)



  }, error = function(e) {
    message(sprintf("Error processing %s: %s", s, e$message))
    NULL  
  })
  
}

parallel::mclapply(chunk_files, process_file, mc.cores = (detectCores()-1) )
lapply(chunk_files[1:3], process_file)
