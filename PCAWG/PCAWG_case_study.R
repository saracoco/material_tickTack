# tarfile <- "clonal_analysis_PCAWG.tar.gz"
# data <- read.delim(file = untar(tarfile,compressed="gzip"),sep="\t")
.libPaths(new="~/R/rstudio_v3/") 
library(tickTack)
library(dplyr)
library(ggplot2)
library(parallel)
  
print(parallel::detectCores()-1)

set.seed(seed=123)
spath <-  "../../data"
sfile <-  "clonal_analysis_PCAWG.tar.gz"

outputdir <- "../../data"
data <- untar(file.path(spath,sfile), exdir = outputdir)

# Rename columns from list
# setNames(old = c(cna), 
#          new = c(segments))
# fit <- fit %>% mutate(segments = fit$cna)

library(tibble)

tolerance = 0.0001
iter = 200
grad_samples = 10
elbo_samples = 100
# samples_metadata <- readRDS("./samples_info.rds")
# load data

vector_names <- list.files("../../data/clonal_analysis_PCAWG/")


# s <- vector_names[1]


parallel::mclapply(vector_names, function(s){
  
  tryCatch({
  
  fit <- readRDS(paste0("../../data/clonal_analysis_PCAWG/",s,"/fit.rds"))
  original_dir <- getwd()
  
  print(s)
  
  
  dir.create(file.path(original_dir, paste0("results_tickTack/",s)), showWarnings = FALSE)
  new_dir = paste0(original_dir, paste0("/results_tickTack/",s))
  dir.create(file.path(new_dir, paste0("results")), showWarnings = FALSE)
  dir.create(file.path(new_dir, paste0("plots")), showWarnings = FALSE)
  
  
  x <- list( mutations = tibble(chr = fit$snvs$chr, 
                                from = fit$snvs$from, 
                                to = fit$snvs$to, 
                                ref = fit$snvs$ref, 
                                alt = fit$snvs$alt, 
                                DP = fit$snvs$DP, 
                                NV = fit$snvs$NV, 
                                VAF = fit$snvs$NV/fit$snvs$DP, 
                                sample = 1), 
             cna = tibble(chr = fit$cna$chr, 
                          from = fit$cna$from, 
                          to = fit$cna$to, 
                          Major = fit$cna$Major, 
                          minor = fit$cna$minor,   
                          CCF = 0, 
                          total_cn = Major + minor), 
             metadata = tibble(purity = fit$purity) 
  )
  
  
  data <- x
  
  x <- tickTack::fit_h(x, 
                       max_attempts=2, 
                       INIT=TRUE, 
                       tolerance = tolerance,
                       initial_iter = iter,
                       grad_samples=grad_samples,
                       elbo_samples=elbo_samples)
  
  
  
  results <- x$results_timing
  saveRDS(x, paste0(new_dir,"/results/x_after_inference.rds"))
  
  results_model_selection <- tickTack::model_selection_h(results, n_components = 0)
  saveRDS(results_model_selection, paste0(new_dir,"/results/results_model_selection.rds"))
  
  best_K <- results_model_selection$best_K
  model_selection_tibble <- results_model_selection$model_selection_tibble
  entropy <- results_model_selection$entropy_list
  
  # posterior_clocks <- tickTack::plot_posterior_clocks_h(results, 2)
  # posterior_weights <- tickTack::plot_posterior_weights_h(results, 2)
  
  K = nrow(results_model_selection$model_selection_tibble)
  
  p_elbo <- list()
  for (i in 1:K){
    p_elbo[[i]] <- tickTack::plot_elbo_h(results$elbo_iterations[[i]]) + ggplot2::ggtitle(paste0("K = ", i))
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
  plot_model_selection_inference <- gridExtra::grid.arrange(grobs = plot_model_selection_inference, nrow = K) #add global title
  ggsave(paste0(new_dir,"/plots/plot_timing_all_K_h.png"),plot = plot_model_selection_inference, width = 25, height = 30)
  
  
  
  
  # Single segment inference # Single segment inference # Single segment inference 
  segments <- data$cna
  mutations <- data$mutations
  purity <- data$metadata$purity
  
  # Run the fit function
  results_single <- fit(
    segments = segments,
    mutations = mutations,
    purity = purity,
    possible_k = c("2:1", "2:2", "2:0"),
    beta_binomial = TRUE
  )
  
  saveRDS(results_single, paste0(new_dir, "/results/results_single.rds"))
  
  
  print(results_single$summarized_results, n= nrow(results_single$summarized_results) )
  
  p <- tickTack::plot_timing(results_single, segments, colour_by = "karyotype")
  ggsave(paste0(new_dir,"/plots/plot_timing.png"),plot = p, , width = 25, height = 5)
  
  
  # mtationtimer inference 
  
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", s, e$message))
    NULL  # Return NULL if there is an error
  })
  
  
}, mc.cores = (parallel::detectCores()-1), mc.allow.recursive = TRUE)









