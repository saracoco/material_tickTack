.libPaths(new="~/R/rstudio_v3/") 
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("../utils.R")

data_path <-  "../../../data/clonal_analysis_PCAWG/"

sample_path <- paste0(data_path,"01658141-8398-4585-9f0f-8355dd9b0604", "/fit.rds")

fit <- tryCatch(readRDS(sample_path), error = function(e) return(NULL))

tolerance = 0.0001
initial_iter = 200
grad_samples = 10
elbo_samples = 100
min_mutations_number = 10



x = list( mutations = fit$snvs, cna = fit$cna, metadata= tibble(purity=fit$purity))

x_fin = fit
x_fin$mutations = fit$snvs
x_fin$metadata = tibble(purity=fit$purity)

x <- tickTack::fit_h(x_fin,
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
saveRDS(summarized_results, paste0("results/",name,".rds"))

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




source("../../tickTack/R/plot_cnaqc.R")

merge_timing_and_segments(x)
