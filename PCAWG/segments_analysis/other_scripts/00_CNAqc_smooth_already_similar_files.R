rm(list=ls())
.libPaths(new="~/R/rstudio_v3/") 
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("../../utils.R")

set.seed(seed=123)
data_path <-  "../../../../data/clonal_analysis_PCAWG"
sample_files = list.files(data_path, full.names = "T")

idx=3
sample_path <- paste0(sample_files[idx], "/fit.rds")
fit = readRDS(sample_path)

fit$mutations <- fit$snvs 
cnaqc_x <- CNAqc::init(mutations = fit$snvs, cna = fit$cna,purity = fit$purity)
cnaqc_x <- CNAqc::smooth_segments(cnaqc_x)

# CNAqc::plot_smoothing(cnaqc_x)
# CNAqc::plot_multisample_CNA(list(`Before` = cnaqc_x$before_smoothing, `After` = cnaqc_x))

cnaqc_x$snvs <- cnaqc_x$mutations
is_fittable(cnaqc_x, min_mutations_number = 10)  


# x = list( mutations = cnaqc_x$snvs, cna = cnaqc_x$cna, metadata= tibble(purity=cnaqc_x$purity))
# x <- tickTack::fit_h(x,
#                      max_attempts=2,
#                      INIT=TRUE,
#                      tolerance = tolerance,
#                      initial_iter = initial_iter,
#                      grad_samples=grad_samples,
#                      elbo_samples=elbo_samples,
#                      min_mutations_number = min_mutations_number)
# 
# 
# 
# results <- x$results_timing
# saveRDS(x, paste0(new_dir,"/results/x_after_inference.rds"))
# 
# results_model_selection <- tickTack::model_selection_h(results, n_components = 0)
# saveRDS(results_model_selection, paste0(new_dir,"/results/results_model_selection.rds"))
# 
# summarized_results <- results_model_selection$best_fit$summarized_results
# saveRDS(summarized_results, paste0("results/",name,".rds"))
# 
# best_K <- results_model_selection$best_K
# model_selection_tibble <- results_model_selection$model_selection_tibble
# entropy <- results_model_selection$entropy_list
# 
# K = nrow(results_model_selection$model_selection_tibble)
# 
# p_elbo <- list()
# for (i in 1:K){
#   p_elbo[[i]] <- tickTack::plot_elbo_h(results$elbo_iterations[[as.character(i)]]) + ggplot2::ggtitle(paste0("K = ", i))
# }
# p_elbo <- gridExtra::grid.arrange(grobs = p_elbo, ncol = 2)  #add global title
# ggplot2::ggsave(paste0(new_dir,"/plots/plot_elbo.png"),plot = p_elbo, width = 30, height = 30)
# 
# 
# 
# p <- tickTack::plot_timing_h(results, best_K)
# ggsave(paste0(new_dir,"/plots/plot_timing_h.png"),plot = p, width = 25, height = 5)
# 
# print(results$draws_and_summary[[best_K]]$summarized_results, n = nrow(results$draws_and_summary[[best_K]]$summarized_results))
# print(results_model_selection$model_selection_tibble)
# 
# 
# plot_model_selection_inference <- list()
# for (i in 1:K){
#   plot_model_selection_inference[[i]] <- tickTack::plot_timing_h(results, i) + ggplot2::ggtitle(paste0("K = ", i))
# }
# plot_model_selection_inference <- gridExtra::grid.arrange(grobs = plot_model_selection_inference, nrow = K) 
# ggsave(paste0(new_dir,"/plots/plot_timing_all_K_h.png"),plot = plot_model_selection_inference, width = 25, height = 30)





# fittable_samples_4_min_mutations <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/data/fittable_samples_4_min_mutations.RDS")
# vector_names <- list.files("../../data/clonal_analysis_PCAWG", full.names = TRUE)
# 
# vector_names_1 <- list.files("../../data/clonal_analysis_PCAWG", full.names = TRUE)
# vector_names_2 <- fittable_samples_4_min_mutations
# 
# ids_1 <- basename(vector_names_1)
# ids_2 <- basename(vector_names_2)
# 
# common_ids <- intersect(ids_1, ids_2)
# 
# filtered_vector_names <- vector_names_1[basename(vector_names_1) %in% common_ids]
# saveRDS(fittable_paths, "data/fittable_samples.RDS")



# 
# sample_files = list.files(data_path, full.names = "T")
# 
# 
# sample_segments_info <- lapply(1:length(sample_files), function(idx) {
#   print(paste0("Completion = ", idx / length(sample_files) * 100, "%"))
#   sample_path <- paste0(sample_files[idx], "/fit.rds")
#   fit = readRDS(sample_path)
#   avg_mut_per_seg(fit)
# })
# 
# saveRDS(sample_segments_info, "data/sample_segments_info.RDS")
# 
# avg_mut_tot <- lapply(1:length(sample_segments_info), function(idx) {
#   s <- sample_segments_info[idx][[1]]
#   if(all(s==FALSE)){
#     NA
#   }else{
#     avg_mut <- sample_segments_info[idx][[1]]$avg_mut_x_seg
#     avg_mut
#   }
# }) %>% unlist()
# 
# 
# median <- median(avg_mut_tot[avg_mut_tot<50], na.rm=TRUE)
# 
# ggplot(data = data.frame(x = avg_mut_tot), aes(y = x), na.rm=TRUE) +  
#   geom_boxplot(fill = "lightblue", outlier.color = "red", width = 0.3) +
#   labs(title = "Boxplot of Values", y = "Values", x = "") +
#   ylim(c(0,50))+
#   theme_minimal()
# 

