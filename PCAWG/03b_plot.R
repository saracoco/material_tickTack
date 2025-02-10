.libPaths(new="~/R/rstudio_v3/") 
library(CNAqc)
library(tickTack)
library(dplyr)
library(stringr)
library(fuzzyjoin)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyr)
source("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/simulations_v2/utils_plot.R")


eg_1_component <- "0a9c9db0-c623-11e3-bf01-24c6515278c0"
eg_2_component <- "2a8d63eb-0174-4213-9214-413f391f512c"
eg_3_component <- "05616329-e7ba-4efd-87b1-d79cd0f7af3d"
name <- eg_3_component

# x_after_inference_sim <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/res_1_5_0.3_20_10_1.rds")
x_after_inference_PCAWG <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole_40mut_AIC/", name, "/results/x_after_inference.rds"))
fit <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/data/clonal_analysis_PCAWG/", name, "/fit.rds"))

fit$mutations <- fit$snvs

# x_after_inference_sim$reference_genome <- "GRCh38"
x_after_inference_PCAWG$reference_genome <- fit$reference_genome
x_after_inference_PCAWG$cna_subclonal <- fit$cna_subclonal
x_after_inference_PCAWG$peaks_analysis <- fit$peaks_analysis
x_after_inference_PCAWG$ploidy <- fit$ploidy
x_after_inference_PCAWG$purity <- fit$purity
x_after_inference_PCAWG$n_mutations <- fit$n_snvs
x_after_inference_PCAWG$n_cna_clonal <- fit$n_cna_clonal
x_after_inference_PCAWG$ttype <- unique(fit$snvs$ttype)[2] 


x <- x_after_inference_PCAWG
results<-x$results_timing

K=3
# plot_segments_tick_tack(x)


merge_timing_and_segments(x, K=K)
ggsave("../../plot_sample_3_CNA.pdf",height = 8,width = 15)


# plot_1 <- plot_segments_tick_tack_data(x)






