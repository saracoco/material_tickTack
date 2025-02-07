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

x_after_inference_sim <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/res_1_5_0.3_20_10_1.rds")
x_after_inference_PCAWG <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole/2a8d63eb-0174-4213-9214-413f391f512c/results/x_after_inference.rds")
fit <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/data/clonal_analysis_PCAWG/2a8d63eb-0174-4213-9214-413f391f512c/fit.rds")

fit$mutations <- fit$snvs

x_after_inference_sim$reference_genome <- "GRCh38"
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
plot_segments_tick_tack(x)


merge_timing_and_segments(x, K=2)
ggsave("plot_sample_2.png")


# plot_1 <- plot_segments_tick_tack_data(x)






