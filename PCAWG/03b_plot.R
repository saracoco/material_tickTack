library(CNAqc)
library(tickTack)
library(dplyr)
library(stringr)
library(fuzzyjoin)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyverse)
library(mobster)
#source("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/simulations_v2/utils_plot.R")
source("../simulations_v2/utils_plot.R")


eg_OV <- "2a8d63eb-0174-4213-9214-413f391f512c"
eg_PRAD <- "0bfd1043-8180-e3e4-e050-11ac0c4860c5"
eg <- eg_OV



name <- eg

# x_after_inference_sim <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/res_1_5_0.3_20_10_1.rds")
x_after_inference_PCAWG <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole/", name, "/results/x_after_inference.rds"))
res_model_selection <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole/", name, "/results/results_model_selection.rds"))
fit <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/data/clonal_analysis_PCAWG/", name, "/fit.rds"))
x_after_inference_PCAWG_40 <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole_40mut_AIC/", name, "/results/x_after_inference.rds"))
res_model_selection_40 <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole_40mut_AIC/", name, "/results/results_model_selection.rds"))

fit$mutations <- fit$snvs

OV <- CNAqc::plot_segments(fit_OV)
PRAD <- CNAqc::plot_segments(fit_PRAD)

ggplot2::ggsave("../../plot_paper_tickTack/OV_segments.pdf",height = 8,width = 20, plot=OV)




x_after_inference_PCAWG <- x_after_inference_PCAWG_40
# x_after_inference_sim$reference_genome <- "GRCh38"
x_after_inference_PCAWG$reference_genome <- fit$reference_genome
x_after_inference_PCAWG$cna_subclonal <- fit$cna_subclonal
x_after_inference_PCAWG$peaks_analysis <- fit$peaks_analysis
x_after_inference_PCAWG$ploidy <- fit$ploidy
x_after_inference_PCAWG$purity <- fit$purity
x_after_inference_PCAWG$n_mutations <- fit$n_snvs
x_after_inference_PCAWG$n_cna_clonal <- fit$n_cna_clonal
x_after_inference_PCAWG$ttype <- fit$snvs$project_code
x_after_inference_PCAWG$K = res_model_selection$best_K


x <- x_after_inference_PCAWG
results<-x$results_timing

OV <- merge_timing_and_segments(x)
ggplot2::ggsave("../../plot_paper_tickTack/OV_40.pdf",height = 8,width = 20, plot=OV)


p1 <- merge_timing_and_segments(x)

p3 <- merge_timing_and_segments(x)
p <- p3 / p2 / p1
ggsave("plot/PCAWG_examples_supplementary_WGD.png", units = "in", dpi = 600, width = 16, height = 8, plot = p3)
ggsave("plot/PCAWG_examples_supplementary.png", units = "in", dpi = 900, width = 12, height = 14, plot = p)

ggsave("../../plot_paper_tickTack/1_component.pdf",height = 8,width = 20)




sum_results <- x$results_timing$draws_and_summary[[as.character(x_after_inference_PCAWG$K)]]
table(sum_results$summarized_results$clock_mean)

tab <- table(sum_results$summarized_results$clock_mean, sum_results$summarized_results$karyotype)
tab
# plot_1 <- plot_segments_tick_tack_data(x)












