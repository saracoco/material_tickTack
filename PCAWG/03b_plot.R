# .libPaths(new="~/R/rstudio_v3/")

rm(list = ls())
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


# 1 cluster 
# 0554ffe5-31f7-43f5-8372-2b73c9cf3527 no
# 09497b9b-6fca-48cb-af97-161a3e434a51 no
# 209a9b10-7129-48fe-a899-d14ba17efe6f si HM ma ploidy 3 
# 8c5f9574-c622-11e3-bf01-24c6515278c0 si HM ma picco vaf nero 
# ee368301-b2f7-428f-9330-67c054c1a09d si HM ma picco vaf nero
# fc9ef456-75a2-5967-e040-11ac0c484477 no
# 96a2896c-1e32-4827-a526-6b7104832f9a no
# 
# 
09cb8bc5-13ac-44ac-9b7d-6de143373570 si WGD
05616329-e7ba-4efd-87b1-d79cd0f7af3d classic 3 comp
0bfd1068-3fdf-a95b-e050-11ac0c4860c3 classsic 4 comp

# 2 cluster no wgd 
# 28e81540-4744-4865-b627-c7c9d8a3c2b8 no
# 28839c75-90a8-493f-b658-8c63e0ebd324 si bassa pi 0.4
# 28e81540-4744-4865-b627-c7c9d8a3c2b8 no






eg <- "2a8d63eb-0174-4213-9214-413f391f512c"



name <- eg

# x_after_inference_sim <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/res_1_5_0.3_20_10_1.rds")
x_after_inference_PCAWG <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole/", name, "/results/x_after_inference.rds"))
res_model_selection <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole/", name, "/results/results_model_selection.rds"))
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
x_after_inference_PCAWG$ttype <- fit$snvs$project_code
x_after_inference_PCAWG$K = res_model_selection$best_K


x <- x_after_inference_PCAWG
results<-x$results_timing

p2 <- merge_timing_and_segments(x)


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












