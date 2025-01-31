
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

results_path <- "/orfeo/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_whole/"

fits_path = "/orfeo/scratch/cdslab/scocomello/data/clonal_analysis_PCAWG/"

IDs <- list.files(results_path)

RES <- lapply(IDs, function(id) {
  results_model_selection_path = paste0(results_path, id, "/results/results_model_selection.rds")
  if (!file.exists(results_model_selection_path)) {
    print(paste0("Skipping sample id ", id))
    return(NULL)
  }
  
  res_model_selection <- readRDS(results_model_selection_path)
  results <- res_model_selection$best_fit$summarized_results
  
  fit = readRDS(paste0(fits_path, id, "/fit.rds"))
  ttype = strsplit(fit$snvs$project_code, "-")[[1]][1]
  
  df = dplyr::tibble(sample_id=id, ttype=ttype)
  
  dplyr::bind_cols(df, parse_summarized_results(results))
}) %>% do.call("bind_rows", .)

saveRDS(RES, "results/summary_all_samples.rds")


# Add Oncogenes and Tumour suppressor annotation to results
load("data/gene_coordinates_hg19.rda")

TSGs <- c("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC", "CDKN2A", "SMAD4", "VHL", "NF1")
Oncogenes <- c("MYC", "KRAS", "BRAF", "EGFR", "HER2", "ALK", "PIK3CA", "ABL1", "CCND1", "NRAS")
genes_of_interest <- c(TSGs, Oncogenes)

gene_coordinates_hg19 <- gene_coordinates_hg19 %>% 
  dplyr::filter(gene %in% genes_of_interest)


library(purrr)

RESfiltered = RES[map_lgl(seq_len(nrow(RES)), function(i) {
  any(gene_coordinates_hg19$chr == RES$chr[i] & gene_coordinates_hg19$from >= RES$from[i] & gene_coordinates_hg19$to <= RES$to[i])
}), ]

genes_annotation = lapply(1:nrow(RESfiltered), function(i) {
  r = RESfiltered[i,]
  gene = gene_coordinates_hg19 %>% 
    dplyr::filter(from >= r$from, chr==r$chr, to <= r$to) %>% 
    dplyr::pull(gene)
  
  if (length(gene) > 1) {
    dplyr::tibble(gene = paste(c(gene), collapse = "-"), type = "Multiple")
  } else if (gene %in% Oncogenes) {
    dplyr::tibble(gene = gene, type = "Oncogene")
  } else {
    dplyr::tibble(gene = gene, type = "TumourSuppressor")
  }
}) %>% do.call("bind_rows", .)

RES_w_genes = bind_cols(RESfiltered, genes_annotation)
saveRDS(RES_w_genes, "results/res_w_onco_and_ts.rds")
