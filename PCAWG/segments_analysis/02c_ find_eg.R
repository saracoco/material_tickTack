
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")
library(RColorBrewer)


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
  best_K <- res_model_selection$best_K
  
  results <- res_model_selection$best_fit$summarized_results
  
  fit = readRDS(paste0(fits_path, id, "/fit.rds"))
  ttype = unique(fit$snvs$ttype)[2]
  
  df = dplyr::tibble(sample_id=id, ttype=ttype, best_K=best_K, purity=fit$purity)
  
  dplyr::bind_cols(df, parse_summarized_results(results))
}) %>% do.call("bind_rows", .)

saveRDS(RES, "../../select_egsample_purity_.rds")

RES <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/select_egsample_purity_.rds")

filtered_res <- RES %>%
  filter(best_K == 1, purity>0.8 ) %>%  
  group_by(sample_id) %>% 
  summarize(n = n(), .groups = "drop")
  ungroup() %>%  
  select(sample_id, best_K, ttype, purity) 

  
  
print(filtered_res, n=200)
saveRDS(filtered_res, "../../select_egsample_many_clusters.rds")




RES_HM <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/data/RES.rds")



WGD <- RES_HM %>% filter(is_HM=="WGD", n_clusters==1)
WGD[6,]$sample_id
HM <- RES_HM %>% filter(is_HM=="HM", n_clusters==2, wgd_status=="no_wgd") #  
HM[4,]$sample_id


HM <- RES_HM %>% filter(is_HM=="Classic", n_clusters==4, wgd_status=="no_wgd") #  


mytable <- table(filtered_res$ttype)
threshold <- 20

tumor_data <- data.frame(
  Category = names(mytable),
  SampleSize = as.numeric(mytable),
  Group = ifelse(mytable < threshold, "Other", names(mytable))
)

tumor_data_grouped <- aggregate(SampleSize ~ Group, data = tumor_data, sum)
tumor_data_grouped$Category <- tumor_data_grouped$Group
tumor_data_grouped <- tumor_data_grouped[, c("Category", "SampleSize")]
tumor_data_sorted <- tumor_data_grouped[order(-tumor_data_grouped$SampleSize), ]
lbls <- paste(tumor_data_sorted$Category, tumor_data_sorted$SampleSize)
colors <- brewer.pal(min(nrow(tumor_data_sorted), 12), "Set3")
pie(tumor_data_sorted$SampleSize, labels = lbls, col = colors, main = "Pie Chart of Tumour Type\n (with Sample Sizes)")



filtered_res <- RES %>%
  group_by(sample_id) %>% 
  slice(1) %>%  
  ungroup() %>%  
  select(sample_id, best_K,ttype)  



ggplot(filtered_res, aes(x = factor(best_K))) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(title = "Distribution of best_K",
       x = "best_K",
       y = "Count") +
  theme_minimal()


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
