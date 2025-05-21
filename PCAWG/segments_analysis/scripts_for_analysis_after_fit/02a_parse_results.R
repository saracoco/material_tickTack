rm(list = ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
library(tickTack)
source("utils.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Please provide PCAWG_INFO_DIR and TICK_TACK_FITS_DIR and OUTPUT_DIR as command line arguments.")
}

PCAWG_INFO_DIR = args[1]
TICK_TACK_FITS_DIR <- args[2]
OUTPUT_DIR <- args[3]

cat("PCAWG info directory:", PCAWG_INFO_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")

if (!dir.exists(PCAWG_INFO_DIR)) {
  stop(paste0("Input dir: ", PCAWG_INFO_DIR," does not exist"))
}

if (!dir.exists(TICK_TACK_FITS_DIR)) {
  stop(paste0("Input dir: ", TICK_TACK_FITS_DIR," does not exist"))
}

if (!dir.exists(OUTPUT_DIR)) {
  message("Creating output dir")
  dir.create(OUTPUT_DIR, recursive = T)
}

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

#results_path <- "/results/"
#fits_path = "/data/clonal_analysis_PCAWG/"

IDs <- list.files(TICK_TACK_FITS_DIR)
IDs = IDs[!grepl("single", IDs)]

id = IDs[1]
RES <- lapply(IDs, function(id) {
  print(which(IDs == id) / length(IDs) * 100)
  
  tryCatch({
    #results = readRDS(paste0(results_path, id))
    tickTack_res = readRDS(file.path(TICK_TACK_FITS_DIR, id))
    pcawg_data = readRDS(paste0(PCAWG_INFO_DIR, "/", unlist(strsplit(id, ".rds")), "/fit.rds"))
    tumour_name = (unique(pcawg_data$snvs$ttype) %>% na.omit())[1]
    ploidy = pcawg_data$ploidy
    ttype = strsplit(pcawg_data$snvs$project_code, "-")[[1]][1]
    df = dplyr::tibble(sample_id=id, ttype=ttype, ploidy=ploidy, tumour_name=tumour_name)
    return(dplyr::bind_cols(df, parse_summarized_results(tickTack_res)))
  }, error = function(e) {
    # Error handling
    print(paste0("An error occurred:", id))
    return(NULL)
  })
  
}) %>% do.call("bind_rows", .)


saveRDS(RES, file.path(OUTPUT_DIR, "summary_all_samples.rds"))
RES = readRDS(file.path(OUTPUT_DIR, "summary_all_samples.rds"))

RES %>% 
  dplyr::select(ttype, tumour_name) %>% 
  dplyr::distinct()

# Add arm-level annotation

# Add Oncogenes and Tumour suppressor annotation to results
load("data/gene_coordinates_hg19.rda")

TSGs <- c("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC", "CDKN2A", "SMAD4", "VHL", "NF1")
Oncogenes <- c("MYC", "KRAS", "BRAF", "EGFR", "HER2", "ALK", "PIK3CA", "ABL1", "CCND1", "NRAS")
DNA_repair = c("RAD51")

# TSGs <- c("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC", "CDKN2A", "SMAD4", "VHL", "NF1",
#          "ARID1A", "ATM", "BAP1", "CDH1", "FANCD2", "MLH1", "NBN", "PBRM1", "RAD51C", "STK11",
#          "MEN1", "TSC1", "TSC2", "DICER1", "CHEK2", "WT1", "EP300", "KMT2D", "KMT2C", "FBXW7",
#          "SETD2", "SMARCA4", "NF2", "CDKN1B", "RUNX1", "SUFU", "SMARCB1")

# Oncogenes <- c("MYC", "KRAS", "BRAF", "EGFR", "HER2", "ALK", "PIK3CA", "ABL1", "CCND1", "NRAS",
#                "GNAS", "IDH1", "JAK2", "KIT", "MDM2", "MET", "NTRK1", "PDGFRA", "RET", "TERT",
#                "FGFR1", "FGFR2", "FGFR3", "HRAS", "CDK4", "CDK6", "NOTCH1", "BCL2", "BCL6", 
#                "FOXA1", "FOXP1", "ETV1", "ETV4", "ETV5", "FOS", "JUN", "SRC", "RSPO2", "RSPO3")

genes_of_interest <- c(TSGs, Oncogenes, DNA_repair)

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
  } else if (gene %in% DNA_repair) {
    dplyr::tibble(gene = gene, type = "DNA_repair")
  }else {
    dplyr::tibble(gene = gene, type = "TumourSuppressor")
  }
}) %>% do.call("bind_rows", .)

RES_w_genes = bind_cols(RESfiltered, genes_annotation)
saveRDS(RES_w_genes, file.path(OUTPUT_DIR, "res_w_onco_and_ts.rds"))
