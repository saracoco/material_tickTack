rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("utils.R")

set.seed(seed=123)
data_path <-  "data/clonal_analysis_PCAWG/"

sample_files = list.files(data_path, full.names = "T")

fittable_flags <- lapply(1:length(sample_files), function(idx) {
  print(paste0("Completion = ", idx / length(sample_files) * 100, "%"))
  sample_path <- paste0(sample_files[idx], "/fit.rds")
  fit = readRDS(sample_path)
  is_fittable(fit)  
}) %>% unlist()

fittable_flags %>% sum()
fittable_paths <- sample_files[fittable_flags]
saveRDS(fittable_paths, "data/fittable_samples.RDS")

