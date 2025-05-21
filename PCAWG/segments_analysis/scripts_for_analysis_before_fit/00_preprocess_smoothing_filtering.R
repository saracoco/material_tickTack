# tarfile <- "clonal_analysis_PCAWG.tar.gz"
# data <- read.delim(file = untar(tarfile,compressed="gzip"),sep="\t")
.libPaths(new="~/R/rstudio_v3/") 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
# source("avg_mut_per_seg.R")
source("../../utils.R")

set.seed(seed=123)

max_distance_CNAqc = 5e6
data_path <-  "../../../../data/clonal_analysis_PCAWG"
sample_files = list.files(data_path, full.names = "T")

fittable_flags <- lapply(1:length(sample_files), function(idx) {
  print(paste0("Completion = ", idx / length(sample_files) * 100, "%"))
  sample_path <- paste0(sample_files[idx], "/fit.rds")
  fit = readRDS(sample_path)
  fit$mutations <- fit$snvs 
  cnaqc_x <- CNAqc::init(mutations = fit$snvs, cna = fit$cna,purity = fit$purity)
  cnaqc_x <- CNAqc::smooth_segments(cnaqc_x, maximum_distance = max_distance_CNAqc)
  
  cnaqc_x$snvs <- cnaqc_x$mutations
  is_fittable(cnaqc_x, min_mutations_number = 15)  
}) %>% unlist()

# fittable_flags %>% sum()
fittable_paths <- sample_files[fittable_flags]
saveRDS(fittable_paths, "data/fittable_paths.rds")

