# tarfile <- "clonal_analysis_PCAWG.tar.gz"
# data <- read.delim(file = untar(tarfile,compressed="gzip"),sep="\t")

.libPaths(new="~/R/rstudio_v3/") 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("avg_mut_per_seg.R")

set.seed(seed=123)
data_path <-  "../../data/clonal_analysis_PCAWG"

sample_files = list.files(data_path, full.names = "T")


sample_segments_info <- lapply(1:length(sample_files), function(idx) {
  print(paste0("Completion = ", idx / length(sample_files) * 100, "%"))
  sample_path <- paste0(sample_files[idx], "/fit.rds")
  fit = readRDS(sample_path)
  avg_mut_per_seg(fit)
})

saveRDS(sample_segments_info, "data/sample_segments_info.RDS")



avg_mut_tot <- lapply(1:length(sample_segments_info), function(idx) {
  s <- sample_segments_info[idx][[1]]
  if(all(s==FALSE)){
    NA
  }else{
    avg_mut <- sample_segments_info[idx][[1]]$avg_mut_x_seg
    avg_mut
    }
  
}) %>% unlist()


median <- median(avg_mut_tot[avg_mut_tot<50], na.rm=TRUE)



ggplot(data = data.frame(x = avg_mut_tot), aes(y = x), na.rm=TRUE) +  
  geom_boxplot(fill = "lightblue", outlier.color = "red", width = 0.3) +
  labs(title = "Boxplot of Values", y = "Values", x = "") +
  ylim(c(0,50))+
  theme_minimal()



data_path <-  "../../data/clonal_analysis_PCAWG/"

sample_files = list.files(data_path, full.names = "T")

fittable_flags <- lapply(1:length(sample_files), function(idx) {
  print(paste0("Completion = ", idx / length(sample_files) * 100, "%"))
  sample_path <- paste0(sample_files[idx], "/fit.rds")
  fit = readRDS(sample_path)
  is_fittable(fit, min_mutations_number=10)  
}) %>% unlist()

fittable_flags %>% sum()
fittable_paths <- sample_files[fittable_flags]

saveRDS(fittable_paths, "data/fittable_paths_min_mut_10.rds")

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