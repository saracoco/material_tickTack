# tarfile <- "clonal_analysis_PCAWG.tar.gz"
# data <- read.delim(file = untar(tarfile,compressed="gzip"),sep="\t")
.libPaths(new="~/R/rstudio_v3/") 
library(tickTack)

library(dplyr)
library(tibble)

vector_names <- list.files("/orfeo/cephfs/scratch/cdslab/scocomello/Simulations_tickTack/clonal_analysis_PCAWG/")

samples <- read.csv("./TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep="\t")
sample_metadata <- dplyr::tibble()
for (f in vector_names) {
  s <- gsub(".rds", "", f)
  t <-samples %>%
    filter(sample_id == s) %>%
    select(ttype) %>%
    unlist() %>%
    unique()
  sample_metadata <- dplyr::bind_rows(sample_metadata, dplyr::tibble(sample_id = s, ttype = t, file_name = f))
}
rm(f,s,t, samples)

# Get unique tumor types
sample_metadata %>% pull(ttype) %>% unique() %>% length()
# saveRDS(sample_metadata, "./clonal_analysis_PCAWG_info.rds")

vector_names_check = tibble(sample_id = vector_names)
info_sample <- right_join(sample_metadata,vector_names_check)
# table(info_sample$ttype)
unique(info_sample$ttype)
sum(is.na(info_sample$ttype))






sample_metadata

unique(info_sample$ttype)
PCAWG_statistics <- lapply(sample_metadata$file_name, function(sam){
    
    fit <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/Simulations_tickTack/clonal_analysis_PCAWG/",sam,"/fit.rds"))
    
     return (dplyr::tibble(
    ttype = sample_metadata$ttype[which(sample_metadata$file_name == sam)],
    n_cna = fit$n_cna,
    n_snvs = fit$n_snvs,
    purity = fit$purity,
    most_prevalent_karyotype = fit$most_prevalent_karyotype,
    l_karyotype = fit$l_karyotype,
    n_karyotype_2_0 = fit$n_karyotype["2:0"][[1]],
    n_karyotype_2_1 = fit$n_karyotype["2:1"][[1]],
    n_karyotype_2_2 = fit$n_karyotype["2:2"][[1]],
    length = mean(fit$cna$length) ))
  
}) %>%bind_rows()

saveRDS(PCAWG_statistics, "./PCAWG_statistics_ICGC.rds")









summary_statistics <- lapply(unique(info_sample$ttype), function(t){
    
   s <- PCAWG_statistics %>% filter(ttype == t)
    
   dplyr::tibble(n_cna = mean(s$n_cna, na.rm = T), 
                  n_cna_sd = sd(s$n_cna, na.rm = T),
                  
                  n_snvs = mean(s$n_snvs, na.rm = T), 
                  n_snvs_sd = sd(s$n_snvs, na.rm = T), 
                  
                  purity = mean(s$purity, na.rm = T), 
                  purity_sd = sd(s$purity, na.rm = T), 
                  
                  n_karyotype_2_0 = (mean(s$n_karyotype_2_0, na.rm = T)/n_snvs)*n_cna,
                  n_karyotype_2_1 = (mean(s$n_karyotype_2_1, na.rm = T)/n_snvs)*n_cna,
                  n_karyotype_2_2 = (mean(s$n_karyotype_2_2, na.rm = T)/n_snvs)*n_cna,
                  
                  length = mean(s$length, na.rm = T),
                  length_sd = sd(s$length, na.rm = T),
                  
                  n_karyotype_2_0_frq = n_karyotype_2_0/(n_karyotype_2_0 + n_karyotype_2_1 + n_karyotype_2_2),
                  n_karyotype_2_1_frq = n_karyotype_2_1/(n_karyotype_2_0 + n_karyotype_2_1 + n_karyotype_2_2),
                  n_karyotype_2_2_frq = n_karyotype_2_2/(n_karyotype_2_0 + n_karyotype_2_1 + n_karyotype_2_2),
                  
                  most_common_karyotype = max(n_karyotype_2_0, n_karyotype_2_1, n_karyotype_2_2),
                  ttype = unique(s$ttype),
                  n_samples = nrow(s))
  
}) %>% bind_rows()

summary <-summary_statistics %>% bind_rows()

summary_statistics <- subset(summary_statistics, select = -c(n_karyotype_2_0,n_karyotype_2_1,n_karyotype_2_2) )

saveRDS(summary_statistics, "./PCAWG_summary_statistics_ICGC.rds")





