
.libPaths(new="~/R/rstudio/") 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("utils.R")

parse_summarized_results <- function(results) {
  results %>% 
    tidyr::separate(segment_name, into = c("chr", "from", "to"), sep = "_", convert = TRUE) %>% 
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>% 
    dplyr::select(chr, from, to, segment_id, karyotype, clock_mean, clock_low, clock_high, clock_rank)
}


ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  select(sample_id, ttype) %>% 
  dplyr::distinct()

results_path <- "/orfeo/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_tickTack_parallel_2/"

IDs <- list.files(results_path)


id = "0624eb1d-3aff-4037-a3c5-fc363a9edd02"
RES <- lapply(IDs, function(id) {
  results_model_selection_path = paste0(results_path, id, "/results/results_model_selection.rds")
  if (!file.exists(results_model_selection_path)) {
    print(paste0("Skipping sample id ", id))
    return(NULL)
  }
  res_model_selection <- readRDS(results_model_selection_path)
  results <- res_model_selection$best_fit$summarized_results
  
  df <- ttypes %>% 
    dplyr::filter(sample_id == id)
  
  if (nrow(df) == 0) {
    df = dplyr::tibble(sample_id=id, ttype="Uknown")
  }
  
  
  dplyr::bind_cols(df, parse_summarized_results(results))
}) %>% do.call("bind_rows", .)

saveRDS(RES, "results/summary_all_samples.rds")

# Look at statistics
RES %>% 
  dplyr::select(ttype, sample_id) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(n) %>% 
  dplyr::mutate(ttype = factor(ttype, levels=ttype)) %>% 
  ggplot2::ggplot(mapping = aes(x=ttype, y=n)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  theme_bw()

RES %>% 
  dplyr::group_by(ttype, karyotype) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::mutate(f = n / sum(n)) %>% 
  ggplot2::ggplot(mapping = aes(x=karyotype, y=n)) +
  ggplot2::geom_col() +
  facet_wrap(.~ttype, scales = "free_y") +
  theme_bw()


