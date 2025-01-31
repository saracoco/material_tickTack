
rm(list = ls())
require(tidyverse)

# Read data ####
main_path = "/orfeo/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results/"
sim_list <- list.files(main_path)

NREP = 5

whole_res <- lapply(sim_list, function(s) {
  info = unlist(strsplit(s, "_"))
  if (length(info) == 7) {
    print(s)
    
    n_clocks <- as.numeric(info[3])
    n_segments <- as.numeric(info[4])
    purity <- as.numeric(info[5])
    coverage = as.numeric(info[6])
    n_muts = as.numeric(info[7])
    
    rr <- lapply(1:NREP, function(i) {
      results_files = list.files(paste0(main_path, s, "/", i, "/results"), full.names = T)
      results_files = results_files[grepl("compare", results_files)]
      
      if (length(results_files) == 1) {
        readRDS(results_files) %>% 
          dplyr::select(segment_original_indx, time_tickTack, time_MutTime, singleTT, real_clocks) %>% 
          dplyr::mutate(iter = i, n_clocks=n_clocks, n_segments=n_segments, purity=purity, coverage=coverage, n_muts=n_muts)    
      }
    }) %>% do.call("bind_rows", .)
    
    return(rr)
  }  
}) %>% do.call("bind_rows", .)

saveRDS(whole_res, "processed_results/simulations_results.rds")