
.libPaths(new="~/R/rstudio_v3/") 
library(tickTack)
library(dplyr)
library(ggplot2)
library(parallel)



get_files_success = function(a, counter) {
  n_files <- length(list.files(paste0("./",a,"/results")))
  if (n_files == 0){
    return(a)
  }
  else { return(NULL)}
}

results_time_parallel <- list.files("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_tickTack_parallel_2")
successes_files_parallel <-(lapply(results_time_parallel, get_files_success)%>%unlist())
