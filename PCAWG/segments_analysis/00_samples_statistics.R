rm(list=ls())
.libPaths(new="~/R/rstudio_v3/") 
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("../utils.R")

data_path <-  "../../../data/clonal_analysis_PCAWG"

sample_files = list.files(data_path, full.names = "T")


library(future.apply)

# Set up parallel execution
options(future.globals.maxSize = 2*1024**3)

plan(multisession, workers = parallel::detectCores() - 1)  # Use all but one core



for (i in c(5e5, 1e6, 1e7)) {
  info_segments <- future_lapply(1:(length(sample_files)), function(idx) {
    
    print(paste0("Completion = ", idx / length(sample_files) * 100, "%"))
    sample_path <- paste0(sample_files[idx], "/fit.rds")
    
    fit <- tryCatch(readRDS(sample_path), error = function(e) return(NULL))
    if (is.null(fit)) return(NULL)  # Skip if reading fails
    
    fit$mutations <- fit$snvs 
    cnaqc_x <- CNAqc::init(mutations = fit$snvs, cna = fit$cna, purity = fit$purity)
    cnaqc_x <- CNAqc::smooth_segments(cnaqc_x, maximum_distance = i)  # Assuming 'i' is the distance
    
    karyotypes <- c("2:0", "2:1", "2:2")
    cna <- tryCatch(cnaqc_x$cna, error = function(e) return(NULL))
    
    if (!is.null(cna)) {
      tryCatch({
        cna$k <- paste(cna$Major, cna$minor, sep = ':')
        cna <- filter(cna, k %in% karyotypes)
        cna_simple <- cna
      }, error = function(e) return(NULL))
    }
    
    if (is.null(cna_simple)) return(NULL)  # Skip if no valid CNA data
    
    
    n_mutations_per_segment <- lapply(1:(length(cna_simple)), function(segment_idx) {
      
      
      segment <- cna_simple[segment_idx, ]
      chr <- segment$chr                                            
      segment_id <- paste(chr, segment$from, segment$to, sep = "_")
      k <- segment$k
      
      segment_mutations <- cnaqc_x$mutations %>%
        dplyr::filter(chr == segment$chr,.data$from > segment$from, .data$to < segment$to) %>%
        tidyr::drop_na(DP)
      
      n_mutations = tryCatch(
        nrow(segment_mutations), 
        error = function(e) as.character(NA)
      )
      
      return(n_mutations)
    })%>% unlist()
    
    
    
    info_segments_single <- list(
      sample = unique(cnaqc_x$mutations$sample),
      n_cna_tot = tryCatch(cnaqc_x$n_cna_clonal, error = function(e) NA),
      n_cna_simple = tryCatch(nrow(cna_simple), error = function(e) NA),
      type = tryCatch(
        strsplit(cnaqc_x$mutations$project_code, "-")[[1]][1], 
        error = function(e) as.character(NA)
      ),
      median_length = tryCatch(median(cnaqc_x$cna %>% filter(segment_id %in% cna_simple$segment_id) %>% pull(length)), error = function(e) NA),
      sd_length = tryCatch(sd(cnaqc_x$cna %>% filter(segment_id %in% cna_simple$segment_id) %>% pull(length)), error = function(e) NA),
      meadian_mutation_number = tryCatch(median(n_mutations_per_segment), error = function(e) NA),
      sd_mutation_number = tryCatch(sd(n_mutations_per_segment), error = function(e) NA),
      min_mutation_number = tryCatch(min(n_mutations_per_segment), error = function(e) NA),
      max_mutation_number = tryCatch(max(n_mutations_per_segment), error = function(e) NA)
      
    )
    
    return(info_segments_single)
    
  }) %>% bind_rows()
  
  saveRDS(info_segments, paste0("./data/info_segments_", i, ".rds"))
}

# Clean up parallel workers
plan(sequential)
