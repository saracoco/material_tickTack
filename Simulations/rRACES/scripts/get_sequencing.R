get_sequencing <- function(phylo_forest, coverage = 80, purity=0.9, chromosomes = c("5","6","10","12","13","22") ){
  
  # phylo_forest <- load_phylogenetic_forest("results/phylo_forest.sff")
  
  cov <- coverage
  # s <- simulate_seq(phylo_forest, coverage = cov, chromosomes = "1", write_SAM = FALSE)
  chromosomes <- chromosomes
  
  seq_results <- parallel::mclapply(chromosomes, function(c) {
    simulate_seq(phylo_forest, coverage = cov, purity = purity, chromosomes = c, write_SAM = FALSE)
  }, mc.cores = parallel::detectCores()) 
  
  
  seq_results_final <- data.frame()
  for(i in 1: length(chromosomes)){
    seq_results_final <- bind_rows(seq_results_final,seq_results[[i]]$mutations)
  }
  
  # %>% do.call("bind_rows", .)
  
  # saveRDS(seq_results, file = paste0('data/seq_', cov, 'X.RDS'))
  
    return(list(seq_results = seq_results, seq_results_final = seq_results_final))
}
