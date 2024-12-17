# source("utils.R")

# source("../simulate_tissue.R")
# source("../place_mutations.R")
# source("../get_sequencing.R")
# source("../plot_races.R")
# source("../results_races.R")


simulate_rRACES <-  function(purity = 0.90, coverage = 80, chromosomes = c("5","6","10","12","13","22"), seed = 123, K_clocks = 3, n_events = 8, simulation_time = 500, epsilon = 0.20, type_length = c("long")){
 
  # seed = 123
  # coverage = 0.1
  # purity = 0.3

  set.seed(seed)
  
  library(rRACES)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  #####
  # K_clocks = 3
  # n_events = 8
  # simulation_time = 500
  # epsilon = 0.20
  # chromosomes = c("5","6","10","12","13", "22")
  # type_length = c("long")

  # forest <- load_samples_forest("samples_forest.sff")
  # phylo_forest <- load_phylogenetic_forest("results/phylo_forest.sff")
  # seq_results <- readRDS(paste0("results/seq_final_", coverage, "X.RDS"))
  # seq_results_old <- readRDS(paste0("results/seq_", coverage, "X.RDS"))
  #
  #
  # source("place_mutations.R")
  #####
  
  vector_karyo <- c("2:0", "2:1", "2:2")
  weights_karyo <- c(0.33, 0.33, 0.33)
  
  # vector_tau = rep(0, K_clocks)
  vector_tau = runif(K_clocks, 0, 1)
  while(min(diff(vector_tau)) <= epsilon) {
    vector_tau = runif(K_clocks, 0, 1)
  }
  
  
  # for (j in 1:K_clocks){
  #   vector_tau[j] = rbeta(1,2, 2)
  #   if (j != 1){
  #     while (!all ( abs(vector_tau[1:j-1] - vector_tau[j]) > epsilon  )   ){
  #       vector_tau[j] = rbeta(1,2, 2)
  #     }
  #   }
  # }
  weights_tau <- rep(1/K_clocks, K_clocks)
  
  
  get_taus_karyo = function (number_events,
                                vector_tau, vector_karyo,
                                weigths_tau, weights_karyo ){

    taus <- sample( vector_tau[1:length(vector_tau)], number_events, replace=TRUE, prob=weigths_tau )
    karyo <- sample( vector_karyo[1:length(vector_karyo)], number_events, replace=TRUE, prob=weights_karyo )
  
    return(list(taus=taus,karyo=karyo))
  
  }
  
  
  
  get_races_tau_conversion = function (taus_pseudo,
                                       simulation_time){
    
    taus_RACES <- rep(0, length(taus_pseudo))

    for(i in 1:length(taus_pseudo)){
      taus_RACES[i] <- taus_pseudo[i]*simulation_time
    }
  
    return(taus_RACES)
    
  }
  
  data_simulation <- get_taus_karyo(n_events, vector_tau, vector_karyo, weights_tau, weights_karyo)
  saveRDS(data_simulation, file = paste0('results/data_simulation', coverage, 'X.RDS'))
  
  
  
  taus_pseudo <- data_simulation$taus
  taus <- get_races_tau_conversion(taus_pseudo, simulation_time)
  karyo = data_simulation$karyo
  
  ############# SIMULATE TISSUE ##################################
   results_simulate <- simulate_tissue(K_clocks = K_clocks, n_events = n_events, taus = sort(unique(taus)), simulation_time = simulation_time)
  
   forest <- results_simulate$samples_forest
   forest$save("results/samples_forest.sff")
   
   pl <- results_simulate$pl
   ggsave("plots/pl_tissue.pdf", plot = pl, height = 20, width = 20)

   ############ PLACE MUTATIONS #################################
   # long_or_short_vec = sample(c("long", "short"), n_events, replace = TRUE)
   
   long_or_short_vec = sample(type_length, n_events, replace = TRUE)
   results_mutations <- place_mutations(forest = forest, karyo = karyo, chromosomes = chromosomes, taus = taus, long_or_short_vec = long_or_short_vec)

   phylo_forest <- results_mutations$phylo_forest
   phylo_forest$save("results/phylo_forest.sff")

   sticks <- results_mutations$sticks
   annot_forest <- results_mutations$annot_forest
   exp_timeline <- results_mutations$exp_timeline
   pl <- results_mutations$pl


   ggplot2::ggsave('plots/sticks.png', plot = sticks, width = 210, height = 297, units = "mm", dpi=300)

   ggplot2::ggsave('plots/annot_forest.png', plot = annot_forest, width = 210, height = 297, units = "mm", dpi=300)
   ggplot2::ggsave('plots/exp_timeline.png', plot = exp_timeline, width = 210, height = 297, units = "mm", dpi=300)

   ggplot2::ggsave('plots/mutations.png', plot = pl, width = 210, height = 297, units = "mm", dpi=300)
   
   # all_SNV = results_mutations$all_SNV

   # ############  GET SEQUENCING ##################################

   result_seq_results <- get_sequencing(phylo_forest = phylo_forest, coverage = coverage, purity = purity, chromosomes = chromosomes)
   
   # saveRDS(result_seq_results$seq_results, file = paste0('results/seq_', coverage, 'X.RDS'))
   # saveRDS(result_seq_results$seq_results_final, file = paste0('results/seq_final_', coverage, 'X.RDS'))
   
   saveRDS(result_seq_results, file = paste0('results/result_seq_results_', coverage, 'X.RDS'))
   

   # ############ PLOT #############################################
   results_plot <- plot_races(result_seq_results)

   ggplot2::ggsave(filename = 'plots/seq_plot_baf_dr.png', dpi = 300, plot = results_plot$gw_plot,  width = 410, height = 297, units = "mm")
   ggplot2::ggsave(filename = 'plots/seq_plot_vaf.png', dpi = 300, plot = results_plot$vaf,  width = 410, height = 297, units = "mm")

   ggplot2::ggsave('plots/hist_results_plot.png', plot = results_plot$hist, width = 210, height = 297, units = "mm", dpi=300)
   ggplot2::ggsave('plots/p_results_plot.pdf', plot = results_plot$p, width = 297, height = 50, units = "mm", dpi=300)


   # ############ RESULTS ##########################################
   # results <- results_races(phylo_forest)
   # saveRDS(results, file = paste0('data/results.RDS'))
   # 
   
}

