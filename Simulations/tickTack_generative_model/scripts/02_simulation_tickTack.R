source("../scripts/convert_to_vcf.R")
source("../scripts/convert_to_granges.R")
source("../scripts/03_get_simulation_tickTack.R")
library(ggplot2)

simulation_tickTack = function (n_clocks=3, 
                                       n_events=8, 
                                       purity=0.9, 
                                       coverage=100, 
                                       epsilon=0.20, 
                                       seed = 123, 
                                       tolerance = 0.001,
                                       max_attempts = 2,
                                       INIT = TRUE,
                                       min_mutations_number = 4) {
  
  print("02_simulation_tickTack")

  
  res_simulate <- get_simulation_tickTack(number_clocks=n_clocks, 
                                               number_events=n_events, 
                                               purity=purity, 
                                               coverage=coverage, 
                                               epsilon=epsilon, 
                                               seed = seed)
  
  data_simulation = as.data.frame(res_simulate$data_simulation) %>% mutate(chr= 1:length(res_simulate$data_simulation$taus))
  x = res_simulate$x
  
  df = x$cna %>% left_join(data_simulation)
  saveRDS(x, "input_data.rds")
  
 # timing inference ticktack hierarchical
  x <- tickTack::fit_h(x, max_attempts=max_attempts, INIT=INIT, tolerance = tolerance)

  results_simulated <- x$results_timing
  results_model_selection <- tickTack::model_selection_h(results_simulated)
  best_K <- results_model_selection$best_K
  results_tickTack <- list( results = results_simulated, results_model_selection = results_model_selection )
  
  cat("Best K =",best_K)
  
  
  k_number <- nrow(results_model_selection$model_selection_tibble)
  
  p <- list()
  for (i in 1:k_number){
    p[[i]] <- tickTack::plot_timing_h(results_simulated, i, split_contiguous_segments = FALSE) + ggplot2::ggtitle(paste0("K = ", i))
    if (i == best_K){
      print(i)
      png(filename="./plots/tickTack.png", height=250, width=700)
      tickTack::plot_timing_h(results_simulated, i, split_contiguous_segments = FALSE) + ggplot2::ggtitle(paste0("K = ", i,"best K: ",best_K))
      dev.off()
    }else{NULL}
  }
  png(filename="./plots/tickTack_multiple.png", height=1000, width=700)
  plot_tickTack <- gridExtra::grid.arrange(grobs = p, nrow=k_number)  #add global title
  dev.off()
  # plot_tickTack <- p[[best_K]]
  
  
  
  # compare assignment for all the 3 methods so that it is easy to handle the results for analysis

  clock_assignment <- results_simulated$draws_and_summary[[best_K]]$summarized_results %>% dplyr::select(segment_original_indx, segment_id, karyotype, chr, clock_mean)
  
  clock_assignment <- clock_assignment %>% mutate(chr = as.integer(segment_original_indx))
  
  compare_assignment <- clock_assignment %>% 
    left_join(df, by = "chr") 
  
  
  compare_assignment %>% dplyr::rename(time_tickTack = clock_mean)


  
  ############## convert data from tickTack to vcf format #############
  
  library(VariantAnnotation)
  library(stringr)
  
  vcf <- convert_to_vcf(x$mutations, x$cna, x$metadata)
  bb <- convert_to_granges(x$cna, x$mutations, x$metadata)
  clusters <- data.frame(cluster=1:1, proportion=c(purity), n_ssms=c(100))

  mt <- MutationTimeR::mutationTime(vcf, bb, clusters=clusters, n.boot=10)
  mcols(bb) <- cbind(mcols(bb),mt$T)
  plot_MutTime <- MutationTimeR::plotSample(vcf,bb)
  MutationTimeR::plotSample(vcf,bb)
  ggsave("plots/plot_Muttime.png", height=5, width=10)
  
  res_MutTime <- list(vcf = vcf, cn_timed = bb)
  cn_timed = bb
  
  cn_timed$chr <- cn_timed$time
  
  df_cn_timed <- as.data.frame(cn_timed$chr)
  colnames(df_cn_timed)[1] <- "time_MutTime" 
  df_cn_timed$chr <- rep(1:length(res_simulate$data_simulation$taus))
  compare_assignment <- compare_assignment %>%
    left_join(df_cn_timed, by = "chr") 
  
  #################################
  
  compare_assignment <- compare_assignment %>% dplyr::mutate(clock_real_factor = as.numeric(factor(taus))) %>%
    dplyr::mutate(tickTack_estimate_factor = as.numeric(factor(clock_mean)))%>%
    dplyr::mutate(MutTime_estimate_factor = as.numeric(factor(time_MutTime)))%>%
    dplyr::rename(real_clocks = taus)%>%
    dplyr::rename(time_tickTack = clock_mean)


  
  
  ################################ SINGLE SEGMENT #######################################################################
  
  library(tidyr)
  
  fit <- tickTack::fit(segments = x$cna, mutations = x$mutations, purity = purity,
                       alpha = .05, min_mutations_number = 3,
                       beta_binomial = F, beta_binomial_disp = 0.01)

  
  plot_single <- tickTack::plot_timing(fit, x$cna)


  compare_assignment <- compare_assignment %>%
    left_join(fit$summarized_results, by = "chr") %>% 
    dplyr::rename(singleTT = tau_mean)
  

  res = list(data = df,
             compare_assignment = compare_assignment,
             plot_MutTime = plot_MutTime, res_MutTime = res_MutTime, 
             plot_tickTack = plot_tickTack, res_tickTack = results_tickTack, 
             plot_SingleTT = plot_single , res_SingleTT = fit )
  
  return(res)
  
}






