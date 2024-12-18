.libPaths(new="~/R/rstudio_v3/") 
library(MutationTimeR)
library(dplyr)
library(tickTack)
library(tibble)
source("../scripts/parse_gerstrung_modified.R")
original_directory = getwd()

source("../scripts/get_simulation_MutationTimeR.R")



simulations_gerstrung_data = function (n_clocks=3, 
                                       n_events=8, 
                                       purity=0.9, 
                                       coverage=100, 
                                       epsilon=0.20, 
                                       seed = 123, 
                                       tolerance = 0.01,
                                       max_attempts = 2,
                                       INIT = TRUE,
                                       min_mutations_number = 3) {
  
  print("simulate gerstrung data")

  print (n_clocks) 
  print (n_events)
  print(purity) 
  print(coverage)
  print(epsilon)
  print(seed) 
  print(tolerance) 
  print(max_attempts) 

  
  
  res_simulate <- get_simulation_MutationTimeR(number_clocks=n_clocks, 
                                               number_events=n_events, 
                                               purity=purity, 
                                               coverage=coverage, 
                                               epsilon=epsilon, 
                                               seed = seed)
  df <- res_simulate$df
  vcf <- res_simulate$vcf
  cn <- res_simulate$cn
  cn_timed <- res_simulate$cn_timed
  clusters <- res_simulate$clusters
  purity <- res_simulate$purity
  ####################  MUTATIONTIMER INFERENCE   #######################################################
  mt <- mutationTime(vcf, cn_timed, clusters=clusters, n.boot=10)
  # saveRDS(mt,"../MutationTimeR/mt_42_3_clocks.rds")
  mcols(cn_timed) <- cbind(mcols(cn_timed),mt$T)
  
  info(header(vcf)) <- rbind(info(header(vcf)),MutationTimeR:::mtHeader())
  info(vcf) <- cbind(info(vcf), mt$V)
  
  plot_MutTime <- plotSample(vcf,cn_timed)
  plotSample(vcf,cn_timed)
  ggsave("plots/plot_Muttime.png", height=5, width=10)
  
  res_MutTime <- list(vcf = vcf, cn_timed = cn_timed)
  
  
  ############## convert data for tickTack from vcf format #############
  
  x_gerstrung <- parse_gerstrung_modified(vcf = vcf, cna = cn)
  # saveRDS(data_gerstrung,"../MutationTimeR/data_gerstrung_42_3_clocks.rds")
  # x_gerstrung$metadata = tibble(purity=purity)
  
  x_gerstrung <- list( mutations = tibble(chr = x_gerstrung$mutations$chr, from = x_gerstrung$mutations$from, to = x_gerstrung$mutations$to, ref = "", alt = "", DP = x_gerstrung$mutations$DP, NV = x_gerstrung$mutations$NV, VAF = x_gerstrung$mutations$NV/x_gerstrung$mutations$DP, sample = 1), 
                       cna = tibble(chr = x_gerstrung$cna$chr, from = x_gerstrung$cna$from, to = x_gerstrung$cna$to , Major = x_gerstrung$cna$Major, minor = x_gerstrung$cna$minor,   CCF = 0, total_cn = Major + minor), 
                       metadata = tibble(purity = purity) 
  )
  
  ##################### TICKTACK INFERENCE ############
  x_gerstrung <- tickTack::fit_h(x_gerstrung, max_attempts=max_attempts, INIT=INIT, tolerance = tolerance)
  
  
  # save results 
  results_simulated <- x_gerstrung$results_timing

  results_model_selection <-  tickTack::model_selection_h(results_simulated)

  best_K <- results_model_selection$best_K
  model_selection_tibble <- results_model_selection$model_selection_tibble
  
  results_tickTack <- list( results = results_simulated, results_model_selection = results_model_selection )
  cat("Best K =",best_K)
  
  k_number <- nrow(results_model_selection$model_selection_tibble)
  
  
  library(ggplot2)
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
  
  compare_assignment

  
  cn_timed$chr <- cn_timed$time
  
   df_cn_timed <- as.data.frame(cn_timed$chr)
   colnames(df_cn_timed)[1] <- "time_MutTime" 
   df_cn_timed$chr <- rep(1:23)
  compare_assignment <- compare_assignment %>%
    left_join(df_cn_timed, by = "chr") 
  
  compare_assignment <- compare_assignment %>% dplyr::mutate(clock_real_factor = as.numeric(factor(taus))) %>%
    dplyr::mutate(tickTack_estimate_factor = as.numeric(factor(clock_mean)))%>%
    dplyr::mutate(MutTime_estimate_factor = as.numeric(factor(time_MutTime)))%>%
    dplyr::rename(real_clocks = taus)%>%
    dplyr::rename(time_tickTack = clock_mean)
  
  ################################ SINGLE SEGMENT #######################################################################
  
  library(tidyr)
  
  fit <- tickTack::fit(segments = x_gerstrung$cna, mutations = x_gerstrung$mutations, purity = purity,
                       alpha = .05, min_mutations_number = 3,
                       beta_binomial = F, beta_binomial_disp = 0.01)

  
  plot_single <- tickTack::plot_timing(fit, x_gerstrung$cna)

  
  compare_assignment <- compare_assignment %>%
    left_join(fit$summarized_results, by = "chr") %>% 
    dplyr::rename(singleTT = tau_mean)
  
  
  res = list(data = df,
             compare_assignment = compare_assignment,
             plot_Muttime = plot_MutTime, res_MutTime = res_MutTime, 
             plot_tickTack = plot_tickTack, res_tickTack = results_tickTack, 
             plot_SingleTT = plot_single , res_SingleTT = fit )
  
  return(res)
  
}






