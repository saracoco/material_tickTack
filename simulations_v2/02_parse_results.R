
rm(list = ls())
require(tidyverse)

res_folder = "results/" 
res_names = list.files(res_folder)

r = res_names[1]

all_res = lapply(res_names, function(r) {
  print(r)
  
  info = unlist(strsplit(r, "_"))
  n_clocks = as.numeric(info[2])
  n_events = as.numeric(info[3])
  purity = as.numeric(info[4])
  coverage = as.numeric(info[5])
  n_mutations = as.numeric(info[6 ])
  
  sub_dirs = list.files(paste0(res_folder, r))
  sub_i = sub_dirs[1]
  
  sub_res = lapply(sub_dirs, function(sub_i) {
    fp = paste0(res_folder, r, "/", sub_i)
    
    if (file.exists(paste0(fp, "/merged_res.rds"))) {
      readRDS(paste0(fp, "/merged_res.rds")) %>% 
        dplyr::mutate(n_clocks=n_clocks, n_events=n_events, purity=purity, coverage=coverage, n_mutations=n_mutations, i.iter = as.numeric(sub_i))    
    }
    
  }) %>% do.call("bind_rows", .)
  
  sub_res  
}) %>% do.call("bind_rows", .)

saveRDS(all_res, "results_summarised/whole_res.RDS")


# Parse cluster rules results
all_res = lapply(res_names, function(r) {
  print(r)
  
  info = unlist(strsplit(r, "_"))
  n_clocks = as.numeric(info[2])
  n_events = as.numeric(info[3])
  purity = as.numeric(info[4])
  coverage = as.numeric(info[5])
  n_mutations = as.numeric(info[6 ])
  
  sub_dirs = list.files(paste0(res_folder, r))
  
  sub_res = lapply(sub_dirs, function(sub_i) {
    fp = paste0(res_folder, r, "/", sub_i)
    
    if (file.exists(paste0(fp, "/sim.rds"))) {
      
      sim = readRDS(paste0(fp, "/sim.rds"))
      true_clocks = sim$taus_clust
      
      fit = readRDS(paste0(fp, "/res_tickTack_h.rds"))
      
      ms_tibble = fit$results_model_selection$model_selection_tibble
      scores = colnames(ms_tibble)[2:6]
      ms_tibble$LOO = -ms_tibble$LOO 
      ms_tibble$Log_lik = -ms_tibble$Log_lik
      
      score_res = lapply(scores, function(s) {
        best_k = unlist(ms_tibble[,s]) %>% which.min()
        clocks = fit$results$draws_and_summary[[best_k]]$summarized_results$clock_mean
        ri = fossil::rand.index(clocks, true_clocks)
        dplyr::tibble(score = s, RandIndex = ri)  
      }) %>% do.call(bind_rows, .) %>% 
        dplyr::bind_rows(
          dplyr::tibble(
            score = "Entropy",
            RandIndex = fossil::rand.index(fit$results_model_selection$best_fit$summarized_results$clock_mean, true_clocks)
          )
        )
      
      score_res = score_res %>% 
        dplyr::bind_cols(dplyr::tibble(n_clocks=n_clocks, n_events=n_events, purity=purity, coverage=coverage, n_mutations=n_mutations, i.iter = as.numeric(sub_i)))
      score_res
      
    }
  }) %>% do.call("bind_rows", .)
  
  sub_res  
}) %>% do.call("bind_rows", .)

saveRDS(all_res, "results_summarised/scores_randIndex.RDS")