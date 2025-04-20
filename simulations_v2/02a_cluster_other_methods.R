.libPaths("~/R/rstudio_v3/")

rm(list = ls())
source("utils.R")

res_folder = "results/" 
res_names = list.files(res_folder)

all_res = lapply(res_names, function(r) {
  print(r)
  
  info = unlist(strsplit(r, "_"))
  n_clocks = as.numeric(info[2])
  n_events = as.numeric(info[3])
  purity = as.numeric(info[4])
  coverage = as.numeric(info[5])
  n_mutations = as.numeric(info[6 ])
  
  sub_dirs = list.files(paste0(res_folder, r))
  # sub_i = sub_dirs[1]
  
  sub_res = lapply(sub_dirs, function(sub_i) {
    fp = paste0(res_folder, r, "/", sub_i)
    
    if (file.exists(paste0(fp, '/sim.rds'))&
        file.exists(paste0(fp, '/res_AmpTimeR.rds'))& 
        file.exists(paste0(fp, '/res_MutTime.rds'))&
        file.exists(paste0(fp, '/res_tickTack_h.rds'))&
        file.exists(paste0(fp, '/res_tickTack_single.rds'))) {
      sim = readRDS(paste0(fp, '/sim.rds'))
      
      
      
      cl_at_valid_indeces = clusterAmpTimeR(fp)
      cl_at = cl_at_valid_indeces$clusters
      cl_mt = clusterMutTimeR(fp)
      cl_tt = cluster_tickTack_single(fp)
      cl_tth = cluster_tickTack_h(fp)
      
      if (length(cl_at_valid_indeces$na_indices)==0){
        cl_at_full <- rep(NA, length(sim$taus_clust))  
        cl_at_full[-cl_at_valid_indeces$na_indices] <- cl_at 
        cl_at <- cl_at_full
      }
      ri_at = catsim::rand_index(cl_at, sim$taus_clust, na.rm=TRUE)
      ri_mt = catsim::rand_index(cl_mt, sim$taus_clust)
      ri_tt = catsim::rand_index(cl_tt, sim$taus_clust)
      ri_tth = catsim::rand_index(cl_tth, sim$taus_clust)
      
      # if (unique(sim$taus_clust)!=1){
      #   ri_at = fossil::adj.rand.index(cl_at, sim$taus_clust)
      #   ri_mt = fossil::adj.rand.index(cl_mt, sim$taus_clust)
      #   ri_tt = fossil::adj.rand.index(cl_tt, sim$taus_clust)
      #   ri_tth = fossil::adj.rand.index(cl_tth, sim$taus_clust)
      # }
      
      # if (unique(sim$taus_clust)!=1){
        dplyr::tibble(n_clocks, n_events, purity,coverage, n_mutations, i.iter = sub_i, ri_at, ri_mt, ri_tt, ri_tth)  
      # } else {
        # dplyr::tibble(n_clocks, n_events, purity,coverage, n_mutations, i.iter = sub_i, ri_at, ri_mt, ri_tt, ri_tth)  
      # }
    }
  }) %>% do.call("bind_rows", .)
  sub_res
  }) %>% do.call("bind_rows", .)

saveRDS(all_res, "results_summarised/clustering_results_test.RDS")







