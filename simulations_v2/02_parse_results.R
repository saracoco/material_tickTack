
rm(list = ls())

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

