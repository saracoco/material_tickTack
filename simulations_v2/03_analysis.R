
rm(list = ls())
require(tidyverse)

all_res = readRDS("results_summarised/whole_res.RDS")
all_res %>% 
  dplyr::group_by(n_clocks, n_events, n_mutations, purity, coverage) %>% 
  dplyr::summarise(n = length(unique(i.iter))) %>% 
  dplyr::pull(n) %>% 
  table()

all_res %>% 
  dplyr::group_by(n_clocks, n_events, n_mutations, purity, coverage, i.iter) %>% 
  dplyr::summarise(
    rmse_AmpTimeR = sqrt(mean((true_tau - tau_AmpTimeR)**2)),
    rmse_MutTimeR = sqrt(mean((true_tau - tau_MutTimeR)**2)),
    rmse_tickTack = sqrt(mean((true_tau - tau_tickTack)**2)),
    rmse_tickTack_h = sqrt(mean((true_tau - tau_tickTack_h)**2))
  ) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(!c(n_clocks, n_events, n_mutations, purity, coverage, i.iter)) %>% 
  dplyr::group_by(n_mutations, purity, coverage, name, n_clocks) %>% 
  dplyr::summarise(mean_rmse = mean(value), sd_rmse = sd(value)) %>% 
  dplyr:::ungroup() %>% 
  ggplot(mapping = aes(x=n_mutations, y=mean_rmse, ymin=mean_rmse-sd_rmse, ymax=mean_rmse+sd_rmse, col=name)) +
  geom_pointrange() +
  geom_line() +
  ggh4x::facet_nested(n_clocks+coverage~purity) +
  theme_bw()

all_res %>% 
  tidyr::pivot_longer(c(tau_AmpTimeR, tau_MutTimeR, tau_tickTack, tau_tickTack_h)) %>% 
  dplyr::group_by(n_clocks, n_events, purity, coverage, n_mutations, i.iter, name) %>% 
  dplyr::summarise(rmse = sqrt(mean(value - true_tau)**2)) %>% 
  ggplot(mapping = aes(x = name, y=rmse, col=name)) +
  geom_boxplot() +
  facet_grid(n_events~n_clocks) +
  theme_bw()

# Clutering results plots
res_clust = readRDS("results_summarised/clustering_results.RDS")
res_clust %>% 
  dplyr::select(n_clocks, n_events, ri_at, ri_mt, ri_tt, ri_tth) %>% 
  tidyr::pivot_longer(c(ri_at, ri_mt, ri_tt, ri_tth)) %>% 
  ggplot(mapping = aes(x = name, y=value, col=name)) +
  geom_boxplot() +
  facet_grid(n_events~n_clocks) +
  theme_bw()
