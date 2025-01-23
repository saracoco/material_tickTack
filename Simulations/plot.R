
rm(list = ls())
require(tidyverse)

# Read data ####
s = readRDS("/orfeo/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results/tickTack_sim_1_20_0.8_60/4/input_data.rds")
s$mutations %>%
  group_by(chr, from, to) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::pull(n) %>% 
  hist()

l = (19999999 - 2)
mu = 1e-4
w = 1e-2
dt = 7

rpois(1, l * mu * w * dt)
l * w * dt * mu

res = readRDS("tickTack_generative_model/data/res_single_inference.rds")
results <- lapply(names(res), function(name){
  info = unlist(strsplit(name, "_"))
  n_clocks <- as.numeric(info[2])
  n_segments <- as.numeric(info[3])
  purity <- as.numeric(info[4])
  coverage = as.numeric(info[5])
  iter = as.numeric(unlist(strsplit(info[6], ".rds")))
  r = res[[name]]
  dplyr::bind_cols(
    dplyr::tibble(n_clocks=n_clocks, n_segments=n_segments, purity=purity, coverage=coverage, iter=iter),
    r$mae_ri_results
  )  
}) %>% do.call("bind_rows", .)

results %>% 
  tidyr::pivot_longer(!c(n_clocks, n_segments, purity, coverage, iter)) %>% 
  ggplot(aes(x=name, y=value, col=name)) +
  geom_boxplot() +
  facet_grid(purity~coverage) +
  theme_bw()


