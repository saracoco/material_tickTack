
rm(list = ls())
require(tidyverse)

# Read data ####
whole_res <- readRDS("processed_results/simulations_results.rds")

whole_res %>% 
  tidyr::pivot_longer(c(time_tickTack, time_MutTime, singleTT)) %>%
  dplyr::mutate(percent_error = abs(value - real_clocks) / real_clocks) %>%
  dplyr::filter(n_muts <= 30) %>% 
  ggplot(mapping = aes(x=name, y=percent_error, col=name)) +
  geom_boxplot() +
  facet_grid(coverage ~ purity) +
  scale_y_continuous(transform = "log10")


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


