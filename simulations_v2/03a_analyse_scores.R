
rm(list = ls())
require(tidyverse)
require(patchwork)
source("utils_plot.R")

all_res = readRDS("results_summarised/scores_randIndex.RDS")
all_res %>% 
  ggplot(mapping = aes(x=score, y=RandIndex, col=score)) +
  geom_boxplot() +
  ggh4x::facet_nested("N clocks"+n_clocks~ "N events"+n_events) +
  theme_bw()

ggsave("plot/scores_ri_dists.pdf", width = 10, height = 10, units = "in", dpi=600)

df_scores = lapply(seq(.01, .99, length = 20) , function(q) {
  all_res %>% 
    dplyr::group_by(score, n_clocks, n_events) %>% 
    dplyr::summarise(m = stats::quantile(RandIndex, q), .groups = "drop") %>% 
    dplyr::group_by(n_clocks, n_events) %>% 
    dplyr::arrange(-m, .by_group = TRUE) %>%
    dplyr::mutate(rank = dense_rank(-m)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(score) %>% 
    dplyr::summarise(avg_rank = mean(rank)) %>% 
    dplyr::arrange(avg_rank) %>% 
    dplyr::mutate(q = q)
}) %>% do.call("bind_rows", .)

df_scores %>% 
  ggplot(mapping = aes(x=q, y=avg_rank, col=score)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Quantile q", y="Average rank")

all_res %>% 
  dplyr::filter(n_clocks > 1) %>% 
  dplyr::group_by(score, n_clocks, n_events, purity, coverage, n_mutations) %>% 
  dplyr::summarise(m = median(RandIndex), .groups = "drop") %>% 
  dplyr::group_by(n_clocks, n_events, purity, coverage, n_mutations) %>% 
  dplyr::filter(m == max(m)) %>% 
  dplyr::pull(score) %>% 
  table()
