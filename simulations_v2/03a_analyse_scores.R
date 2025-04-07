.libPaths("~/R/rstudio_v3/")
rm(list = ls())
require(tidyverse)
require(patchwork)
source("utils_plot.R")

all_res = readRDS("results_summarised/scores_randIndex_test.RDS")
all_res %>% 
  ggplot(mapping = aes(x=score, y=RandIndex, col=score)) +
  geom_boxplot() +
  ggh4x::facet_nested("N clocks"+n_clocks~ "N events"+n_events) +
  theme_bw()

#ggsave("plot/scores_ri_dists.pdf", width = 10, height = 10, units = "in", dpi=600)

df_scores = lapply(seq(.01, .99, length = 20) , function(q) {
  all_res %>% 
    dplyr::group_by(score, n_clocks, n_events) %>% 
    dplyr::summarise(m = stats::quantile(RandIndex, q), .groups = "drop") %>% 
    dplyr::group_by(n_clocks, n_events) %>% 
    dplyr::arrange(-m, .by_group = TRUE) %>%
    dplyr::mutate(rank = dense_rank(-m)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(score) %>% 
    dplyr::filter(score=="AIC"| score=="ICL" | score=="BIC" | score=="LOO" | score=="Log_lik")%>%
    dplyr::summarise(avg_rank = mean(rank)) %>% 
    dplyr::arrange(avg_rank) %>% 
    dplyr::mutate(q = q)
}) %>% do.call("bind_rows", .)

p1 = df_scores %>% 
  ggplot(mapping = aes(x=q, y=avg_rank, col=score)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Quantile q", y="Average rank") +
  ggtitle("Cluster 1 to 5") +
  labs(col = "Criterion")
p1
ggsave("plot/ranks_1_to_5.pdf", width = 10, height = 10, units = "in", dpi=600)


df_scores = lapply(seq(.01, .99, length = 20) , function(q) {
  all_res %>% 
    dplyr::filter(n_clocks > 2) %>% 
    dplyr::group_by(score, n_clocks, n_events) %>% 
    dplyr::summarise(m = stats::quantile(RandIndex, q), .groups = "drop") %>% 
    dplyr::group_by(n_clocks, n_events) %>% 
    dplyr::arrange(-m, .by_group = TRUE) %>%
    dplyr::mutate(rank = dense_rank(-m)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(score) %>% 
    dplyr::filter(score=="AIC"| score=="ICL" | score=="BIC" | score=="LOO" | score=="Log_lik")%>%
    dplyr::summarise(avg_rank = mean(rank)) %>% 
    dplyr::arrange(avg_rank) %>% 
    dplyr::mutate(q = q)
}) %>% do.call("bind_rows", .)

p2 = df_scores %>% 
  ggplot(mapping = aes(x=q, y=avg_rank, col=score)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Quantile q", y="Average rank") +
  ggtitle("Cluster 3 to 5") +
  labs(col = "Criterion")
p2
ggsave("plot/ranks_3_to_5.pdf", width = 10, height = 10, units = "in", dpi=600)


p12 = p1 + p2 +
  plot_layout(ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold")
  )
ggsave("plot/model_selection_ranks.pdf", width = 10, height = 10, units = "in", dpi=600, plot=p12)
