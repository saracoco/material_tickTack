
rm(list = ls())
require(tidyverse)
require(patchwork)
source("utils_plot.R")

all_res = readRDS("results_summarised/whole_res.RDS") %>%
  dplyr::mutate(setting = paste(n_clocks, n_events, purity, coverage, n_mutations, i.iter, sep = "_"))

res_clust = readRDS("results_summarised/clustering_results.RDS")

# Plots
error_over_nmuts = plot_over_nmutations(all_res)
error_over_clocks_and_muts = plot_nclocks_v_nevents(all_res)
rand_index_plot = plot_rand_index(res_clust)

design = "
A
A
A
B
B
C
C
"

sim_plot = rand_index_plot + error_over_nmuts + error_over_clocks_and_muts +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold")
  )
sim_plot
ggsave("plot/sim_plot.pdf", units = "in", dpi = 600, width = 12, height = 12, plot = sim_plot)
