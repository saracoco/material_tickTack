.libPaths("~/R/rstudio_v3/")
rm(list = ls())
require(tidyverse)
require(patchwork)
require(scales)
source("utils_plot.R")

all_res = readRDS("results_summarised/whole_res_test.RDS") %>%
  dplyr::mutate(setting = paste(n_clocks, n_events, purity, coverage, n_mutations, i.iter, sep = "_"))

# res_clust = readRDS("results_summarised/clustering_results_test_multi_test.RDS")
res_clust = readRDS("results_summarised/clustering_results_test.RDS")



# Define filter criteria 
purity_values <- c(0.1, 0.2, 0.3, 0.6, 0.9)
coverage_values <- c(100, 50, 20, 10, 5)
n_clocks_values <- c(1, 2, 3, 4, 5)
n_events_values <- c(5, 10, 20)
n_mutations_values <- c(10, 20, 50, 100)


# purity_values= c(0.1, 0.2, 0.3, 0.6, 0.9)
# coverage_values= c(100, 50, 20, 10, 5)
# n_clocks_values= c(6, 8, 10, 15 )
# n_events_values= c(15, 20, 30, 40)
# n_mutations_values= c(10, 20, 50, 100)


# 
# purity_values=c(0.1, 0.2, 0.3, 0.6, 0.9)
# coverage_values=c(100, 50, 20, 10, 5)
# n_clocks_values=c(6, 8, 10, 15, 20)
# n_events_values=c(6, 8, 10, 15, 20)
# n_mutations_values=c(10, 20, 50, 100)

# 


# 
# # Apply filtering
filtered_res <- all_res %>%
  filter(
    purity %in% purity_values,
    coverage %in% coverage_values,
    n_clocks %in% n_clocks_values,
    n_events %in% n_events_values,
    n_mutations %in% n_mutations_values
  )
# filtered_res <- filtered_res %>%
#     filter((n_clocks == n_events))

all_res <- filtered_res

filtered_clust <- res_clust %>%
  filter(
    purity %in% purity_values,
    coverage %in% coverage_values,
    n_clocks %in% n_clocks_values,
    n_events %in% n_events_values,
    n_mutations %in% n_mutations_values
  )
# filtered_clust <- filtered_clust %>%
#     filter((n_clocks == n_events))
res_clust <- filtered_clust


# Plots
error_over_nmuts = plot_over_nmutations(all_res)
error_over_clocks_and_muts = plot_nclocks_v_nevents(all_res)
error_over_purity = plot_over_purity(all_res)
error_over_coverage = plot_over_coverage(all_res)

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
ggsave("plot/sim_plot_3.pdf", units = "in", dpi = 600, width = 12, height = 12, plot = sim_plot)
