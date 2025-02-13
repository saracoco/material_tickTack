
rm(list = ls())
require(tidyverse)
library(patchwork)

#pA = readRDS("plot/scatter_Ncna_v_Nclusters_with_HM.rds") + theme(legend.position = "right")
pA = readRDS("plot/HM_scatter_with_smooth.rds")
pB = readRDS("plot/distribution_of_HM_fraction_per_ttype.rds")
pC = readRDS("plot/GACA_ESAD_gene_incidence_fraction.rds") + theme(legend.position = 'right')
pE = readRDS("plot/scatters/gene_level_events.rds")

des = "
AABCC
DDEEE
FFEEE
FFEEE"

p = free(pA) + free(pB) + free(pC) + plot_spacer()  + free(pE) + plot_spacer() +
  plot_layout(design = des) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold")
  )
p
ggsave("plot/main2.pdf", plot = p, width = 10, height = 10, units = "in", dpi = 900)
