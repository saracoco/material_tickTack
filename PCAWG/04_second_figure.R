
rm(list = ls())
require(tidyverse)
library(patchwork)

pA = readRDS("plot/scatter_Ncna_v_Nclusters_with_HM.rds") + theme(legend.position = "right")
pB = readRDS("plot/HM_scatter_with_smooth.rds")
pC = readRDS("plot/distribution_of_HM_fraction_per_ttype.rds")
pE = readRDS("plot/scatters/gene_level_events.rds")

des = "
AABC
AABC
DDEE
FFEE
FFEE"

p = free(pA) + free(pB) + free(pC) + plot_spacer()  + free(pE) + plot_spacer() +
  plot_layout(design = des) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold")
  )
p
ggsave("plot/main2.pdf", plot = p, width = 12, height = 14, units = "in", dpi = 900)
