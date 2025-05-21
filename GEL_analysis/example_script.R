library(tickTack)
require(dplyr)

# View example dataset components
mutations <- tickTack::pcawg_example$mutations
cna <- tickTack::pcawg_example$cna
metadata <- tickTack::pcawg_example$metadata

# Extract input data
segments <- tickTack::pcawg_example$cna
mutations <- tickTack::pcawg_example$mutations
purity <- tickTack::pcawg_example$metadata$purity

# Run the fit function
results <- fit(
  segments = segments,
  mutations = mutations,
  purity = purity,
  possible_k = c("2:1", "2:2", "2:0"),
  beta_binomial = TRUE
)

plot = tickTack::plot_timing(results, segments, colour_by = "karyotype")
ggplot2::ggsave(plot, filename='home/example_plot.png')
