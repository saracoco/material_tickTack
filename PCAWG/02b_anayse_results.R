
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

RES = readRDS("results/summary_all_samples.rds")

# Look at statistics
RES %>% 
  dplyr::select(ttype, sample_id) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(n) %>% 
  dplyr::mutate(ttype = factor(ttype, levels=ttype)) %>% 
  ggplot2::ggplot(mapping = aes(x=ttype, y=n)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  theme_bw()

ggsave("plot/ttype_distributions.png", width = 10, height =10, units="in", dpi=300)

RES$ttype %>% unique()

RES %>% 
  dplyr::group_by(ttype, karyotype) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::mutate(f = n / sum(n)) %>% 
  ggplot2::ggplot(mapping = aes(x=karyotype, y=n)) +
  ggplot2::geom_col() +
  facet_wrap(.~ttype, scales = "free_y") +
  theme_bw()

# Classify each segment by Gene
res_w_drivers = readRDS("results/res_w_onco_and_ts.rds")

driver_order_by_clock = res_w_drivers %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(median_clock = median(clock_mean)) %>% 
  dplyr::arrange(-median_clock) %>% 
  dplyr::pull(gene)

res_w_drivers %>% 
  dplyr::mutate(gene = factor(gene, levels = c(driver_order_by_clock))) %>% 
  ggplot(mapping = aes(x=gene, y=clock_mean, col=type)) +
  geom_boxplot() +
  facet_wrap(~ttype) +
  labs(x = "Driver gene", y = "Mean clock") +
  theme_bw()
  
# Select tumour types
tumour_type = "CLLE"

res_w_drivers %>% 
  dplyr::filter(ttype == tumour_type) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(n_clock = length(unique(clock_mean))) %>% 
  dplyr::pull(n_clock) %>% 
  table()

# Filter drivers seen in at least a certain number of samples
frequent_drivers <- res_w_drivers %>% 
  dplyr::filter(ttype == tumour_type) %>% 
  #dplyr::group_by(ttype, driver) %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= 3) %>% 
  dplyr::pull(gene)
frequent_drivers

# Create all pairs of drivers
driver_pairs <- expand.grid(first_driver = frequent_drivers, second_driver = frequent_drivers) %>%
  dplyr::mutate(first_driver = as.character(first_driver), second_driver=as.character(second_driver)) %>% 
  dplyr::filter(first_driver <= second_driver) # Avoid duplicate/reverse pairs

# Precompute necessary data from res_w_drivers
res_processed <- res_w_drivers %>%
  group_by(sample_id) %>%
  mutate(driver_in_pair = gene %in% frequent_drivers) %>%
  filter(driver_in_pair) %>%
  ungroup()

# Initialize results as a list for better performance
results <- vector("list", nrow(driver_pairs))

# Calculate scores for each pair
for (k in seq_len(nrow(driver_pairs))) {
  print(k)
  pair <- driver_pairs[k, ]
  
  score <- res_processed %>%
    filter(gene %in% c(pair$first_driver, pair$second_driver)) %>%
    group_by(sample_id) %>%
    mutate(n = n()) %>%
    filter(n != 1) %>%
    summarise(score = mean(ifelse(clock_rank[gene == pair$first_driver] < clock_rank[gene == pair$second_driver], 
                                  -1, 
                                  ifelse(clock_rank[gene == pair$first_driver] == clock_rank[gene == pair$second_driver], 0, 1)))) %>%
    pull(score)
  
  n_samples = length(score)
  score = mean(score)
  
  # Store result
  results[[k]] <- tibble(first_driver = pair$first_driver, second_driver = pair$second_driver, score = score, n_samples=n_samples)
}

# Combine results into a single dataframe
scores_df <- bind_rows(results)

scores_df %>% 
  na.omit() %>% 
  ggplot(mapping = aes(x=first_driver, y=second_driver, fill=score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white") +
  theme_bw() +
  labs(x = "Driver 1", y = 'Driver 2') +
  ggtitle(tumour_type)


ggsave(paste0("plot/matrices/", tumour_type, ".png"), width = 10, height =10, units="in", dpi=300)
