
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

class_colors = list(
  'Classic' = "#91cf60",
  'Hopeful Monster' = "#fee08b"
)

k_colors = list(
  '2:0' = 'turquoise4',
  '2:1' = ggplot2::alpha('orange', .8),
  '2:2' = 'firebrick3'  
)

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

RES = readRDS("results/summary_all_samples.rds")

RES %>% 
  dplyr::group_by(ttype, sample_id) %>% 
  dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>% 
  dplyr::select(sample_id, n_events, n_clusters) %>% 
  dplyr::distinct() %>% 
  ggplot(mapping = aes(x=n_events, y=n_clusters)) +
  geom_point() +
  theme_bw()

RES %>% 
  dplyr::group_by(ttype, sample_id) %>% 
  dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>% 
  dplyr::select(sample_id, n_events, n_clusters) %>% 
  dplyr::distinct() %>% 
  ggplot(mapping = aes(x=n_events, y=n_clusters)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ttype)

q9 = data %>% 
  dplyr::mutate(x = n_events / n_clusters) %>% 
  dplyr::pull(x) %>% 
  stats::quantile(.8)

rr = RES %>% 
  dplyr::group_by(ttype, sample_id) %>% 
  dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>% 
  dplyr::select(sample_id, n_events, n_clusters) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(f = n_events / n_clusters) %>% 
  dplyr::mutate(class = ifelse(f >= q9, "Hopeful Monster", "Classic")) 

rr %>% 
  ggplot(mapping = aes(x=n_events, y=n_clusters, col=class)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ttype) +
  theme(legend.position = "bottom") +
  labs(x = "N events", y="N clusters", col="") +
  scale_color_manual(values = class_colors)

rr %>% 
  dplyr::group_by(ttype, class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::mutate(n = n / sum(n)) %>% 
  ggplot(mapping = aes(x=reorder(ttype, -n), y=n, fill=class)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  coord_flip()
  