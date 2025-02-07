setwd('~/Documents/GitHub/material_tickTack/PCAWG')
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("utils.R")

RES = readRDS("results/summary_all_samples.rds")
res_w_drivers = readRDS("results/res_w_onco_and_ts.rds")
arm_events = read.csv('data/annotated_cnas.csv')

# Drivers timing distribution 
p1=ggplot(res_w_drivers %>% group_by(ttype, gene) %>% 
         mutate(counts = n(), gene_with_counts = paste0(gene, '(', n(), ')')) %>%
         filter(counts > 10))+
  geom_boxplot(aes(x =  reorder(gene_with_counts, clock_mean, median), y = clock_mean, fill = type)) +
  #geom_jitter(aes(x =  reorder(gene_with_counts, clock_mean, median), y = clock_mean), size=.1, alpha=.5) +
  facet_wrap(~ttype, scales = 'free_x') + xlab('Drivers')+ylab('Clock')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = 'bottom')

ggsave(p1, 
       filename='~/Documents/GitHub/material_tickTack/PCAWG/plot/drivers_timing_distribution.pdf',
       height = 20, width = 20)

# Arm events timing distribution 
p2=ggplot(arm_events %>% rowwise() %>% 
            mutate(class = strsplit(class, '_')[[1]][1]) 
          %>% group_by(ttype, class) %>% 
         mutate(counts = n(), 
                gene_with_counts = paste0(class, '(', n(), ')')) %>%
         filter(counts > 10))+
  geom_boxplot(aes(x =  reorder(gene_with_counts, clock_mean, median), y = clock_mean)) +
  #geom_jitter(aes(x =  reorder(gene_with_counts, clock_mean, median), y = clock_mean), size=.1, alpha=.5) +
  facet_wrap(~ttype, scales = 'free_x') + xlab('Drivers')+ylab('Clock')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = 'bottom')
ggsave(p2, 
       filename='~/Documents/GitHub/material_tickTack/PCAWG/plot/arm_events_timing_distribution.pdf',
       height = 20, width = 30)





