rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

k_colors = list(
  '2:0' = 'turquoise4',
  '2:1' = ggplot2::alpha('orange', .8),
  '2:2' = 'firebrick3'  
)

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

RES = readRDS("results/summary_all_samples.rds")

color = 
  c("LIRI"=rgb(148/255, 131/255, 204/255, alpha = 1),
    "PACA"=rgb(247/255, 216/255, 133/255, alpha = 1),
    "ESAD"=rgb(136/255, 181/255, 215/255, alpha = 1),
    "MALY"=rgb(144/255, 252/255, 206/255, alpha = 1),
    "MELA"=rgb(255/255, 255/255, 167/255, alpha = 1),
    "BRCA"=rgb(246/255, 196/255, 205/255, alpha = 1),
    "BOCA"=rgb(203/255, 204/255, 250/255, alpha = 1),
    'PAEN'='#a7c9e4ff',
    'PBCA'='#db828eff',
    'PRAD'='#9ec6b3ff', 
    'RECA'='#e7ac7fff',
    'OV' = '#d86a89ff',
    'Other (<30)' = 'gainsboro'
    )

# Pie chart 
pie_chart= RES %>% distinct(sample_id, ttype) %>% 
  group_by(ttype) %>% 
  summarise(n=n()) %>% ungroup() %>%
  mutate(ttype = ifelse(n< 30, 'Other (<30)',ttype)) %>% 
  group_by(ttype) %>% 
  summarise(n=sum(n)) %>% ungroup() %>% 
  mutate(label= as.character(n)) %>%
  ggplot(aes(x = "", y = n, fill = ttype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Converts bar chart to pie chart
  theme_void() +  # Removes unnecessary grid and axes
  labs(fill = "", title = "") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5),
            #nudge_x = 0.2,
            color='#4b4d47ff',
            size=2.5)+
  scale_fill_manual(values = color)+ 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 7),  
        legend.title = element_text(size = 10, face = "bold"),  
        legend.key.size = unit(.3, "cm")  ) +
  guides(fill = guide_legend(ncol = 3)) 
pie_chart

# Number of events by ttype
p_karyotype_per_sample_dist <- RES %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::mutate(n_samples = length(unique(sample_id))) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ttype, karyotype) %>% 
  dplyr::mutate(nevents = n()) %>% 
  dplyr::select(ttype, nevents, n_samples, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(events_per_sample = nevents / n_samples) %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::mutate(s = sum(events_per_sample)) %>%  
  ggplot(mapping = aes(x=reorder(ttype, +s), y=events_per_sample, fill = karyotype)) +
  geom_bar(position="dodge", stat="identity") +
  ggplot2::coord_flip() +
  theme_bw() +
  labs(x = "", y = "N events per sample", fill="", size= 5) +
  ggtitle("Average number of simple CNA by tyumour type" ) +
  scale_fill_manual(values = k_colors)+
  xlab('')+ylab("N events per sample")+
  #labs(size=5) +
  theme(
    axis.text = element_text(size=7),
    legend.position = 'right',
    axis.title = element_text(size=8),
    legend.text = element_text(size = 7),  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.key.size = unit(.3, "cm"),
    plot.title = element_text(size= 10)
  ) 
p_karyotype_per_sample_dist

# Figura esempio
eg_1_component <- "0a9c9db0-c623-11e3-bf01-24c6515278c0"
eg_2_component <- "2a8d63eb-0174-4213-9214-413f391f512c"
eg_3_component <- "05616329-e7ba-4efd-87b1-d79cd0f7af3d"
name <- eg_2_component

# x_after_inference_sim <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/Simulations/tickTack_generative_model/results_summary/res_1_5_0.3_20_10_1.rds")
x_after_inference_PCAWG <- readRDS(paste0("~/dati_Orfeo/scocomello/material_tickTack/PCAWG/results_whole_40mut_AIC/", name, "/results/x_after_inference.rds"))
fit <- readRDS(paste0("~/dati_Orfeo/scocomello/data/clonal_analysis_PCAWG/", name, "/fit.rds"))

fit$mutations <- fit$snvs

# x_after_inference_sim$reference_genome <- "GRCh38"
x_after_inference_PCAWG$reference_genome <- fit$reference_genome
x_after_inference_PCAWG$cna_subclonal <- fit$cna_subclonal
x_after_inference_PCAWG$peaks_analysis <- fit$peaks_analysis
x_after_inference_PCAWG$ploidy <- fit$ploidy
x_after_inference_PCAWG$purity <- fit$purity
x_after_inference_PCAWG$n_mutations <- fit$n_snvs
x_after_inference_PCAWG$n_cna_clonal <- fit$n_cna_clonal
x_after_inference_PCAWG$ttype <- unique(fit$snvs$ttype)[2] 


x <- x_after_inference_PCAWG
results<-x$results_timing

K=2
# plot_segments_tick_tack(x)
source("~/Documents/GitHub/material_tickTack/simulations_v2/utils_plot.R")
example = merge_timing_and_segments(x, K=K)+
  theme(
    axis.text = element_text(size=7),
    legend.position = 'bottom',
    legend.text = element_text(size = 7),  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.key.size = unit(.3, "cm"),
    axis.title.y = element_text(size=7)
  )+
  labs(size=7)


######### Assembled figure
fig1 = patchwork::wrap_plots(
  pie_chart,
  p_karyotype_per_sample_dist,
  example,
  design = 
    'AABBB
     AABBB
     CCCCC
     CCCCC
     CCCCC'
) + patchwork::plot_annotation(
  tag_levels = "a"                             
)
fig1 %>% ggsave(filename='/Users/aliceantonello/Documents/GitHub/material_tickTack/PCAWG/plot/PCAWG1.pdf', 
                height = 13, width = 9)








