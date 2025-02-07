
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

# Look at statistics
p_sample_distribution_over_ttype <- RES %>% 
  dplyr::select(ttype, sample_id) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(ttype) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(n) %>% 
  dplyr::mutate(ttype = factor(ttype, levels=ttype)) %>% 
  ggplot2::ggplot(mapping = aes(x=ttype, y=n)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  theme_bw() +
  labs(x = "Tumour type", y = "N sample") +
  ggtitle(paste0("Total samples = ", RES$sample_id %>% unique() %>% length()))
p_sample_distribution_over_ttype
ggsave("plot/ttype_distributions.pdf", width = 10, height =10, units="in", dpi=300)

p_karyotype_distribution_over_ttype <- RES %>% 
  dplyr::group_by(ttype, karyotype) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(mapping = aes(x=ttype, y=n, fill = karyotype)) +
  geom_bar(position="dodge", stat="identity") +
  ggplot2::coord_flip() +
  theme_bw() +
  labs(x = "Tumour type", y = "N events", fill="") +
  ggtitle("Number of simple CNA events by karyotype and tyumour type") +
  scale_fill_manual(values = k_colors)
p_karyotype_distribution_over_ttype
ggsave("plot/karyotype_distribution_over_ttype.pdf", width = 10, height =10, units="in", dpi=300)

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
  labs(x = "Tumour type", y = "N events per sample", fill="") +
  ggtitle("Number of simple CNA events by karyotype and tyumour type") +
  scale_fill_manual(values = k_colors)
p_karyotype_per_sample_dist
ggsave("plot/karyotype_per_sample_dist.pdf", width = 10, height =10, units="in", dpi=300)

# Classify each segment by Gene
res_w_drivers = readRDS("results/res_w_onco_and_ts.rds")

# Plot type of karyotypes by type of event
p_k_per_gene_dist <- res_w_drivers %>% 
  dplyr::filter(type != "Multiple") %>% 
  dplyr::group_by(type, karyotype) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::select(type, karyotype, n) %>% 
  dplyr::distinct() %>% 
  ggplot(mapping = aes(x=type, y=n, fill=karyotype)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = k_colors) +
  labs(x = "Gene type", y = "N events", fill = "CN") + 
  ggtitle("N events per gene type")
p_k_per_gene_dist
ggsave("plot/p_k_per_gene_dist.pdf", width = 10, height =10, units="in", dpi=300)

p_k_per_gene_dist_frac <- res_w_drivers %>% 
  dplyr::filter(type != "Multiple") %>% 
  dplyr::group_by(type, karyotype) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::select(type, karyotype, n) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(type) %>% 
  dplyr::mutate(f = n / sum(n)) %>% 
  ggplot(mapping = aes(x=type, y=f, fill=karyotype)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = k_colors) +
  labs(x = "Gene type", y = "F events", fill = "CN") + 
  ggtitle("Fraction of events per gene type")
p_k_per_gene_dist_frac
ggsave("plot/p_k_per_gene_dist_frac.pdf", width = 10, height =10, units="in", dpi=300)

# Select tumour types

p_cluster_v_events_ttype <- RES %>% 
  dplyr::group_by(ttype, sample_id) %>% 
  dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>% 
  dplyr::select(sample_id, n_events, n_clusters) %>% 
  dplyr::distinct() %>% 
  ggplot(mapping = aes(x=n_events, y=n_clusters)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ttype)
p_cluster_v_events_ttype
ggsave("plot/p_cluster_v_events_ttype.pdf", width = 10, height =10, units="in", dpi=300)

p_drivers_cluster_v_events_ttype <- res_w_drivers %>% 
  dplyr::group_by(ttype, sample_id) %>% 
  dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>% 
  dplyr::select(sample_id, n_events, n_clusters) %>% 
  dplyr::distinct() %>% 
  ggplot(mapping = aes(x=n_events, y=n_clusters)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ttype) +
  ggtitle("Drivers only")
p_drivers_cluster_v_events_ttype
ggsave("plot/p_drivers_cluster_v_events_ttype.pdf", width = 10, height =10, units="in", dpi=300)

p_tau_dist_driver_events_ttype <- res_w_drivers %>% 
  #dplyr::filter(ttype == tumour_type) %>% 
  ggplot2::ggplot(mapping = aes(x=type, y=clock_mean)) +
  geom_violin() +
  theme_bw() +
  facet_wrap(~ttype) +
  coord_flip()
p_tau_dist_driver_events_ttype
ggsave("plot/p_tau_dist_driver_events_ttype.pdf", width = 10, height =10, units="in", dpi=300)

res_w_drivers %>% 
  #dplyr::filter(ttype == tumour_type) %>% 
  ggplot2::ggplot(mapping = aes(x=type, y=clock_mean)) +
  geom_violin() +
  theme_bw() +
  facet_wrap(~ttype) +
  coord_flip()

res_w_drivers %>% 
  #dplyr::filter(ttype == tumour_type) %>% 
  ggplot2::ggplot(mapping = aes(x=karyotype, y=clock_mean)) +
  geom_violin() +
  theme_bw() +
  facet_wrap(~ttype)
  

# Filter drivers seen in at least a certain number of samples
unique(res_w_drivers$ttype)

tumour_type = "PRAD"
k = 3

for (tumour_type in unique(res_w_drivers$ttype)) {
  frequent_drivers <- res_w_drivers %>% 
    dplyr::filter(ttype == tumour_type) %>% 
    #dplyr::group_by(ttype, driver) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= 3) %>% 
    dplyr::pull(gene)
  
  if (length(frequent_drivers) > 0) {
    # Create all pairs of drivers
    driver_pairs <- expand.grid(first_driver = frequent_drivers, second_driver = frequent_drivers) %>%
      dplyr::mutate(first_driver = as.character(first_driver), second_driver=as.character(second_driver)) %>% 
      dplyr::filter(first_driver <= second_driver) # Avoid duplicate/reverse pairs
    
    # Precompute necessary data from res_w_drivers
    res_processed <- res_w_drivers %>%
      dplyr::filter(ttype == tumour_type) %>% 
      group_by(sample_id) %>%
      mutate(driver_in_pair = gene %in% frequent_drivers) %>%
      filter(driver_in_pair) %>%
      ungroup()
    
    # Initialize results as a list for better performance
    results <- vector("list", nrow(driver_pairs))
    
    # Calculate scores for each pair
    scores_df = lapply(1:nrow(driver_pairs), function(k) {
      print(k)
      pair <- driver_pairs[k, ]
      
      co_occurr_df = res_processed %>%
        filter(gene %in% c(pair$first_driver, pair$second_driver)) %>%
        group_by(sample_id) %>%
        mutate(n = n()) %>%
        filter(n != 1) %>% 
        dplyr::select(sample_id, clock_rank, gene) %>% 
        tidyr::pivot_wider(values_from = clock_rank, names_from = gene)
      
      if (nrow(co_occurr_df) > 3) {
        n_pre = sum(co_occurr_df[,pair$first_driver] < co_occurr_df[,pair$second_driver])
        n_coo = sum(co_occurr_df[,pair$first_driver] == co_occurr_df[,pair$second_driver])
        n_post = sum(co_occurr_df[,pair$first_driver] > co_occurr_df[,pair$second_driver])
        
        # Simulated data: Replace these counts with your actual observed data
        O_counts <- c(n_pre, n_coo, n_post)  # Observed frequencies for (-1, 0, 1)
        #DescTools::MultinomCI(O_counts)
        
        # is the mean different from zero?
        t_test <- t.test(c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post)), mu = 0)
        t_test$p.value
        
        score = mean(c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post)))
        n_samples = sum(O_counts)
        
        # Store result
        return(tibble(first_driver = pair$first_driver, second_driver = pair$second_driver, score = score, n_samples=n_samples, p.value = t_test$p.value))
      }
    }) %>% do.call("bind_rows", .)
    
    if (nrow(scores_df) > 0) {
      scores_df = scores_df %>% 
        dplyr::filter(p.value <= .05) %>% 
        na.omit()
      
      if (nrow(scores_df) > 0) {
        mat = scores_df %>% 
          dplyr::filter(p.value <= .05) %>% 
          na.omit() %>% 
          ggplot(mapping = aes(x=first_driver, y=second_driver, fill=score)) +
          geom_tile() +
          geom_text(aes(label=n_samples)) +
          scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
          theme_bw() +
          labs(x = "Driver 1", y = 'Driver 2') +
          ggtitle(tumour_type)
        
        ggsave(paste0("plot/matrices/", tumour_type, ".png"), width = 10, height =10, units="in", dpi=300, plot = mat)    
        
        scatterp = scores_df %>% 
          dplyr::filter(p.value <= .05) %>% 
          na.omit() %>% 
          dplyr::mutate() %>% 
          ggplot(mapping = aes(x = score, y=first_driver, fill=score, size=n_samples, col=score, label=second_driver)) +
          geom_point() +
          theme_bw() +
          ggrepel::geom_label_repel(col="black", size=4, fill=alpha('white', .1), min.segment.length = 0, box.padding = .5) +
          lims(x = c(-1,1)) +
          scale_color_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
          scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          labs(x = "Score", y="Driver", fill="Score", size="N samples", col="Score")
        
        ggsave(paste0("plot/scatters/", tumour_type, ".png"), width = 10, height =10, units="in", dpi=300, plot = scatterp)    
        
        
      }
    }
  }
}


scores_df %>% 
  dplyr::filter(p.value <= .05) %>% 
  na.omit() %>% 
  dplyr::mutate() %>% 
  ggplot(mapping = aes(x = score, y=first_driver, fill=score, size=n_samples, col=score, label=second_driver)) +
  geom_point() +
  theme_bw() +
  ggrepel::geom_label_repel(col="black", size=4, box.padding = unit(-0.5, "lines")) +
  lims(x = c(-1,1)) +
  scale_color_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
  scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed")
