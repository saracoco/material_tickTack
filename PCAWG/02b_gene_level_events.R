
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
saveRDS(p_sample_distribution_over_ttype, "plot/ttype_distributions.rds")
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
saveRDS(p_karyotype_distribution_over_ttype, "plot/p_karyotype_distribution_over_ttype.rds")
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
saveRDS(p_karyotype_per_sample_dist, "plot/p_karyotype_per_sample_dist.rds")
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
saveRDS(p_k_per_gene_dist, "plot/p_k_per_gene_dist.rds")
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
saveRDS(p_k_per_gene_dist_frac, "plot/p_k_per_gene_dist_frac.rds")
ggsave("plot/p_k_per_gene_dist_frac.pdf", width = 10, height =10, units="in", dpi=300)

# Select tumour types
unique(res_w_drivers$ttype)

scores_df_all = dplyr::tibble()
for (tumour_type in unique(res_w_drivers$ttype)) {
  print(tumour_type)
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
    
    max_samples = length(unique(res_processed$sample_id))
    
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
        samples = c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post))
        if (all(samples == 0)) {
          p_value = 1
        } else if (all(samples == -1) | all(samples == 1)) {
          p_value = .001
        } else {
          t_test <- t.test(c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post)), mu = 0)
          p_value = t_test$p.value  
        }
        
        score = mean(samples)
        n_samples = sum(O_counts)
        
        # Store result
        return(tibble(first_driver = pair$first_driver, second_driver = pair$second_driver, score = score, n_samples=n_samples, p.value = p_value))
      }
    }) %>% do.call("bind_rows", .) %>% 
      dplyr::mutate(ttype = tumour_type, max_samples=max_samples)
    
    scores_df_all <- dplyr::bind_rows(scores_df_all, scores_df)
    
    if (nrow(scores_df) > 0) {
      scores_df = scores_df %>% 
        dplyr::filter(p.value <= .05) %>% 
        na.omit()
      
      if (nrow(scores_df) > 0) {
        # mat = scores_df %>% 
        #   dplyr::filter(p.value <= .05) %>% 
        #   na.omit() %>% 
        #   ggplot(mapping = aes(x=first_driver, y=second_driver, fill=score)) +
        #   geom_tile() +
        #   geom_text(aes(label=n_samples)) +
        #   scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
        #   theme_bw() +
        #   labs(x = "Driver 1", y = 'Driver 2') +
        #   ggtitle(tumour_type)
        # 
        # ggsave(paste0("plot/matrices/", tumour_type, ".png"), width = 10, height =10, units="in", dpi=300, plot = mat)
        
        scores_df <- dplyr::bind_rows(
          scores_df,
          dplyr::tibble(first_driver=scores_df$first_driver[1], second_driver="", score=0.0, n_samples = .0, p.value = .01)
        )
        
        scatterp = scores_df %>% 
          dplyr::filter(p.value <= .05) %>% 
          #na.omit() %>% 
          dplyr::mutate() %>% 
          ggplot(mapping = aes(x = score, y=first_driver, fill=score, size=n_samples, col=score, label=second_driver)) +
          geom_point() +
          theme_bw() +
          ggrepel::geom_label_repel(col="black", size=4, fill=alpha('white', .1), min.segment.length = 0, box.padding = .5) +
          lims(x = c(-1,1)) +
          scale_color_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
          scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          labs(x = "Score", y="Driver", fill="Score", size="N samples", col="Score") +
	        ggtitle(tumour_type)

        saveRDS(scatterp, paste0("plot/scatters/", tumour_type, ".rds"))
  	    ggsave(paste0("plot/scatters/", tumour_type, ".png"), width = 5, height =5, units="in", dpi=300, plot = scatterp)
  	    
      }
    }
  }
}

TSGs <- c("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC", "CDKN2A", "SMAD4", "VHL", "NF1")
Oncogenes <- c("MYC", "KRAS", "BRAF", "EGFR", "HER2", "ALK", "PIK3CA", "ABL1", "CCND1", "NRAS")
genes_of_interest <- c(TSGs, Oncogenes)

p = scores_df_all %>% 
  dplyr::mutate(class = ifelse(second_driver %in% Oncogenes, "Oncogene", "TSG")) %>% 
  dplyr::mutate(fraction = n_samples / max_samples) %>% 
  dplyr::filter(p.value <= .05) %>% 
  #dplyr::mutate() %>% 
  ggplot(mapping = aes(x = score, y=first_driver, size=fraction, col=class, label=second_driver)) +
  geom_point() +
  theme_bw() +
  ggrepel::geom_label_repel(col="black", size=4, fill=alpha('white', .1), min.segment.length = 0, box.padding = .5) +
  lims(x = c(-1,1)) +
  #scale_color_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
  #scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Score", y="Driver", fill="Score", size="Sample fraction", col="") +
  facet_grid(ttype~., scales = "free", space="free") +
  scale_color_manual(values = list("TSG" = alpha("#006663", .7), "Oncogene"=alpha("#b4662a", .7)))
p
ggsave(filename = "plot/scatters/gene_level_events.pdf", width = 8, height = 10, units = 'in', plot = p)
saveRDS(p, "plot/scatters/gene_level_events.rds")

# scores_df_all %>% 
#   dplyr::mutate(fraction = n_samples / max_samples) %>% 
#   #dplyr::filter(fraction >= .1) %>% 
#   dplyr::mutate(score = abs(score)) %>% 
#   dplyr::mutate(class = ifelse(p.value >= .05, "Not significant", "Significant")) %>% 
#   ggplot(mapping = aes(x=score, fill = class)) +
#   geom_histogram()

p <- scores_df_all %>% 
  dplyr::mutate(fraction = n_samples / max_samples) %>% 
  #dplyr::filter(fraction >= .1) %>% 
  #dplyr::mutate(score = abs(score)) %>% 
  dplyr::mutate(class = ifelse(p.value >= .05, " ", ttype)) %>% 
  dplyr::mutate(lab = if_else(p.value <= .05, paste0(first_driver, "-", second_driver), "")) %>% 
  ggplot(mapping = aes(x=score, y = -log10(p.value), size=-log10(p.value), col=class, label=lab)) +
  geom_point() +
  geom_vline(xintercept = c(.1, -.1), linetype="dashed") +
  geom_hline(yintercept = -log10(.05), linetype="dashed") +
  ggrepel::geom_label_repel(col="black", size=4, fill=alpha('white', .99), min.segment.length = 0, box.padding = .5) +
  theme_bw() +
  #scale_color_manual(values = list("Not signif."= "grey70", "Signif."="indianred")) +
  labs(x = "Score", y=bquote(-log[10] ~ pvalue), col="") +
  scale_color_manual(
    values = list(
     " " = "gray90",
     "BOCA" = rgb(203/255, 204/255, 250/255, alpha = 1),
     "BRCA" = rgb(246/255, 196/255, 205/255, alpha = 1),
     "ESAD" = rgb(136/255, 181/255, 215/255, alpha = 1),
     "MELA" = rgb(255/255, 255/255, 167/255, alpha = 1),
     "PACA" = rgb(247/255, 216/255, 133/255, alpha = 1),
     "PRAD" = '#9ec6b3ff'
    )
  ) +
  guides(size="none")
p  

ggsave(filename = "plot/volcano_plot_scores.pdf", width = 8, height = 8, units = 'in', plot = p)
saveRDS(p, "plot/volcano_plot_scores.rds")
