
.libPaths(new="~/R/rstudio/") 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("utils.R")

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  select(sample_id, ttype) %>% 
  dplyr::distinct()

results_path <- "/orfeo/scratch/cdslab/scocomello/material_tickTack/PCAWG/results_tickTack_parallel_2/"

IDs <- list.files(results_path)

RES <- lapply(IDs, function(id) {
  results_model_selection_path = paste0(results_path, id, "/results/results_model_selection.rds")
  if (!file.exists(results_model_selection_path)) {
    print(paste0("Skipping sample id ", id))
    return(NULL)
  }
  res_model_selection <- readRDS(results_model_selection_path)
  results <- res_model_selection$best_fit$summarized_results
  
  df <- ttypes %>% 
    dplyr::filter(sample_id == id)
  
  if (nrow(df) == 0) {
    df = dplyr::tibble(sample_id=id, ttype="Uknown")
  }
  
  
  dplyr::bind_cols(df, parse_summarized_results(results))
}) %>% do.call("bind_rows", .)

saveRDS(RES, "results/summary_all_samples.rds")

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
add_TSGs_and_Oncogenes_to_results <- function(RES) {
  load("data/gene_coordinates_hg19.rda")
  
  TSGs <- c("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC", "CDKN2A", "SMAD4", "VHL", "NF1")
  Oncogenes <- c("MYC", "KRAS", "BRAF", "EGFR", "HER2", "ALK", "PIK3CA", "ABL1", "CCND1", "NRAS")
  genes_of_interest <- c(TSGs, Oncogenes)
  
  i = 1
  new_res <- lapply(1:nrow(RES), function(i) {
    print(i)
    seg <- RES[i,]  
    
    genes_found <- gene_coordinates_hg19 %>% 
      dplyr::filter(chr == seg$chr) %>% 
      dplyr::filter(from >= seg$from, to <= seg$to) %>% 
      dplyr::pull(gene)
    
    driver_genes <- genes_found[genes_found %in% genes_of_interest]
    
    if (length(driver_genes) == 0) {
      driver_genes <- c("None")
    }
    
    dplyr::bind_cols(seg, dplyr::tibble(driver = driver_genes))
  }) %>% do.call("bind_rows", .)
  
  new_res
  
}

add_drivers_per_ttype_to_results <- function(RES) {
  load("data/gene_coordinates_hg19.rda")
  
  metadata_samples <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t")
  
  drivers_per_ttype = list()
  for (t in unique(metadata_samples$ttype)) {
    drivers <- metadata_samples %>% 
      dplyr::filter(ttype == t) %>% 
      dplyr::pull(gene) %>% 
      unique()
    drivers_per_ttype[[t]] = drivers
  }  
  
  new_res <- lapply(1:nrow(RES), function(i) {
    print(i)
    seg <- RES[i,]  
    
    genes_found <- gene_coordinates_hg19 %>% 
      dplyr::filter(chr == seg$chr) %>% 
      dplyr::filter(from >= seg$from, to <= seg$to) %>% 
      dplyr::pull(gene)
    
    driver_genes <- genes_found[genes_found %in% drivers_per_ttype[[seg$ttype]]]
    
    if (length(driver_genes) == 0) {
      driver_genes <- c("None")
    }
    
    dplyr::bind_cols(seg, dplyr::tibble(driver = driver_genes))
  }) %>% do.call("bind_rows", .)
  
  new_res
}


#res_w_drivers = add_drivers_to_results(RES)
res_w_drivers = add_TSGs_and_Oncogenes_to_results(RES)
res_w_drivers = res_w_drivers %>% dplyr::filter(driver != "None")

driver_order_by_clock = res_w_drivers %>% 
  dplyr::group_by(driver) %>% 
  dplyr::summarise(median_clock = median(clock_mean)) %>% 
  dplyr::arrange(-median_clock) %>% 
  dplyr::pull(driver)

res_w_drivers %>% 
  dplyr::mutate(driver = factor(driver, levels = c(driver_order_by_clock))) %>% 
  ggplot(mapping = aes(x=driver, y=clock_mean)) +
  geom_boxplot() +
  labs(x = "Driver gene", y = "Mean clock") +
  theme_bw()
  

# Select tumour types
# tumour_type = "Liver-HCC"
# res_w_drivers <- res_w_drivers %>% 
#   dplyr::filter(driver != "None", ttype == tumour_type)

# Filter drivers seen in at least a certain number of samples
frequent_drivers <- res_w_drivers %>% 
  #dplyr::group_by(ttype, driver) %>% 
  dplyr::group_by(driver) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= 3) %>% 
  dplyr::pull(driver)

# Create all pairs of drivers
driver_pairs <- expand.grid(first_driver = frequent_drivers, second_driver = frequent_drivers) %>%
  dplyr::mutate(first_driver = as.character(first_driver), second_driver=as.character(second_driver)) %>% 
  dplyr::filter(first_driver <= second_driver) # Avoid duplicate/reverse pairs

# Precompute necessary data from res_w_drivers
res_processed <- res_w_drivers %>%
  group_by(sample_id) %>%
  mutate(driver_in_pair = driver %in% frequent_drivers) %>%
  filter(driver_in_pair) %>%
  ungroup()

# Initialize results as a list for better performance
results <- vector("list", nrow(driver_pairs))

# Calculate scores for each pair
for (k in seq_len(nrow(driver_pairs))) {
  print(k)
  pair <- driver_pairs[k, ]
  
  score <- res_processed %>%
    filter(driver %in% c(pair$first_driver, pair$second_driver)) %>%
    group_by(sample_id) %>%
    mutate(n = n()) %>%
    filter(n != 1) %>%
    summarise(score = mean(ifelse(clock_rank[driver == pair$first_driver] < clock_rank[driver == pair$second_driver], 
                                  -1, 
                                  ifelse(clock_rank[driver == pair$first_driver] == clock_rank[driver == pair$second_driver], 0, 1))),
              .groups = "drop") %>%
    pull(score)
  
  # Store result
  results[[k]] <- tibble(first_driver = pair$first_driver, second_driver = pair$second_driver, score = score)
}

# Combine results into a single dataframe
scores_df <- bind_rows(results)

scores_df %>% 
  na.omit() %>% 
  ggplot(mapping = aes(x=first_driver, y=second_driver, fill=score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white") +
  theme_bw()
