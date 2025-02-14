
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

N_CHR = 11
class_colors = list(
  'Classic' = "#f1a340",
  'HM' = "#998ec3",
  "WGD" = "#a6dba0"
)

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>%
  dplyr::select(sample_id, ttype) %>%
  dplyr::distinct()

META = readRDS("data/metadata_all.rds") %>% dplyr::select(sample_id, wgd_status)
RES = readRDS("results/summary_all_samples.rds")
RES$sample_id = lapply(RES$sample_id, function(s) {
  str_replace_all(s, ".rds", "")
}) %>% unlist()
RES = RES %>%
  dplyr::left_join(META, by="sample_id")

RES = lapply(unique(RES$sample_id), function(s) {
  RES %>%
    dplyr::filter(sample_id == s) %>%
    dplyr::select(sample_id, ttype, ploidy, ttype, karyotype, clock_rank, wgd_status, chr) %>%
    dplyr::mutate(n_cna = n(), n_clusters = max(clock_rank)) %>%
    dplyr::group_by(clock_rank) %>%
    dplyr::mutate(n_chr_affected = length(unique(chr))) %>%
    dplyr::select(sample_id, ttype, ploidy, wgd_status, n_cna, n_clusters, n_chr_affected, clock_rank) %>%
    dplyr::distinct()
}) %>% do.call("bind_rows", .)

RES = RES %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd" | ploidy == 2) & any(n_chr_affected >= N_CHR), "HM", ifelse(wgd_status == "wgd", "WGD", "Classic"))) %>%
  #dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd") & any(n_chr_affected >= N_CHR), "HM", "Classic")) %>%
  na.omit()

RES = RES %>%
  dplyr::select(!c(n_chr_affected, clock_rank)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct()

# Check GENE in HM and non HM ####
GENE = "NF1"
gene_analysis_df = lapply(genes_of_interest, function(GENE) {
  print(GENE)
  driver_res = readRDS("results/res_w_onco_and_ts.rds") %>%
    dplyr::filter(gene == GENE) %>%
    dplyr::mutate(sample_id = str_replace_all(sample_id, ".rds", ""))

  df_driver = RES %>%
    na.omit() %>%
    dplyr::group_by(is_HM, ttype) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(sample_id %in% driver_res$sample_id) %>%
    dplyr::group_by(is_HM, ttype, n) %>%
    dplyr::summarise(nhm = n()) %>%
    dplyr::mutate(frac = nhm / n) %>%
    dplyr::filter(is_HM %in% c("HM", "Classic"))

  tt = "BRCA"
  gene_res = lapply(unique(df_driver$ttype), function(tt) {
    print(tt)
    d = df_driver %>%
      dplyr::filter(ttype == tt)

    if (nrow(d) == 2) {
      # Create a contingency table
      table_data <- matrix(c(d$nhm, d$n - d$nhm), nrow = 2)
      fish.test = fisher.test(table_data)
      pval = fish.test$p.value

      fHM = d %>% dplyr::filter(is_HM == "HM") %>% dplyr::pull(frac)
      fClassic = d %>% dplyr::filter(is_HM == "Classic") %>% dplyr::pull(frac)

      dplyr::tibble(class = c("HM", "Classic"), frac = c(fHM, fClassic), gene = GENE,
                    ttype = tt,
                    p_value = pval)
    }
  }) %>% do.call("bind_rows", .)

  gene_res
}) %>% do.call("bind_rows", .)

p = gene_analysis_df %>%
  tidyr::pivot_wider(names_from = class, values_from = frac) %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(p_value = p.adjust(p_value, method="BH")) %>%
  dplyr::filter(p_value <= .05) %>%
  tidyr::pivot_longer(!c(gene, ttype, p_value)) %>%
  ggplot(mapping = aes(x = gene, y=value, fill=name)) +
  geom_col(position = "dodge") +
  facet_grid(~ttype, space = "free_x", scales = "free") +
  theme_bw() +
  scale_fill_manual(values = class_colors) +
  labs(fill="", x='Gene', y="Incidence fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5))
p
ggsave("plot/gene_incidence_fraction.pdf", width = 14, height = 5, dpi = 600, units = "in", plot = p)


p = gene_analysis_df %>%
  tidyr::pivot_wider(names_from = class, values_from = frac) %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(p_value = p.adjust(p_value, method="BH")) %>%
  dplyr::filter(p_value <= .05) %>%
  tidyr::pivot_longer(!c(gene, ttype, p_value)) %>%
  dplyr::filter(ttype %in% c("ESAD")) %>%
  ggplot(mapping = aes(x = gene, y=value, fill=name)) +
  geom_col(position = "dodge") +
  facet_grid(~ttype, space = "free_x", scales = "free") +
  theme_bw() +
  scale_fill_manual(values = class_colors) +
  labs(fill="", x='Gene', y="Incidence fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  theme(legend.position = "right")
p
saveRDS(p, "plot/GACA_ESAD_gene_incidence_fraction.rds")
ggsave("plot/GACA_ESAD_gene_incidence_fraction.pdf", width = 14, height = 5, dpi = 600, units = "in", plot = p)



gene_analysis_df %>%
  tidyr::pivot_wider(names_from = class, values_from = frac) %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(p_value = p.adjust(p_value, method="BH")) %>%
  dplyr::filter(p_value <= .05) %>%
  tidyr::pivot_longer(!c(gene, ttype, p_value)) %>%
  dplyr::filter(ttype == "BRCA") %>%
  dplyr::filter(gene %in% c("BRCA1", "BRCA2"))
