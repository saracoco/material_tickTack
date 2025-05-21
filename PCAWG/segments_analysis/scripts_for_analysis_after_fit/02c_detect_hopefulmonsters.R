
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide INPUT_DIR as command line argument.")
}

INPUT_DIR <- args[1]
RDS_OUTPUT = file.path(INPUT_DIR, "rds")
PLOT_OUTPUT = file.path(INPUT_DIR, "pdfs")

dir.create(RDS_OUTPUT, recursive = TRUE)
dir.create(PLOT_OUTPUT, recursive = TRUE)
#INPUT_DIR = "tmp_output/"

class_colors = list(
  'Classic' = "#f1a340",
  'HM' = "#998ec3",
  "WGD" = "#a6dba0"
)

k_colors = list(
  '2:0' = 'turquoise4',
  '2:1' = ggplot2::alpha('orange', .8),
  '2:2' = 'firebrick3'
)

load("data/gene_coordinates_hg19.rda")

ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>%
  dplyr::select(sample_id, ttype) %>%
  dplyr::distinct()

META = readRDS("data/metadata_all.rds") %>% dplyr::select(sample_id, wgd_status)
RES = readRDS(file.path(INPUT_DIR, "summary_all_samples.rds"))
RES$sample_id = lapply(RES$sample_id, function(s) {
  str_replace_all(s, ".rds", "")
}) %>% unlist()
RES = RES %>%
  dplyr::left_join(META, by="sample_id")

#N_CHR = 11

# Add fraction of "Modified genome" per clock
total_bps_in_genome = sum(gene_coordinates_hg19$to - gene_coordinates_hg19$from + 1)
RES = RES %>% 
  dplyr::group_by(sample_id, clock_rank) %>% 
  dplyr::mutate(frac_genome_affected = sum(to - from + 1) / total_bps_in_genome) %>% 
  dplyr::ungroup()


RES = lapply(unique(RES$sample_id), function(s) {
  print(s)
  RES %>%
    dplyr::filter(sample_id == s) %>%
    dplyr::select(sample_id, ttype, ploidy, ttype, karyotype, clock_rank, clock_mean, wgd_status, chr, frac_genome_affected) %>%
    dplyr::mutate(n_cna = n(), n_clusters = max(clock_rank)) %>%
    dplyr::group_by(clock_rank) %>%
    dplyr::mutate(n_chr_affected = length(unique(chr))) %>%
    dplyr::select(sample_id, ttype, ploidy, wgd_status, n_cna, n_clusters, n_chr_affected, clock_rank, clock_mean, frac_genome_affected) %>%
    dplyr::distinct() #%>% 
    #dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd" | ploidy == 2) & any(n_chr_affected >= N_CHR), "HM", ifelse(wgd_status == "wgd", "WGD", "Classic")))
}) %>% do.call("bind_rows", .)


# Definition by chromosome ####
N_CHR = 11

res = RES %>%
  dplyr::group_by(sample_id, clock_rank) %>%
  dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd" | ploidy == 2) & any(n_chr_affected >= N_CHR), "HM", ifelse(wgd_status == "wgd", "WGD", "Classic"))) %>%
  #dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd") & any(n_chr_affected >= N_CHR), "HM", "Classic")) %>%
  na.omit()

res = res %>%
  dplyr::select(!c(n_chr_affected, clock_rank)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct()

p = res %>%
  #dplyr::filter(wgd_status == "no_wgd") %>%
  dplyr::filter(is_HM %in% c("Classic", "HM")) %>%
  ggplot(mapping = aes(x=n_cna, y=n_clusters, col=is_HM, label=is_HM)) +
  geom_point(alpha = .8) +
  geomtextpath::geom_labelsmooth(fill="white", method = "lm", formula = y~x) +
  scale_color_manual(values = class_colors) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "N CNA", y= "N clusters") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_x_continuous(transform = "log10")
p

saveRDS(p, file.path(RDS_OUTPUT, "nchr_HM_scatter_with_smooth.rds"))
ggsave(file.path(PLOT_OUTPUT, "nchr_HM_scatter_with_smooth.pdf"), width = 8, height = 8, units = "in", plot = p)
#saveRDS(p, "plot/HM_scatter_with_smooth.rds")
#ggsave("plot/HM_scatter_with_smooth.pdf", width = 8, height = 8, units = "in", plot = p)

p = res %>%
  #dplyr::filter(wgd_status == "no_wgd") %>%
  ggplot(mapping = aes(x=n_cna, y=n_clusters, col=is_HM, label=is_HM)) +
  geom_point(alpha = .8) +
  geomtextpath::geom_labelsmooth(fill="white", method = "lm", formula = y~x) +
  scale_color_manual(values = class_colors) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "N CNA", y= "N clusters") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_x_continuous(transform = "log10")
p

saveRDS(p, file.path(RDS_OUTPUT, "nchr_WGD_HM_scatter_with_smooth.rds"))
ggsave(file.path(PLOT_OUTPUT, "nchr_WGD_HM_scatter_with_smooth.pdf"), width = 8, height = 8, units = "in", plot = p)
# saveRDS(p, "plot/WGD_HM_scatter_with_smooth.rds")
# ggsave("plot/WGD_HM_scatter_with_smooth.pdf", width = 8, height = 8, units = "in", plot = p)

p = res %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(nsamples = n()) %>%
  dplyr::filter(is_HM %in% c("HM", "Classic")) %>%
  ggplot(mapping = aes(x=n_cna, y=n_clusters, col=is_HM)) +
  geom_point() +
  facet_wrap(~ttype) +
  scale_color_manual(values = class_colors) +
  theme_bw() +
  labs(x = "N CNA", y= "N clusters", col="") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  theme(legend.position = "bottom")
p
saveRDS(p, file.path(RDS_OUTPUT, "nchr_all_scatter_Ncna_v_Nclusters_with_HM.rds"))
ggsave(file.path(PLOT_OUTPUT, "nchr_all_scatter_Ncna_v_Nclusters_with_HM.pdf"), width = 8, height = 8, units = "in", plot = p)
#saveRDS(p, "plot/all_scatter_Ncna_v_Nclusters_with_HM.rds")
#ggsave("plot/all_scatter_Ncna_v_Nclusters_with_HM.pdf", width = 10, height = 10, units = "in", plot = p)

res %>%
  dplyr::group_by(ttype, is_HM) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::pull(n) %>% sum()

res %>%
  dplyr::group_by(ttype, is_HM) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(sn = sum(n)) %>%
  #dplyr::filter(is_HM %in% c("HM", "WGD")) %>%
  dplyr::mutate(f = n / sn) %>%
  dplyr::select(f, is_HM, ttype) %>%
  ggplot(mapping = aes(x = ttype, y=f, fill=is_HM)) +
  geom_col(position = "stack")

df_incidence = res %>%
  dplyr::left_join(readRDS("data/metadata_all.rds") %>%
                     dplyr::group_by(tumor_type) %>%
                     dplyr::summarise(nPCAWG=n()) %>%
                     dplyr::rename(ttype=tumor_type), by="ttype") %>%
  dplyr::group_by(ttype, is_HM, nPCAWG) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(sn = sum(n)) %>%
  dplyr::filter(is_HM %in% c("HM", "WGD")) %>%
  dplyr::mutate(f = n / sn) %>%
  dplyr::mutate(fPCAWG = n / nPCAWG) %>%
  dplyr::select(fPCAWG, f, is_HM, ttype)

#df_incidence$ttype = factor(df_incidence$ttype, levels=df_incidence$ttype)

p = df_incidence %>%
  dplyr::select(ttype, fPCAWG, is_HM) %>%
  dplyr::filter(is_HM == "HM") %>%
  ggplot(mapping = aes(x=reorder(ttype, fPCAWG), y=fPCAWG, fill=is_HM)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  coord_flip() +
  #scale_y_continuous(limits = c(0,.2)) +
  labs(x = "Tumour type", y="HM fraction") +
  theme(legend.position = "none")
p
saveRDS(p, file.path(RDS_OUTPUT, "nchr_distribution_of_HM_fraction_per_ttype.rds"))
ggsave(file.path(PLOT_OUTPUT, "nchr_distribution_of_HM_fraction_per_ttype.pdf"), width = 8, height = 8, units = "in", plot = p)
#saveRDS(p, "plot/distribution_of_HM_fraction_per_ttype.rds")
#ggsave("plot/distribution_of_HM_fraction_per_ttype.pdf", width = 8, height = 8, units = "in", plot = p)

p = df_incidence %>%
  dplyr::select(ttype, fPCAWG, is_HM) %>%
  ggplot(mapping = aes(x=reorder(ttype, fPCAWG), y=fPCAWG, fill=is_HM)) +
  geom_col(position="dodge") +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  coord_flip() +
  #scale_y_continuous(limits = c(0,.2)) +
  labs(x = "Tumour type", y="Incidence fraction") +
  theme(legend.position = "none")
p

saveRDS(p, file.path(RDS_OUTPUT, "nchr_distribution_of_HM_and_WGD_fraction_per_ttype.rds"))
ggsave(file.path(PLOT_OUTPUT, "nchr_distribution_of_HM_and_WGD_fraction_per_ttype.pdf"), width = 8, height = 8, units = "in", plot = p)

p = res %>%
  # Remove samples with WGD
  dplyr::filter(is_HM != "WGD") %>% 
  # Group by sample_id and keep only those that have both Classic and HM
  group_by(sample_id) %>%
  filter(all(c("Classic", "HM") %in% is_HM)) %>%
  # For samples with multiple HM rows, keep only the one with highest frac_genome_affected
  mutate(rank = if_else(is_HM == "HM", rank(-frac_genome_affected, ties.method = "first"), NA_real_)) %>%
  filter(is_HM != "HM" | rank == 1) %>%
  select(-rank) %>%
  ungroup() %>% 
  ggplot(mapping = aes(x=ttype, col=is_HM, y=clock_mean)) +
  geom_boxplot() +
  #geom_jitter() +
  theme_bw() +
  labs(x = "Tumour type", y = "Pseudo-time", col = "")

saveRDS(p, file.path(RDS_OUTPUT, "nchr_psuedotimes_HM_vs_Classical.rds"))
ggsave(file.path(PLOT_OUTPUT, "nchr_psuedotimes_HM_vs_Classical.pdf"), width = 8, height = 8, units = "in", plot = p)

#saveRDS(p, "plot/distribution_of_HM_and_WGD_fraction_per_ttype.rds")
#ggsave("plot/distribution_of_HM_and_WGD_fraction_per_ttype.pdf", width = 8, height = 8, units = "in", plot = p)

# Definition by Fraction ####
MIN_FRAC = .25

res = RES %>%
  dplyr::group_by(sample_id, clock_rank) %>%
  dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd" | ploidy == 2) & (frac_genome_affected >= MIN_FRAC), "HM", ifelse(wgd_status == "wgd", "WGD", "Classic"))) %>%
  #dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd") & any(n_chr_affected >= N_CHR), "HM", "Classic")) %>%
  na.omit()

res = res %>%
  dplyr::select(!c(n_chr_affected, clock_rank)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct()

p = res %>%
  #dplyr::filter(wgd_status == "no_wgd") %>%
  dplyr::filter(is_HM %in% c("Classic", "HM")) %>%
  ggplot(mapping = aes(x=n_cna, y=n_clusters, col=is_HM, label=is_HM)) +
  geom_point(alpha = .8) +
  geomtextpath::geom_labelsmooth(fill="white", method = "lm", formula = y~x) +
  scale_color_manual(values = class_colors) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "N CNA", y= "N clusters") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_x_continuous(transform = "log10")
p

saveRDS(p, file.path(RDS_OUTPUT, "frac_HM_scatter_with_smooth.rds"))
ggsave(file.path(PLOT_OUTPUT, "frac_HM_scatter_with_smooth.pdf"), width = 8, height = 8, units = "in", plot = p)
#saveRDS(p, "plot/HM_scatter_with_smooth.rds")
#ggsave("plot/HM_scatter_with_smooth.pdf", width = 8, height = 8, units = "in", plot = p)

p = res %>%
  #dplyr::filter(wgd_status == "no_wgd") %>%
  ggplot(mapping = aes(x=n_cna, y=n_clusters, col=is_HM, label=is_HM)) +
  geom_point(alpha = .8) +
  geomtextpath::geom_labelsmooth(fill="white", method = "lm", formula = y~x) +
  scale_color_manual(values = class_colors) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "N CNA", y= "N clusters") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_x_continuous(transform = "log10")
p

saveRDS(p, file.path(RDS_OUTPUT, "frac_WGD_HM_scatter_with_smooth.rds"))
ggsave(file.path(PLOT_OUTPUT, "frac_WGD_HM_scatter_with_smooth.pdf"), width = 8, height = 8, units = "in", plot = p)
# saveRDS(p, "plot/WGD_HM_scatter_with_smooth.rds")
# ggsave("plot/WGD_HM_scatter_with_smooth.pdf", width = 8, height = 8, units = "in", plot = p)

p = res %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(nsamples = n()) %>%
  dplyr::filter(is_HM %in% c("HM", "Classic")) %>%
  ggplot(mapping = aes(x=n_cna, y=n_clusters, col=is_HM)) +
  geom_point() +
  facet_wrap(~ttype) +
  scale_color_manual(values = class_colors) +
  theme_bw() +
  labs(x = "N CNA", y= "N clusters", col="") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  theme(legend.position = "bottom")
p
saveRDS(p, file.path(RDS_OUTPUT, "frac_all_scatter_Ncna_v_Nclusters_with_HM.rds"))
ggsave(file.path(PLOT_OUTPUT, "frac_all_scatter_Ncna_v_Nclusters_with_HM.pdf"), width = 8, height = 8, units = "in", plot = p)
#saveRDS(p, "plot/all_scatter_Ncna_v_Nclusters_with_HM.rds")
#ggsave("plot/all_scatter_Ncna_v_Nclusters_with_HM.pdf", width = 10, height = 10, units = "in", plot = p)

res %>%
  dplyr::group_by(ttype, is_HM) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::pull(n) %>% sum()

res %>%
  dplyr::group_by(ttype, is_HM) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(sn = sum(n)) %>%
  #dplyr::filter(is_HM %in% c("HM", "WGD")) %>%
  dplyr::mutate(f = n / sn) %>%
  dplyr::select(f, is_HM, ttype) %>%
  ggplot(mapping = aes(x = ttype, y=f, fill=is_HM)) +
  geom_col(position = "stack")

df_incidence = res %>%
  dplyr::left_join(readRDS("data/metadata_all.rds") %>%
                     dplyr::group_by(tumor_type) %>%
                     dplyr::summarise(nPCAWG=n()) %>%
                     dplyr::rename(ttype=tumor_type), by="ttype") %>%
  dplyr::group_by(ttype, is_HM, nPCAWG) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(ttype) %>%
  dplyr::mutate(sn = sum(n)) %>%
  dplyr::filter(is_HM %in% c("HM", "WGD")) %>%
  dplyr::mutate(f = n / sn) %>%
  dplyr::mutate(fPCAWG = n / nPCAWG) %>%
  dplyr::select(fPCAWG, f, is_HM, ttype)

#df_incidence$ttype = factor(df_incidence$ttype, levels=df_incidence$ttype)

p = df_incidence %>%
  dplyr::select(ttype, fPCAWG, is_HM) %>%
  dplyr::filter(is_HM == "HM") %>%
  ggplot(mapping = aes(x=reorder(ttype, fPCAWG), y=fPCAWG, fill=is_HM)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  coord_flip() +
  #scale_y_continuous(limits = c(0,.2)) +
  labs(x = "Tumour type", y="HM fraction") +
  theme(legend.position = "none")
p
saveRDS(p, file.path(RDS_OUTPUT, "frac_distribution_of_HM_fraction_per_ttype.rds"))
ggsave(file.path(PLOT_OUTPUT, "frac_distribution_of_HM_fraction_per_ttype.pdf"), width = 8, height = 8, units = "in", plot = p)
#saveRDS(p, "plot/distribution_of_HM_fraction_per_ttype.rds")
#ggsave("plot/distribution_of_HM_fraction_per_ttype.pdf", width = 8, height = 8, units = "in", plot = p)

p = df_incidence %>%
  dplyr::select(ttype, fPCAWG, is_HM) %>%
  ggplot(mapping = aes(x=reorder(ttype, fPCAWG), y=fPCAWG, fill=is_HM)) +
  geom_col(position="dodge") +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  coord_flip() +
  #scale_y_continuous(limits = c(0,.2)) +
  labs(x = "Tumour type", y="Incidence fraction") +
  theme(legend.position = "none")
p

saveRDS(p, file.path(RDS_OUTPUT, "frac_distribution_of_HM_and_WGD_fraction_per_ttype.rds"))
ggsave(file.path(PLOT_OUTPUT, "frac_distribution_of_HM_and_WGD_fraction_per_ttype.pdf"), width = 8, height = 8, units = "in", plot = p)

p = res %>%
  # Remove samples with WGD
  dplyr::filter(is_HM != "WGD") %>% 
  # Group by sample_id and keep only those that have both Classic and HM
  group_by(sample_id) %>%
  filter(all(c("Classic", "HM") %in% is_HM)) %>%
  # For samples with multiple HM rows, keep only the one with highest frac_genome_affected
  mutate(rank = if_else(is_HM == "HM", rank(-frac_genome_affected, ties.method = "first"), NA_real_)) %>%
  filter(is_HM != "HM" | rank == 1) %>%
  select(-rank) %>%
  ungroup() %>% 
  ggplot(mapping = aes(x=ttype, col=is_HM, y=clock_mean)) +
  geom_boxplot() +
  #geom_jitter() +
  theme_bw() +
  labs(x = "Tumour type", y = "Pseudo-time", col = "")
p

saveRDS(p, file.path(RDS_OUTPUT, "frac_psuedotimes_HM_vs_Classical.rds"))
ggsave(file.path(PLOT_OUTPUT, "frac_psuedotimes_HM_vs_Classical.pdf"), width = 8, height = 8, units = "in", plot = p)
