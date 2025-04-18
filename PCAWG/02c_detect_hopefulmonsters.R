
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("utils.R")

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

N_CHR = 11

RES = lapply(unique(RES$sample_id), function(s) {
  RES %>%
    dplyr::filter(sample_id == s) %>%
    dplyr::select(sample_id, ttype, ploidy, ttype, karyotype, clock_rank, wgd_status, chr) %>%
    dplyr::mutate(n_cna = n(), n_clusters = max(clock_rank)) %>%
    dplyr::group_by(clock_rank) %>%
    dplyr::mutate(n_chr_affected = length(unique(chr))) %>%
    dplyr::select(sample_id, ttype, ploidy, wgd_status, n_cna, n_clusters, n_chr_affected, clock_rank) %>%
    dplyr::distinct() %>%
    dplyr::mutate(is_HM = ifelse((wgd_status == "no_wgd" | ploidy == 2) & any(n_chr_affected >= N_CHR), "HM", ifelse(wgd_status == "wgd", "WGD", "Classic")))
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

p = RES %>%
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
saveRDS(p, "plot/HM_scatter_with_smooth.rds")
ggsave("plot/HM_scatter_with_smooth.pdf", width = 8, height = 8, units = "in", plot = p)

p = RES %>%
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
saveRDS(p, "plot/WGD_HM_scatter_with_smooth.rds")
ggsave("plot/WGD_HM_scatter_with_smooth.pdf", width = 8, height = 8, units = "in", plot = p)


p = RES %>%
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
saveRDS(p, "plot/all_scatter_Ncna_v_Nclusters_with_HM.rds")
ggsave("plot/all_scatter_Ncna_v_Nclusters_with_HM.pdf", width = 10, height = 10, units = "in", plot = p)

RES %>%
  dplyr::group_by(ttype, is_HM) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::pull(n) %>% sum()

RES %>%
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


df_incidence = RES %>%
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
df_incidence$ttype = factor(df_incidence$ttype, levels=df_incidence$ttype)

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
saveRDS(p, "plot/distribution_of_HM_fraction_per_ttype.rds")
ggsave("plot/distribution_of_HM_fraction_per_ttype.pdf", width = 8, height = 8, units = "in", plot = p)

# p = RES %>%
#   dplyr::group_by(ttype, is_HM) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(ttype) %>%
#   dplyr::mutate(sn = sum(n)) %>%
#   dplyr::filter(is_HM %in% c("HM")) %>%
#   dplyr::mutate(f = n / sn) %>%
#   dplyr::arrange(-f) %>%
#   ggplot(mapping = aes(x=reorder(ttype, f), y=f, fill=is_HM)) +
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_manual(values = class_colors) +
#   theme_bw() +
#   coord_flip() +
#   #scale_y_continuous(limits = c(0,.2)) +
#   labs(x = "Tumour type", y="HM fraction") +
#   theme(legend.position = "none")
# p
# saveRDS(p, "plot/distribution_of_HM_fraction_per_ttype.rds")
# ggsave("plot/distribution_of_HM_fraction_per_ttype.pdf", width = 8, height = 8, units = "in", plot = p)

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

saveRDS(p, "plot/distribution_of_HM_and_WGD_fraction_per_ttype.rds")
ggsave("plot/distribution_of_HM_and_WGD_fraction_per_ttype.pdf", width = 8, height = 8, units = "in", plot = p)

# p = RES %>%
#   dplyr::group_by(sample_id) %>%
#   dplyr::group_by(ttype, sample_id) %>%
#   dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>%
#   dplyr::mutate(n_chr = length(unique(chr))) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(sample_id, ttype, n_events, n_clusters, n_chr) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(class = ifelse(n_chr > N_CHR, "HM", "Classic")) %>%
#   ggplot(mapping = aes(x=n_events, y=n_clusters, fill=class, col=class)) +
#   stat_density_2d(geom = "polygon", aes(alpha = ..level..)) +
#   geom_point() +
#   scale_x_continuous(transform = "log10") +
#   scale_y_continuous(transform = "log10") +
#   scale_color_manual(values = class_colors) +
#   scale_fill_manual(values = class_colors) +
#   theme_bw()
# p
#
# saveRDS(p, "plot/density_HM_dist.rds")
# ggsave("plot/density_HM_dist.pdf", width = 8, height = 8, units = "in", plot = p)

# q9 = data %>%
#   dplyr::mutate(x = n_events / n_clusters) %>%
#   dplyr::pull(x) %>%
#   stats::quantile(.8)
#
# rr = RES %>%
#   dplyr::group_by(ttype, sample_id) %>%
#   dplyr::mutate(n_events=n(), n_clusters=length(unique(clock_mean))) %>%
#   dplyr::select(sample_id, n_events, n_clusters) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(f = n_events / n_clusters) %>%
#   dplyr::mutate(class = ifelse(f >= q9, "Hopeful Monster", "Classic"))
#
# rr %>%
#   ggplot(mapping = aes(x=n_events, y=n_clusters, col=class)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~ttype) +
#   theme(legend.position = "bottom") +
#   labs(x = "N events", y="N clusters", col="") +
#   scale_color_manual(values = class_colors)
#
# rr %>%
#   dplyr::group_by(ttype, class) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(ttype) %>%
#   dplyr::mutate(n = n / sum(n)) %>%
#   ggplot(mapping = aes(x=reorder(ttype, -n), y=n, fill=class)) +
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_manual(values = class_colors) +
#   theme_bw() +
#   coord_flip()
#
