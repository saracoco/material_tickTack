.libPaths(new="~/R/rstudio_v3/") 
rm(list=ls())
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
require(tidyverse)
source("../utils.R")


# generate a new RES file including complex copy number
ttypes <- read.delim("../data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

IDs <- list.files("../../../data/clonal_analysis_PCAWG/")
IDs = IDs[!grepl("single", IDs)]

fits_path = "../../../data/clonal_analysis_PCAWG/"





parse_summarized_results <- function(fit) {
  fit %>% 
    tidyr::separate(segment_name, into = c("chr", "from", "to"), sep = "_", convert = TRUE) %>% 
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>% 
    dplyr::select(chr, from, to, segment_id, karyotype, clock_mean, clock_low, clock_high, clock_rank)
}


id = IDs[1]
RES <- lapply(IDs, function(id) {
  print(which(IDs == id) / length(IDs) * 100)
  
  tryCatch({
    fit = readRDS(paste0(fits_path, unlist(strsplit(id, ".rds")), "/fit.rds"))
    tumour_name = (unique(fit$snvs$ttype) %>% na.omit())[1]
    ploidy = fit$ploidy
    ttype = strsplit(fit$snvs$project_code, "-")[[1]][1]

    single_data <- fit$snvs%>% 
      dplyr::select(chr, from, to, segment_id, karyotype)
    
    df = dplyr::tibble(sample_id=id, ttype=ttype, ploidy=ploidy, tumour_name=tumour_name)
    return(dplyr::bind_cols(df, single_data))
  }, error = function(e) {
    # Error handling
    print(paste0("An error occurred:", id))
    return(NULL)
  })
  
}) %>% do.call("bind_rows", .)

saveRDS(RES, "./RES_karyotypes.rds")














k_colors = list(
  '2:0' = 'turquoise4',
  '2:1' = ggplot2::alpha('orange', .8),
  '2:2' = 'firebrick3'  
)

ttypes <- read.delim("../data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  dplyr::select(sample_id, ttype) %>% 
  dplyr::distinct()

RES = readRDS("../results/summary_all_samples.rds")


# data_path <-  "../../data/clonal_analysis_PCAWG"
# sample_files = list.files(data_path, full.names = "T")



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