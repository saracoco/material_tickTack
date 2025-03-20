rm(list=ls())
.libPaths(new="~/R/rstudio_v3/")
# source("./expectations_subclonal.R")
# source("./simulate_complex.R")
source("./simulate_functions.R")
x <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/data/clonal_analysis_PCAWG/0009b464-b376-4fbc-8a56-da538269a02f/fit.rds")
x$mutations <- x$snvs


# understand from subclonal analysis how to generate evolving data for the clonal case 
min_VAF = 0
n_min = 10
x$mutations <- x$mutations %>% mutate(VAF = NV/DP)
candidates = x$mutations%>% dplyr::filter(VAF > min_VAF) %>% dplyr::pull(karyotype) %>% unique
candidates = setdiff(candidates,  c("1:1", "1:0", "2:0", "2:1", "2:2", "NA:NA"))
n_cand = x$n_karyotype[candidates] >= n_min

analysis = names(n_cand)[n_cand]
cna = x$cna %>% filter(n > n_min)

#expectations_subclonal
starting_state = '1:1'
i = 100
CCF_1 = cna$CCF[i]
# Get karyotype
Major <- cna$Major[i]
minor <- cna$minor[i]
purity = x$purity
Major <- cna$Major[i]
minor <- cna$minor[i]
karyotype <- paste(Major, minor, sep=':')
karyotype_1 = karyotype
karyotype_2 = karyotype


#branching evolution
#input
starting = starting_state


left = karyotype
right = karyotype
target=starting



b_model = suppressWarnings(branching_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
  Reduce(f = bind_rows)


l1_model = NULL
if((grepl("0", karyotype_1) & grepl("0", karyotype_2)) | !grepl("0", karyotype_1))
  l1_model = suppressWarnings(linear_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
  Reduce(f = bind_rows)




# expected_peaks = lapply(analysis,
#                         function(k)
#                           expectations_generalised(
#                             m = NULL,
#                             M = NULL,
#                             p = x$purity,
#                             karyotype = k
#                           )) %>%
#   Reduce(f = bind_rows)












