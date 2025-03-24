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




first_child = karyotype_1
second_child = karyotype_2



l1_model = NULL
if((grepl("0", karyotype_1) & grepl("0", karyotype_2)) | !grepl("0", karyotype_1))
  l1_model = suppressWarnings(linear_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
  Reduce(f = bind_rows)

# error in evolve 
#inputs
copy_state = initial_state(starting)
target = first_child



cap = (target %>% strsplit(split = ':'))[[1]] %>% as.numeric %>% sum
cap = 2 * cap
copy_state_tmp = copy_state
copy_state = copy_state_tmp
current_state = original_state = list(copy_state)
current_state_tmp = current_state
copy_state = current_state


mutations = copy_state[[1]] %>% pull(mutations) %>% unlist()
mutations = copy_state[[1]] %>% filter(allele == !!"A1") %>% pull(mutations) %>% unlist()



get_mutations = function(copy_state, allele = NULL){
  if(is.null(allele))
    mutations = copy_state %>% pull(mutations) %>% unlist()
  else
    mutations = copy_state %>% filter(allele == !!allele) %>% pull(mutations) %>% unlist()
  
  return(mutations)
}


# Amplify allele
amplify =  function(copy_state){
  augment = function(which_allele)
  {
    mutations_allele = copy_state[[1]] %>% get_mutations(allele = which_allele)
    
    maj_min = substr(which_allele, 0, 1)
    n_allele = gsub(x = which_allele, 'A', "") %>% gsub(pattern = 'B', replacement = "") %>%
      as.numeric()
    
    # New allele can be +1, unless there are other alleles with larger number
    new_allele = paste(maj_min, n_allele + 1, sep = '')
    
    while(new_allele %in% (copy_state %>% get_alleles())) {
      n_allele = n_allele + 1
      new_allele = paste(maj_min, n_allele + 1, sep = '')
    }
    
    new_entry = copy_state[[1]] %>% filter(allele == which_allele)
    new_entry$allele = new_allele
    
    copy_state %>% bind_rows(new_entry) %>% dplyr::arrange(allele)
  }
  
  copy_state %>%
    get_alleles() %>%
    lapply(augment)
}

lapply(current_state, amplify)
del_new_state = lapply(current_state, delete)

# expected_peaks = lapply(analysis,
#                         function(k)
#                           expectations_generalised(
#                             m = NULL,
#                             M = NULL,
#                             p = x$purity,
#                             karyotype = k
#                           )) %>%
#   Reduce(f = bind_rows)












