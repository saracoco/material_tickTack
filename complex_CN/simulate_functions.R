expected_vaf_fun = function(m, M, mut.allele, p)
{
  P = m + M
  
  expected_mutant_reads = mut.allele * p
  expected_sequencing_depth = 2 * (1 - p) + p * P
  expected_mutant_reads / expected_sequencing_depth
}

# Expected VAF for a general peak, using expected_vaf_fun
expectations_generalised = function(m, M, p, karyotype = NULL)
{
  if(!is.null(karyotype))
  {
    karyotype = strsplit(karyotype, ':')[[1]]
    m = karyotype[2] %>% as.numeric()
    M = karyotype[1] %>% as.numeric()
  }
  
  minor_peaks = lapply(1:m, function(i){
    data.frame(
      minor = m,
      Major = M,
      ploidy = m + M,
      multiplicity = i,
      purity = p,
      peak = expected_vaf_fun(m, M, mut.allele = i, p)
    )
  })
  
  Major_peaks = lapply(1:M, function(i){
    data.frame(
      minor = m,
      Major = M,
      ploidy = m + M,
      multiplicity = i,
      purity = p,
      peak = expected_vaf_fun(m, M, mut.allele = i, p)
    )
  })
  
  return(
    bind_rows(
      Reduce(bind_rows, minor_peaks),
      Reduce(bind_rows, Major_peaks)
    ) %>%
      distinct() %>%
      mutate(karyotype = paste(Major, minor, sep = ':')) %>%
      filter(multiplicity > 0)
  )
}

# exampe use :
# 
# x$mutations <- x$mutations %>% mutate(VAF = NV/DP)
# candidates = x$mutations%>% dplyr::filter(VAF > min_VAF) %>% dplyr::pull(karyotype) %>% unique
# candidates = setdiff(candidates,  c("1:1", "1:0", "2:0", "2:1", "2:2", "NA:NA"))
# n_cand = x$n_karyotype[candidates] >= n_min
# 
# analysis = names(n_cand)[n_cand]
#
# expected_peaks = lapply(analysis,
#                         function(k)
#                           expectations_generalised(
#                             m = NULL,
#                             M = NULL,
#                             p = x$purity,
#                             karyotype = k
#                           )) %>%
#   Reduce(f = bind_rows)










###################### functions for branching evolution : initial state #####################

# Add mutations to all alleles
mutation = function(copy_state){
  for(i in 1:nrow(copy_state))
  {
    muts = copy_state$mutations[[i]] %>% unlist()
    muts = list(c(muts, copy_state %>% new_mutation))
    copy_state$mutations[[i]] = muts
  }
  copy_state
}

# generate a new mutation ID
new_mutation = function(copy_state){
  used = copy_state %>% get_mutations()
  id = NULL
  repeat{
    id = sample(LETTERS, 8, replace = TRUE) %>% paste(collapse = "")
    if(!id %in% used ) break
  }
  
  return(id)
}

# getters
get_mutations = function(copy_state, allele = NULL){
  if(is.null(allele))
    mutations = copy_state %>% pull(mutations) %>% unlist()
  else
    mutations = copy_state %>% filter(allele == !!allele) %>% pull(mutations) %>% unlist()
  
  return(mutations)
}


# initial state
initial_state = function(target)
{
  copy_state = data.frame(allele = c("A1", "B1")) %>% as_tibble()
  copy_state$mutations = NULL
  copy_state$mutations[[1]] = list()
  copy_state$mutations[[2]] = list()
  
  copy_state = copy_state %>% mutation()
  
  if(target == "1:1") return(copy_state)
  else evolve(copy_state, target)[[1]]
}

###########################################################################################

###################### functions for branching evolution : evolve #####################

# Get stast
as_karyotype = function(copy_state){
  tab_counts = copy_state %>% get_alleles() %>% substr(0, 1) %>% table() %>% sort(decreasing = TRUE)
  state = as.numeric(tab_counts) %>% paste(collapse = ':')
  if(!grepl(":", state)) state = paste0(state, ':0')
  
  return(state)
}


get_alleles = function(copy_state){
  copy_state$allele
}




# Amplify allele
amplify =  function(copy_state){
  augment = function(which_allele)
  {
    mutations_allele = copy_state %>% get_mutations(allele = which_allele)
    
    maj_min = substr(which_allele, 0, 1)
    n_allele = gsub(x = which_allele, 'A', "") %>% gsub(pattern = 'B', replacement = "") %>%
      as.numeric()
    
    # New allele can be +1, unless there are other alleles with larger number
    new_allele = paste(maj_min, n_allele + 1, sep = '')
    
    while(new_allele %in% (copy_state %>% get_alleles())) {
      n_allele = n_allele + 1
      new_allele = paste(maj_min, n_allele + 1, sep = '')
    }
    
    new_entry = copy_state %>% filter(allele == which_allele)
    new_entry$allele = new_allele
    
    copy_state %>% bind_rows(new_entry) %>% dplyr::arrange(allele)
  }
  
  copy_state %>%
    get_alleles() %>%
    lapply(augment)
}


# Delete allele
delete =  function(copy_state){
  cancel = function(which_allele)
  {
    mutations_allele = copy_state %>% get_mutations(allele = which_allele)
    
    copy_state %>% filter(allele != which_allele) %>% arrange(allele)
  }
  
  copy_state %>%
    get_alleles() %>%
    lapply(cancel)
}




# genome_double
genome_double =  function(copy_state){
  copy_of = copy_state
  
  for(i in 1:nrow(copy_of))
  {
    ni_allele = substr(copy_state$allele[i], 2, nchar(copy_state$allele[i])) %>% as.numeric
    ci_allele = substr(copy_state$allele[i], 0, 1)
    
    # New allele can be +1, unless there are other alleles with larger number
    new_t = paste0(ci_allele, ni_allele + 1)
    
    while(new_t %in% (copy_state %>% bind_rows(copy_of) %>%  get_alleles()) %>% unique) {
      ni_allele = ni_allele + 1
      new_t = paste0(ci_allele, ni_allele + 1)
    }
    
    copy_of$allele[i] = new_t
    # copy_of$allele[i] = paste0(ci_allele, ni_allele + 1)
  }
  
  copy_state %>% bind_rows(copy_of) %>% arrange(allele) %>%  list()
}


as_ploidy = function(copy_state){
  copy_state %>% nrow()
}


# Compute peaks
get_peaks = function(clone_1, clone_2, CCF_1, purity)
{
  c(
    clone_1 %>% get_mutations(),
    clone_2 %>% get_mutations()
  ) %>%
    table()
  
  m_c1 = clone_1 %>% get_mutations() %>% table() %>% as_tibble() %>%
    mutate(x = n * CCF_1,
           karyotype_1 = clone_1 %>% as_karyotype(),
           genotype_1 = clone_1 %>% get_alleles() %>% sort() %>% paste(collapse = '')
    )
  colnames(m_c1)[1] = 'mutation'
  
  m_c2 = clone_2 %>% get_mutations() %>% table() %>% as_tibble() %>%
    mutate(x = n * (1 - CCF_1),
           karyotype_2 = clone_2 %>% as_karyotype(),
           genotype_2 = clone_2 %>% get_alleles() %>% sort() %>% paste(collapse = '')
    )
  colnames(m_c2)[1] = 'mutation'
  
  denominator = 2 * (1-purity) +
    purity * ( CCF_1 * (clone_1 %>% as_ploidy()) + (1 - CCF_1) * (clone_2 %>% as_ploidy()))
  
  m_c1 %>%
    dplyr::full_join(m_c2, by = 'mutation', suffix = c('.clone_1', '.clone_2')) %>%
    tidyr::replace_na(list(n.clone_1 = 0, x.clone_1 = 0, n.clone_2 = 0, x.clone_2 = 0)) %>%
    dplyr::mutate(peak = (x.clone_1 + x.clone_2) * purity) %>%
    dplyr::distinct(peak, .keep_all = TRUE) %>%
    dplyr::mutate(peak = peak/denominator) %>%
    dplyr::select(mutation, karyotype_1, genotype_1, karyotype_2, genotype_2, n.clone_1, n.clone_2, peak) %>%
    dplyr::arrange(peak)
}




# Evolution models
evolve = function(copy_state, target)
{
  if(copy_state %>% as_karyotype() == target) return(copy_state %>% mutation() %>% list())
  
  cap = (target %>% strsplit(split = ':'))[[1]] %>% as.numeric %>% sum
  cap = 2 * cap
  
  current_state = original_state = list(copy_state)
  
  repeat{
    amp_new_state = lapply(current_state, amplify) %>% unlist(recursive = FALSE)
    del_new_state = lapply(current_state, delete) %>% unlist(recursive = FALSE)
    wgs_new_state = lapply(current_state, genome_double) %>% unlist(recursive = FALSE)
    
    # Filter by capping
    amp_new_state = amp_new_state[sapply(amp_new_state, as_ploidy) <= cap]
    del_new_state = del_new_state[sapply(del_new_state, as_ploidy) <= cap]
    wgs_new_state = wgs_new_state[sapply(wgs_new_state, as_ploidy) <= cap]
    
    current_state = amp_new_state %>%
      append(del_new_state) %>%
      append(wgs_new_state)
    
    what_we_have = current_state %>% sapply(as_karyotype)
    
    # print(what_we_have)
    
    if(target %in% what_we_have) {
      target_state = current_state[which(target == what_we_have)]
      break;
    }
  }
  
  # Remove duplicates - based on allele identities
  identities = sapply(target_state, function(x) x$allele %>% sort() %>% paste(collapse = ''))
  target_state = target_state[which(!duplicated(identities))]
  
  return(target_state %>% lapply(mutation))
}

########################################################################################


###################### Branching evolution (maybe not needed) ############################
branching_evolution = function(starting, left, right, CCF_1, purity)
{
  start = initial_state(starting)
  branch_left = start %>% evolve(left)
  branch_right = branch_left
  
  solutions = tidyr::expand_grid(L = seq_along(branch_left), R = seq_along(branch_left))
  solutions = lapply(1:nrow(solutions), function(i) {
    get_peaks(
      clone_1 = branch_left[[solutions$L[i]]],
      clone_2 = branch_right[[solutions$R[i]]],
      CCF_1,
      purity)
  })
  
  solutions_id = lapply(solutions, function(x) x$peak %>% paste(collapse= ';'))
  solutions = solutions[!duplicated(solutions_id)]
  
  lapply(solutions, function(x){
    geno_1 = x$genotype_1
    geno_1 = geno_1[!is.na(geno_1)] %>% unique
    geno_2 = x$genotype_2
    geno_2 = geno_2[!is.na(geno_2)] %>% unique
    
    x$genotype_initial = start %>% get_alleles() %>% sort() %>% paste(collapse = '')
    x$model = 'branching'
    x$model_id = paste0(x$genotype_initial[1], " -> ", geno_1, " | ", geno_2)
    x %>% mutate(role = ifelse(is.na(karyotype_1) | is.na(karyotype_2), "private", 'shared'))
  })
  
}

#############################################################


################## Linear evolution #######################################
linear_evolution = function(starting, first_child, second_child, CCF_1, purity)
{
  start = initial_state(starting)
  first_children = start %>% evolve(first_child) # generate tibble with possible first children
  second_children = lapply(first_children, function(x) evolve(x, second_child))
  
  solutions = lapply(first_children %>% seq_along, function(i) {
    second_children[[i]] %>% lapply(function(y) {
      get_peaks(clone_1 = first_children[[i]],
                clone_2 = y,
                CCF_1,
                purity)
    })
  }) %>% unlist(recursive=FALSE)
  
  solutions_id = lapply(solutions, function(x) x$peak %>% paste(collapse= ';'))
  solutions = solutions[!duplicated(solutions_id)]
  
  
  lapply(solutions, function(x){
    geno_1 = x$genotype_1
    geno_1 = geno_1[!is.na(geno_1)] %>% unique
    geno_2 = x$genotype_2
    geno_2 = geno_2[!is.na(geno_2)] %>% unique
    
    x$genotype_initial = start %>% get_alleles() %>% sort() %>% paste(collapse = '')
    x$model = 'linear'
    x$model_id = paste0(x$genotype_initial[1], " -> ", geno_1, " -> ", geno_2)
    x %>% mutate(role = ifelse(is.na(karyotype_1) | is.na(karyotype_2), "private", 'shared'))
  })
  #
}









