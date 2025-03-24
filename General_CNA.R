# General CNA algorithm: given a complex CNA state, simulate all possible trajectories that 
# could have led to that state respecting the rules:
# 1. Allowed moves are:
#         - Aplification of one allele
#         - Deletion of one allele 
#         - WGD
# 2. If I have a 0, the allele cannot duplicate
# 3. Maximum parsimony: the final state is reached in the minimum number of steps, 
#                       i.e. the maximum number of steps necessary to reach the final step 
#                       proceeding one step at a time
# Input:
# - starting karyotype 
# - number (proportion) of mutations associated to each peak-multiplicity
# Output:
# - set of steps (ex A1B1->A1A2B1->...) 
# - psudotimes corresponding to each model


compute_pseudotimes = function(starting='1:1', target=)
{
  # Are we starting from a karyo that lacks one of the alleles? (in which case I cannot amplify that allele)
  any_LOH = (strsplit(starting, split = ':')[[1]] == "0") %>% any
  # Are there LOHs in the final karyotypes?
  no_LOH1 = (strsplit(karyotype_1, split = ':')[[1]] != "0") %>% all
  no_LOH2 = (strsplit(karyotype_2, split = ':')[[1]] != "0") %>% all
  # If there is an LOH in the starting karyotype but the final karyotypes don't have LOHs stop. There is
  # no way in which is possible to reach those karytypes from an LOH
  if(any_LOH & (no_LOH1 | no_LOH2))
  {
    cli::cli_abort("No evolution model can reach {.field {karyotype_1}} /
    {.field {karyotype_2}} from {.field {starting}}.
                   Rerun with with different parameters, aborting!")
  }
  
  ####### Utility functions: skip this and go to the applications of these functions,
  ####### come back when you encounter each of them to understand what they do 
  # getters
  # outputs the vector of mutations of the dataframe copy state or referred 
  # to a single allele
  get_mutations = function(copy_state, allele = NULL){
    if(is.null(allele))
      mutations = copy_state %>% pull(mutations) %>% unlist()
    else
      mutations = copy_state %>% filter(allele == !!allele) %>% pull(mutations) %>% unlist()
    
    return(mutations)
  }
  
  get_alleles = function(copy_state){
    copy_state$allele
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
  
  # Add mutations to all alleles
  # adds a new mutation (unique string of characters) to 
  # each row of the dataframe copystate
  mutation = function(copy_state){
    for(i in 1:nrow(copy_state))
    {
      muts = copy_state$mutations[[i]] %>% unlist()
      muts = list(c(muts, copy_state %>% new_mutation))
      copy_state$mutations[[i]] = muts
    }
    copy_state
  }
  
  # Amplify allele
  amplify =  function(copy_state){
    augment = function(which_allele)
    {
      mutations_allele = copy_state %>% get_mutations(allele = which_allele)
      
      maj_min = substr(which_allele, 0, 1)
      n_allele = gsub(x = which_allele, 'A', "") %>% gsub(pattern = 'B', replacement = "") %>%
        as.numeric() # counts the multiplicity of the allele
      
      # New allele can be +1, unless there are other alleles with larger number
      new_allele = paste(maj_min, n_allele + 1, sep = '')
      # If the new allele has already been amplified, change its index 
      while(new_allele %in% (copy_state %>% get_alleles())) {
        n_allele = n_allele + 1
        new_allele = paste(maj_min, n_allele + 1, sep = '')
      }
      
      new_entry = copy_state %>% filter(allele == which_allele)
      new_entry$allele = new_allele # the new allele inherits mutations from the allele it originated from,
      # might need to change this 
      
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
  
  # Genome_double
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
  
  # Multiplicity
  multiplicities = function(copy_state){
    copy_state %>% get_mutations() %>% table()
  }
  
  # Get stast
  as_karyotype = function(copy_state){
    tab_counts = copy_state %>% get_alleles() %>% substr(0, 1) %>% table() %>% sort(decreasing = TRUE)
    state = as.numeric(tab_counts) %>% paste(collapse = ':')
    if(!grepl(":", state)) state = paste0(state, ':0')
    
    return(state)
  }
  
  as_ploidy = function(copy_state){
    copy_state %>% nrow()
  }
  
  # Evolution models - this is the key function 
  # takes in input a copy_state df and a target, and creates a list of dataframes 
  # with all the possible 1-step evolutions that leads to the target
  evolve = function(copy_state, target)
  {
    # If the copystate is already the same as the target, return the mutated copystate
    if(copy_state %>% as_karyotype() == target) return(copy_state %>% mutation() %>% list())
    # Compute the maximum number of steps
    cap = (target %>% strsplit(split = ':'))[[1]] %>% as.numeric %>% sum
    cap = 2 * cap
    
    current_state = original_state = list(copy_state)
    
    repeat{
      # Evolve the current state(s) computing the possible configuration(s) for:
      # - amplification of one allele
      # - deletion of one allele 
      # - WGD
      amp_new_state = lapply(current_state, amplify) %>% unlist(recursive = FALSE)
      del_new_state = lapply(current_state, delete) %>% unlist(recursive = FALSE)
      wgs_new_state = lapply(current_state, genome_double) %>% unlist(recursive = FALSE)
      
      # Filter by capping - make sure none of the new states exceeds the 
      # maximum allowed ploidy
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
  
  ## Initial state
  # When target = starting state (1:1), creates a dataframe with the starting alleles and initialises a list of mutations for each
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
  
  # Models 
  # Creates the branching evolution model, left and right are initialised as karyotype_1 and karyotype_2
  branching_evolution = function(starting, left, right, CCF_1, purity)
  {
    # Initialise the evolution - creates a dataframe with alleles and corresponding mutations 
    start = initial_state(starting)
    # Since the evolution is branched, the two branches will develop independently 
    # The evolve function creates a list of dataframes, each storing 
    # the allele and mutation information of one of the corresponding trajectory
    branch_left = start %>% evolve(left) 
    branch_right = start %>% evolve(right)
    
    solutions = tidyr::expand_grid(L = seq_along(branch_left), R = seq_along(branch_right))
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
  
  
  ############### Apply the models
  
  b_model = suppressWarnings(branching_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
    Reduce(f = bind_rows)
  
  l1_model = NULL
  if((grepl("0", karyotype_1) & grepl("0", karyotype_2)) | !grepl("0", karyotype_1))
    l1_model = suppressWarnings(linear_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
    Reduce(f = bind_rows)
  
  l2_model = NULL
  if((grepl("0", karyotype_2) & grepl("0", karyotype_1)) | !grepl("0", karyotype_2))
    l2_model = suppressWarnings(linear_evolution(starting, karyotype_2, karyotype_1, 1-CCF_1, purity)) %>%
    Reduce(f = bind_rows)
  
  return(bind_rows(b_model, l1_model, l2_model))
}
