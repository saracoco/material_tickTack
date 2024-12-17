place_mutations <- function(forest, karyo, chromosomes = c("5","6","10","12","13", "22"), taus, long_or_short_vec = long_or_short_vec){
  
  ###### DRIVERS ######
  
  
  ##################################
  
  
  generate_segment = function(coords, short_or_long, current_segments = NULL) {
    
    status = ""
    while (status != "PASS") {
      
      status = 'PASS'
      sampled_chr <- sample(coords$chr, 1)
      chr_info = coords %>% dplyr::filter(.data$chr == sampled_chr)
      
      if (short_or_long == "short") {
        l = runif(1, 100, 10000)
      } else {
        l = runif(1, 10000, min(1e7, chr_info$length))
      }
      
      start = runif(1, 1, chr_info$length - l)
      
      if (!is.null(current_segments)) {
        status = "PASS"
        
        d <- current_segments %>%
          dplyr::filter(chr == sampled_chr)
        
        if (nrow(d) > 0) {
          
          previous_starts = d$start
          previous_ends = d$start + d$length
          
          new_start = start
          new_end = start + l
          
          if (any((previous_starts >= new_start) & (previous_ends <= new_end))) status == ""
          if (any((new_start >= previous_starts) & (new_start <= previous_ends))) status == ""
          if (any((new_end >= previous_starts) & (new_end <= previous_ends))) status == ""
        }
      }
      
      res = dplyr::tibble(chr = sampled_chr, start = start, length = l)
    }
    return(res)
  }
  
  
  
  
  
  generate_segments = function(long_or_short_vec, chromosomes = chromosomes) {
    
    curr_segments = NULL
    chromosomes = unlist(lapply(chromosomes, function(x) {paste0("chr", x)}))
    coords = tickTack::chr_coordinates_GRCh38 %>% dplyr::filter(chr %in% chromosomes)
    l_or_s = long_or_short_vec[3]
    for (l_or_s in long_or_short_vec) {
      new_segment <- generate_segment(coords, l_or_s, current_segments = curr_segments)
      curr_segments <- dplyr::bind_rows(
        curr_segments,
        new_segment
      ) %>% dplyr::mutate(chr = gsub("chr", "", chr))
    }
    curr_segments
  }
  
  
  segs <- generate_segments(long_or_short_vec, chromosomes = chromosomes)
  
  d_taus <- dplyr::tibble(tau = taus, karyotype = karyo) %>%
    dplyr::mutate(tau_factor = as.numeric(factor(tau))) %>%
    dplyr::mutate(Clone = paste0("Clone " , tau_factor)) %>%
    dplyr::arrange(tau)
  
  events <- cbind(segs, d_taus)
  saveRDS(events, "results/true_events.rds")
  
  
  ###### MUTATION ENGINE ######
  
  dir <- getwd()
  setwd("/orfeo/scratch/cdslab/shared/races/GRCh38/")
  m_engine <- MutationEngine(setup_code = "GRCh38")
  setwd(dir)
  nAllele = 2
  
  # clone = "Clone 1"
  for (clone in unique(events$Clone)) {
    tmp_events <- events %>% filter(Clone == clone)
    # nAllele = 2
    
    list_CNA <- list()
    for (i in 1:nrow(tmp_events)){
      curr_event <- tmp_events[i,]
      
      if (curr_event$karyotype=="2:1"){
        
        chr <-  curr_event$chr
        len <-  curr_event$length
        chr_pos <-  curr_event$start
      
        list_CNA <- c(list_CNA, rRACES::CNA(type = "A", chr = chr, chr_pos = chr_pos, len = len, src_allele = 0, allele = nAllele))
        nAllele = nAllele + 1# 2:1
        
      }else if(curr_event$karyotype=="2:2"){
        chr <-  curr_event$chr
        len <-  curr_event$length
        chr_pos <-  curr_event$start
        
        list_CNA <- append(list_CNA,  rRACES::CNA(type = "A", chr = chr, chr_pos = chr_pos, len = len, src_allele = 0, allele = nAllele)) # 2:2
        nAllele = nAllele + 1
        list_CNA <- append(list_CNA,  rRACES::CNA(type = "A", chr = chr, chr_pos = chr_pos, len = len, src_allele = 1, allele = nAllele)) # 2:2
        nAllele = nAllele + 1
        
        
      }else{ #2:0
        
        chr <-  curr_event$chr
        len <-  curr_event$length
        chr_pos <-  curr_event$start
        
        list_CNA <- append(list_CNA,  rRACES::CNA(type = "A", chr = chr, chr_pos = chr_pos, len = len, src_allele = 1, allele = nAllele)) # 2:0 
        nAllele = nAllele + 1
        
        list_CNA <- append(list_CNA,  rRACES::CNA(type = "D", chr = chr, chr_pos = chr_pos, len = len, allele = 0)) # 2:0
        
      }
    }
    
    mu_SNV <- 5e-8
    mu_CNA <- 0
    
    m_engine$add_mutant(mutant_name =  paste0(clone),
                          passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                          drivers = list_CNA)
  }
  
  list_CNA = list()
  m_engine$add_mutant(mutant_name =  paste0("Clone ", 0),
                      passenger_rates = c(SNV=5e-8, CNA=0),
                      drivers = list_CNA)
  
  
  m_engine$add_mutant(mutant_name =  paste0("Final clone"),
                      passenger_rates = c(SNV=5e-8, CNA=0),
                      drivers = list_CNA)
  
  
  
  # events <- place_events(possible_positions, possible_len, chromosomes)
  # associate rank dplyr to each event 
  ###############################################
  # transform taus in the index that tau has in vector_tau
  
  
  ###### SIGNATURE ######
  # treatment_info <- readRDS("data/treatment_info.rds")
  m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5)) # ID1 is missing
  # m_engine$add_exposure(time = treatment_info$treatment_start, c(SBS5 = 0.4, SBS1 = 0.4, SBS25 = 0.2))
  # m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5))
  # 
  ###### PHYLO FOREST ######
  # forest <- load_samples_forest("data/samples_forest.sff")                      # errore in SPN04
  phylo_forest <- m_engine$place_mutations(forest,
                                           num_of_preneoplatic_SNVs = 800,        # need to set to zero?
                                           num_of_preneoplatic_indels = 200)
  # phylo_forest$save("data/phylo_forest.sff")
  
  
  
  
  all_SNV <- phylo_forest$get_sampled_cell_mutations() %>% as_tibble()
  all_SNV %>%
    group_by(cause) %>%
    summarise(nPos = n_distinct(chr_pos)) %>%
    print()
  
  all_SNV %>%
    group_by(class) %>%
    summarise(nPos = n_distinct(chr_pos)) %>%
    print()
  
  annot_forest <- plot_forest(forest) %>%
    annotate_forest(phylo_forest,
                    samples = T,
                    MRCAs = T,
                    exposures = T,
                    drivers=T,
                    add_driver_label = T)
  
  exp_timeline <- plot_exposure_timeline(phylo_forest)
  
  
  labels <- get_relevant_branches(forest)
  sticks <- plot_sticks(forest, labels)
  
  # ggplot2::ggsave('plots/sticks.png', plot = sticks, width = 210, height = 297, units = "mm", dpi=300)
  # 
  # ggplot2::ggsave('plots/annot_forest.png', plot = annot_forest, width = 210, height = 297, units = "mm", dpi=300)
  # ggplot2::ggsave('plots/exp_timeline.png', plot = exp_timeline, width = 210, height = 297, units = "mm", dpi=300)
  # 
  
  
  
  
  pl <- annot_forest + sticks + exp_timeline + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
  
  # ggplot2::ggsave('plots/mutations.png', plot = pl, width = 210, height = 297, units = "mm", dpi=300)
  # ggplot2::ggsave('plots/mutations.pdf', plot = pl, width = 210, height = 297, units = "mm", dpi=300)
  # 
  
  results = list(sticks = sticks, exp_timeline = exp_timeline, pl = pl, phylo_forest = phylo_forest, all_SNV = all_SNV)
}
