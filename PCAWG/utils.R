
TSGs <- c("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC", "CDKN2A", "SMAD4", "VHL", "NF1")
Oncogenes <- c("MYC", "KRAS", "BRAF", "EGFR", "HER2", "ALK", "PIK3CA", "ABL1", "CCND1", "NRAS")
DNA_repair = c("RAD51")
#genes_of_interest <- c(TSGs, Oncogenes)
genes_of_interest <- c(TSGs, Oncogenes, DNA_repair)

is_fittable = function(fit, possible_k = c("2:1", "2:2", "2:0"), alpha = .05, min_mutations_number = 2) {
  get_clonal_peaks = function(k, purity) {
    multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
    major <- multiplicities[1]
    n_tot <- sum(multiplicities)
    # get only Major and 1
    multiplicities <- c(1, major)
    peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
    return(sort(peaks))
  }
  
  segments = fit$cna
  mutations = fit$snvs
  purity = fit$purity
  
  segments <- segments %>%
    tidyr::drop_na(Major, minor)
  
  n_segments <- nrow(segments)
  
  accepted_segment_idx <- 0
  for (segment_idx in 1:n_segments) {
    
    # Segments
    segment <- segments[segment_idx, ]
    chr <- segment$chr
    
    segment_id <- paste(chr, segment$from, segment$to, sep = "_")
    
    # Get karyotype
    Major <- segment$Major
    minor <- segment$minor
    
    k <- paste(Major, minor, sep=':')
    
    peaks <- get_clonal_peaks(k, purity)
    
    if (k %in% possible_k & (segment$to-segment$from)>1e6 & purity >= 0.4) {
      # Get info for mutations
      segment_mutations <- mutations %>%
        dplyr::filter(chr == segment$chr,.data$from > segment$from, .data$to < segment$to) %>%
        tidyr::drop_na(DP)
      
      
      
      accepted_mutations <- data.frame()
      if (nrow(segment_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV
        
        accepted_idx <- lapply(1:length(DP), function(i) {
          #fisso i picchi per il segmento e vedo se tutte le mutazioni ricadono in almeno uno dei due intervalli intorno ai picchi
          # che sono diversi a seconda del valore di DP per la specifica mutazione
          for (p in peaks) {
            quantiles <- stats::qbinom(probs, DP[i], p)
            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
        }) %>% unlist()
        
        # Get only good mutations
        accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
        
      }
      
      if (nrow(accepted_mutations) >= min_mutations_number) {
        return(TRUE)
        
      }
    }
  }
  return(FALSE)
}

parse_summarized_results <- function(results) {
  results %>% 
    tidyr::separate(segment_name, into = c("chr", "from", "to"), sep = "_", convert = TRUE) %>% 
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>% 
    dplyr::select(chr, from, to, segment_id, karyotype, clock_mean, clock_low, clock_high, clock_rank)
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
