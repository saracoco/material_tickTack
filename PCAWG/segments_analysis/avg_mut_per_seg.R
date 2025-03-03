
avg_mut_per_seg = function(fit, possible_k = c("2:1", "2:2", "2:0"), alpha = .05, min_mutations_number = 2) {
  
  get_clonal_peaks = function(k, purity) {
    multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
    major <- multiplicities[1]
    n_tot <- sum(multiplicities)
    multiplicities <- c(1, major)
    peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
    return(sort(peaks))
  }
  
  ttype = unique(fit$snvs$ttype)[2]
  segments = fit$cna
  mutations = fit$snvs
  purity = fit$purity
  
  segments <- segments %>%
    tidyr::drop_na(Major, minor)
  
  n_segments <- nrow(segments)
  
  accepted_segments <- 0
  mutations_per_segment <- c()
  
  for (segment_idx in 1:n_segments) {
    
    segment <- segments[segment_idx, ]
    chr <- segment$chr
    
    segment_id <- paste(chr, segment$from, segment$to, sep = "_")
    
    Major <- segment$Major
    minor <- segment$minor
    k <- paste(Major, minor, sep=':')
    
    peaks <- get_clonal_peaks(k, purity)
    
    if (k %in% possible_k) {
      segment_mutations <- mutations %>%
        dplyr::filter(chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
        tidyr::drop_na(DP)
      
      accepted_mutations <- data.frame(DP = numeric(), NV = numeric())  # Ensure it initializes correctly
      
      if (nrow(segment_mutations) > 0) {
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV
        
        accepted_idx <- lapply(1:length(DP), function(i) {
          for (p in peaks) {
            if (is.na(p) || is.na(DP[i])) return(NULL)  # Handle NA values
            
            quantiles <- stats::qbinom(probs, DP[i], p)
            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
          return(NULL)
        }) %>% unlist()
        
        if (!is.null(accepted_idx) && length(accepted_idx) > 0) {
          accepted_segments <- accepted_segments + 1
          accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
          mutations_per_segment <- c(mutations_per_segment, nrow(accepted_mutations))
        }
      }
    }
  }
  
  if (accepted_segments > 0) {
    return(list(avg_mut_x_seg = mean(mutations_per_segment, na.rm = TRUE), accepted_segments = accepted_segments, cancer_type = ttype))
  } else {
    return(FALSE)
  }
}
