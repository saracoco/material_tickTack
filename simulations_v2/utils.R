.libPaths("~/R/rstudio_v3/")
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(VariantAnnotation)
library(stringr)
require(transport)
library(cluster)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(tickTack)


generate_multiplicities = function(k, tau, N_mutations, m=1) {
  if (k == "2:0") {
    # 2:0
    n1 = 2 * m * (1 - tau)
    n2 = m * tau
  } else if (k == "2:1") {
    # 2:1
    n1 = m * (3 - 2*tau)
    n2 = m * tau
  } else if (k == "2:2") {
    # 2:2
    n1 = 4 * m * (1 - tau)
    n2 = 2 * m * tau
  } else {
    stop("k not recognized")
  }
  muts <- c(rmultinom(n = 1, size = N_mutations, prob = c(n1, n2)))
  list(n1 = muts[1], n2 = muts[2])
}

simulate_dataset = function(N_events, N_clocks, N_mutations, pi, coverage, sigma_tau = .01, min_dist = .1) {
  taus = runif(N_clocks, .01, .99)
  while (!all(diff(sort(taus)) >= min_dist) & length(unique(taus)) != N_clocks) {
    taus = runif(N_clocks, .01, .99)  
  }
  taus_clust = taus = sample(taus, N_events, replace = TRUE)
  taus = lapply(taus, function(t) {rnorm(1, t, sigma_tau)}) %>% unlist()
  
  karyotypes = sample(c("2:0", "2:1", "2:2"), N_events, replace = TRUE)
  
  sim_data = lapply(1:N_events, function(idx) {
    m = 1
    tau = taus[idx]
    k = karyotypes[idx]
    
    mults = generate_multiplicities(k = k, tau = tau, N_mutations = N_mutations, m = 1)
    
    start_off = 10000
    
    Major = as.numeric(unlist(strsplit(k, ":"))[1])
    minor = as.numeric(unlist(strsplit(k, ":"))[2])
    
    cn = dplyr::tibble(
      chr = 1,
      startpos = 1 + (idx - 1) * N_mutations,
      endpos = (idx) * N_mutations,
      nMaj1_A=Major,
      nMin1_A=minor) %>%
      dplyr::mutate(startpos = startpos + start_off, endpos=endpos+start_off)
    
    mult = dplyr::tibble(
      chr = 1,
      end = cn$startpos:cn$endpos,
      no.chrs.bearing.mut = c(rep(1,mults$n1), rep(2,mults$n2))
    )
    
    muts = dplyr::tibble(
      chr = 1,
      start = cn$startpos:cn$endpos,
      end = cn$startpos:cn$endpos,
      ref = "G"
    )
    muts$alt = "A"
    
    # Add DP and
    peaks <- get_clonal_peaks(k, pi)
    
    muts$DP = coverage
    muts$NV = lapply(1:nrow(muts), function(j) {
      p = peaks[mult[j,]$no.chrs.bearing.mut]
      rbinom(1, size = coverage, prob = p)
    }) %>% unlist()
    muts$VAF = muts$NV / muts$DP
    
    list(cn = cn, mult=mult, muts=muts)
  })
  
  sim_data
  
  sim_cn <- lapply(1:N_events, function(i) { sim_data[[i]]$cn }) %>% do.call("rbind", .)
  sim_mult <- lapply(1:N_events, function(i) { sim_data[[i]]$mult }) %>% do.call("rbind", .)
  sim_muts <- lapply(1:N_events, function(i) { sim_data[[i]]$muts }) %>% do.call("rbind", .)
  
  list(cn=sim_cn, mult=sim_mult, muts=sim_muts, true_taus = taus, taus_clust=taus_clust)
}

fit_AmpTimeR = function(sim) {
  N_events = nrow(sim$cn)
  lapply(1:N_events, function(idx) {
    cn = sim$cn[idx, ]
    mult = sim$mult
    muts = sim$muts %>% dplyr::select(chr, start, end, ref, alt)
    
    segment_time <- time_amplification(
      cn_data = cn %>% as.data.frame(),
      multiplicity_data = mult %>% dplyr::filter(chr == cn$chr, end >= cn$startpos, end <= cn$endpos) %>% as.data.frame(),
      mutation_data = muts %>% dplyr::filter(chr == cn$chr, end >= cn$startpos, end <= cn$endpos) %>% as.data.frame(),
      muts_type = "All",
      sample_id = "test_sample",
      amplification_chrom = cn$chr,
      amplification_start = cn$startpos,
      amplification_stop = cn$endpos,
      is_WGD = TRUE,
      genome = "hg19"
    )
    
    dplyr::tibble(
      segment_idx = idx, 
      tau = segment_time$t_1_mean_bootstrap, 
      tau_low = segment_time$t_1_lower_ci, 
      tau_high = segment_time$t_1_upper_ci, 
      model = "AmpTimeR"
    )
  }) %>% do.call("bind_rows", .)
}

fit_tickTack_single = function(sim, pi, min_mutations) {
  N_events = nrow(sim$cn)
  pcawg_example$cna
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end)
  
  tickTack_single_res <- tickTack::fit(cn, mutations = muts, purity = pi, alpha = .05, min_mutations_number = min_mutations, beta_binomial = F)
  tickTack_single_res
}


fit_tickTack_h = function(sim, pi, min_mutations, tolerance, INIT, max_attempts) {
  N_events = nrow(sim$cn)
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos) %>% 
    dplyr::mutate(CCF = 1)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end) %>% 
    dplyr::mutate(CCF = 1)
  
  x = list(
    cna = cn, 
    mutations = muts,
    metadata = dplyr::tibble(sample = "sample", purity=pi)
  )
  
  x <- tickTack::fit_h(x, max_attempts=max_attempts, INIT=INIT, tolerance = tolerance)
  
  results_simulated <- x$results_timing
  results_model_selection <- tickTack::model_selection_h(results_simulated)
  best_K <- results_model_selection$best_K
  results_tickTack <- list( results = results_simulated, results_model_selection = results_model_selection)
  
  results_tickTack
}

fit_MutTimeR <- function(sim, pi) {
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos) %>% 
    dplyr::mutate(CCF = 1)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end) %>% 
    dplyr::mutate(CCF = 1)
  
  x = list(
    cna = cn, 
    mutations = muts,
    metadata = dplyr::tibble(sample = "sample", purity=pi)
  )
  
  
  vcf <- convert_to_vcf(x$mutations, x$cna, x$metadata)
  bb <- convert_to_granges(x$cna, x$mutations, x$metadata)
  clusters <- data.frame(cluster=1:1, proportion=c(pi), n_ssms=c(100))
  
  mt <- MutationTimeR::mutationTime(vcf, bb, clusters=clusters)
  mcols(bb) <- cbind(mcols(bb),mt$T)
  
  res_MutTime <- list(vcf = vcf, cn_timed = bb)  
  res_MutTime
}




convert_to_vcf <- function(mutations, cna, metadata) {
  # Validate input
  if(missing(mutations) || missing(cna) || missing(metadata)) {
    stop("mutations, cna, and metadata are required")
  }
  
  # Create VRanges object
  v <- VRanges(
    seqnames = unlist(lapply(mutations$chr, function(c) {str_replace_all(c, "chr", "")})),
    ranges = IRanges(start = mutations$from, end = mutations$to),
    ref = mutations$ref,
    # alt = mutations$alt,
    totalDepth = mutations$DP,
    altDepth = mutations$NV,
    QUAL = rep(NA_real_, nrow(mutations)),
    FILTER = rep("PASS", nrow(mutations))
  )
  
  # Set sample name
  sampleNames(v) <- metadata$sample
  
  # Convert to VCF
  vcf <- as(v, "VCF")
  
  # Create custom geno header
  geno_header <- DataFrame(
    Number = c("2", "1", "1"),
    Type = c("Integer", "Integer", "String"),
    Description = c(
      "Allelic depths (number of reads in each observed allele)",
      "Total read depth",
      "Variant filters"
    ),
    row.names = c("AD", "DP", "FT")
  )
  
  # Add custom info headers
  info(header(vcf)) <- rbind(
    info(header(vcf)),
    DataFrame(
      Number = 1,
      Type = rep("Integer", 2),
      Description = c("Tumour ref count", "Tumour alt count"),
      row.names = c("t_ref_count", "t_alt_count")
    )
  )
  
  # Set geno header
  geno_header <- DataFrame(
    Number = c("2", "1", "1"),
    Type = c("Integer", "Integer", "String"),
    Description = c(
      "Allelic depths of reference and alternate alleles",
      "Total read depth",
      "Variant filters"
    ),
    row.names = c("AD", "DP", "FT")
  )
  geno(vcf)
  geno(header(vcf)) <- geno_header
  
  # Prepare genotype data
  geno(vcf)$AD <- as.matrix(rep(NA, nrow(mutations)))
  geno(vcf)$DP <- as.matrix(mutations$DP)
  geno(vcf)$FT <- as.matrix(rep(NA, nrow(mutations)))
  
  # Add info fields
  info(vcf)$t_alt_count <- altDepth(v)
  info(vcf)$t_ref_count <- totalDepth(v) - altDepth(v)
  
  return(vcf)
}

convert_to_granges <- function(cna, mutations, metadata) {
  # Validate input
  if(missing(cna) || missing(mutations) || missing(metadata)) {
    stop("All three inputs (cna, mutations, metadata) are required")
  }
  
  # Create GRanges
  cna_granges <- GRanges(
    seqnames = gsub("^chr", "", cna$chr),
    ranges = IRanges(start = cna$from, end = cna$to),
    strand = "*",
    major_cn = cna$Major,
    minor_cn = cna$minor,
    clonal_frequency = metadata$purity
  )
  
  return(cna_granges)
}


plot_heatmap <- function(matrix_data, 
                         x_label = "Columns", 
                         y_label = "Rows", 
                         fill_label = "Value", 
                         low_color = "blue", 
                         high_color = "red", 
                         title = "Heatmap") {
  # Check if input is a matrix
  if (!is.matrix(matrix_data)) {
    stop("Input must be a matrix.")
  }
  
  # Convert matrix to a long-format data frame
  df <- melt(matrix_data)
  colnames(df) <- c("Row", "Column", "Value")
  
  # Create the heatmap
  ggplot(df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = low_color, high = high_color, name = fill_label) +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}

compute_metrics <- function(posterior_draws) {
  # Input: 
  # posterior_draws: A list where each element contains posterior draws for a sample
  
  # Output:
  # A matrix with Wasserstein distance, KL divergence, and overlap-based similarity
  
  n <- length(posterior_draws)
  results <- matrix(NA, nrow = n, ncol = n, 
                    dimnames = list(paste0("Sample", 1:n), paste0("Sample", 1:n)))
  
  wasserstein <- results
  kl_divergence <- results
  overlap_similarity <- results
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Extract draws for sample i and j
      draws_i <- posterior_draws[[i]]
      draws_j <- posterior_draws[[j]]
      
      # Compute Wasserstein distance
      wasserstein[i, j] <- wasserstein1d(draws_i, draws_j)
      wasserstein[j, i] <- wasserstein[i, j]
      
      # Compute KL divergence using density estimates
      density_i <- density(draws_i)
      density_j <- density(draws_j)
      
      # Interpolate densities to match the support
      common_x <- sort(unique(c(density_i$x, density_j$x)))
      p <- approx(density_i$x, density_i$y, xout = common_x, rule = 2)$y
      q <- approx(density_j$x, density_j$y, xout = common_x, rule = 2)$y
      
      # Normalize to avoid NaN issues
      p <- p / sum(p)
      q <- q / sum(q)
      
      kl_divergence[i, j] <- sum(ifelse(p > 0 & q > 0, p * log(p / q), 0))
      kl_divergence[j, i] <- sum(ifelse(q > 0 & p > 0, q * log(q / p), 0))
      
      # Compute overlap-based similarity
      overlap_similarity[i, j] <- sum(pmin(p, q))
      overlap_similarity[j, i] <- overlap_similarity[i, j]
    }
  }
  
  list(
    Wasserstein = wasserstein,
    KL_Divergence = kl_divergence,
    Overlap_Similarity = overlap_similarity
  )
}

compute_WD <- function(posterior_draws) {
  # Input: 
  # posterior_draws: A list where each element contains posterior draws for a sample
  
  # Output:
  # with Wasserstein distance
  
  n <- length(posterior_draws)
  results <- matrix(NA, nrow = n, ncol = n, 
                    dimnames = list(paste0("Sample", 1:n), paste0("Sample", 1:n)))
  
  wasserstein <- results
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Extract draws for sample i and j
      draws_i <- posterior_draws[[i]]
      draws_j <- posterior_draws[[j]]
      
      # Compute Wasserstein distance
      wasserstein[i, j] <- wasserstein1d(draws_i, draws_j)
      wasserstein[j, i] <- wasserstein[i, j]
    }
  }
  
  
  wasserstein
}

hierarchical_optimal_clusters <- function(metrics_matrix, max_clusters = 10, method = "ward.D2") {
  # Check if the input is a valid matrix
  if (!is.matrix(metrics_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Convert the metrics matrix to a distance object
  dist_matrix <- as.dist(metrics_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = method)
  
  # Initialize variables to store evaluation results
  silhouette_scores <- numeric(max_clusters - 1)
  
  # Loop over the number of clusters (from 2 to max_clusters)
  for (k in 2:max_clusters) {
    # Cut the dendrogram into k clusters
    cluster_assignments <- cutree(hc, k = k)
    
    # Calculate silhouette scores for the current k
    silhouette_result <- silhouette(cluster_assignments, dist_matrix)
    silhouette_scores[k - 1] <- mean(silhouette_result[, 3]) # Mean silhouette width
  }
  
  # Find the optimal number of clusters (max silhouette score)
  optimal_clusters <- which.max(silhouette_scores) + 1
  
  # Return results
  return(list(
    optimal_clusters = optimal_clusters,
    silhouette_scores = silhouette_scores,
    hc = hc # Return the hierarchical clustering object for plotting or further analysis
  ))
}

cluster_hierarchical <- function(metrics_matrix, num_clusters) {
  # Check if the input is a valid matrix
  if (!is.matrix(metrics_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Convert the metrics matrix to a distance object
  dist_matrix <- as.dist(metrics_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  # Cut the dendrogram into clusters
  clusters <- cutree(hc, k = num_clusters)
  
  # Return the cluster assignments
  return(clusters)
}

compute_overlap_matrix <- function(intervals) {
  N <- nrow(intervals)
  overlap_matrix <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        overlap_matrix[i, j] <- 1  # An interval fully overlaps with itself
      } else {
        a1 <- intervals[i, 1]
        b1 <- intervals[i, 2]
        a2 <- intervals[j, 1]
        b2 <- intervals[j, 2]
        
        intersection_start <- max(a1, a2)
        intersection_end <- min(b1, b2)
        intersection_length <- max(0, intersection_end - intersection_start)
        
        length1 <- b1 - a1
        length2 <- b2 - a2
        
        if (length1 > 0 && length2 > 0) {
          overlap_matrix[i, j] <- intersection_length / min(length1, length2)
        } else {
          overlap_matrix[i, j] <- 0  # Avoid division by zero
        }
      }
    }
  }
  
  return(overlap_matrix)
}
                      
# Cluster AmpTimeR
clusterAmpTimeR = function(fp) {
  res = readRDS(paste0(fp, '/res_AmpTimeR.rds'))
  intervals = cbind(res$tau_low, res$tau_high)
  similarity_matrix = compute_overlap_matrix(intervals)
  best_k = hierarchical_optimal_clusters(1-similarity_matrix, max_clusters = nrow(intervals) - 1)$optimal_clusters
  clusters = cluster_hierarchical(1 - similarity_matrix, num_clusters = best_k)  
  clusters
}

# Cluster MutTimeR
clusterMutTimeR = function(fp) {
  res = readRDS(paste0(fp, '/res_MutTime.rds'))
  intervals = cbind(res$cn_timed$time.lo, res$cn_timed$time.up)
  similarity_matrix = compute_overlap_matrix(intervals)
  best_k = hierarchical_optimal_clusters(1 - similarity_matrix, max_clusters = nrow(similarity_matrix) - 1)$optimal_clusters
  clusters = cluster_hierarchical(1 - similarity_matrix, num_clusters = best_k)  
  clusters  
}

# Cluster tickTack single
cluster_tickTack_single = function(fp) {
  res = readRDS(paste0(fp, '/res_tickTack_single.rds'))
  draws = lapply(res$inference_results$segment %>% unique(), function(idx) {
    res$inference_result %>% 
      dplyr::filter(segment == idx) %>% 
      dplyr::filter(tau <= 1, tau >= 0) %>% 
      pull(tau)
  })
  similarity_matrix = compute_WD(draws)
  similarity_matrix[is.na(similarity_matrix)] = 0
  best_k = hierarchical_optimal_clusters(similarity_matrix, max_clusters = nrow(similarity_matrix) - 1)$optimal_clusters
  clusters = cluster_hierarchical(similarity_matrix, num_clusters = best_k)  
  clusters    
}

cluster_tickTack_h = function(fp) {
  res = readRDS(paste0(fp, '/res_tickTack_h.rds'))
  res$results_model_selection$best_fit$summarized_results$clock_mean
}
