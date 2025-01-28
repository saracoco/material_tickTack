
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(VariantAnnotation)
library(stringr)
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
  while (!all(diff(sort(taus)) >= min_dist)) {
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
      tau_high = segment_time$t_1_lower_ci, 
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
