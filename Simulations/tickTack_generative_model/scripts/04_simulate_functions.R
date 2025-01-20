get_major = function (k, get_minor = FALSE){
  Major = 0
  minor = 0
  
  k1 = stringr::str_split_1(k, pattern = ":")[1]
  k2 = stringr::str_split_1(k, pattern = ":")[2]
  
  if (k1 > k2){
    Major = k1
    minor = k2} else{
      Major = k2
      minor = k1}
  major_minor = list(Major=Major, minor = minor)
  
  if (get_minor == TRUE){
    return(minor)
  }else{
    return(Major)}
}




#' get_taus_karyo Function
#'
#' This function allows you to obtain simulated data for copy number events (tau and karyo input for get_simulation or simulate_and_fit) inputting the number of events and the proportions of taus and karyos.
#' @param number_events number of events to be simulated
#' @param vector_taus unique tau to be present in the simulation
#' @param vector_karyo unique karyotype to be present in the simulation
#' @param weignths_taus proportion of tau to generate for each unique value
#' @param weights_karyo proportion of karyo to generate for each unique karyo
#' @keywords simulation
#' @export
#' @examples
#' get_taus_karyo()

get_taus_karyo = function (number_events,
                                vector_tau, vector_karyo,
                                weigths_tau, weights_karyo ){

  taus <- sample( vector_tau[1:length(vector_tau)], number_events, replace=TRUE, prob=weigths_tau )
  karyo <- sample( vector_karyo[1:length(vector_karyo)], number_events, replace=TRUE, prob=weights_karyo )
  

  return(list(taus=taus,karyo=karyo))

}



#' get_taus_karyo gerstrung
#'
#' This function allows you to obtain simulated data for copy number events (tau and karyo input for get_simulation or simulate_and_fit) inputting the number of events and the proportions of taus and karyos.
#' @param number_events number of events to be simulated
#' @param vector_taus unique tau to be present in the simulation
#' @param vector_karyo unique karyotype to be present in the simulation
#' @param weignths_taus proportion of tau to generate for each unique value
#' @param weights_karyo proportion of karyo to generate for each unique karyo
#' @keywords simulation
#' @export
#' @examples
#' get_taus_karyo()

get_taus_karyo_gerstrung = function (number_events,
                           vector_tau, vector_karyo,
                           weigths_tau, weights_karyo, chromosomes=c(1:22) ){
  
  taus <- sample( vector_tau[1:length(vector_tau)], number_events, replace=TRUE, prob=weigths_tau )
  karyo <- sample( vector_karyo[1:length(vector_karyo)], number_events, replace=TRUE, prob=weights_karyo )
  chr <- sample(chromosomes, number_events, replace = FALSE)
  
  return(list(taus=taus,karyo=karyo, chr=chr))
  
}




#' get_simulation Function
#'
#' This function allows you to obtain simulated data for copy number events.
#' @param taus = a vector of doubles between [0,1]
#' @param karyotypes  = a vector of the same length of taus of strings of the type "2:1", "2:0", "2:2"
#' @param purity sample purity
#' @param time_interval time interval not in 0-1, real number
#' @param l length of the segment considered
#' @param mu mutation rate
#' @param w cell division rate
#' @param coverage average number of reads that align to a reference base
#' @keywords simulation
#' @export
#' @examples
#' get_simulation()

get_simulation = function (taus, karyotypes, purity = 0.9, time_interval = 20, l = 1e7, mu = 1e-4, w = 1e-2, coverage = 100){

  S = length(karyotypes)
  names <- paste("segment", 1:S,sep = " ")
  data_all_segments <- dplyr::tibble()
  

  for (j in 1:S) {
    tau <- taus[j]
    karyotype <- karyotypes[j]
    #print(tau)
    data_single_segment <- simulate_mutations(karyotype, time_interval = time_interval, tau = tau, l = l, mu = mu, w = w, segment_id = j)
    data_single_segment <- add_DP_and_NV(karyotype, data_single_segment, coverage = coverage, purity = purity)

    data_single_segment$tau = tau
    data_single_segment$segment_name = names[j]
    data_all_segments <- dplyr::bind_rows(data_all_segments, data_single_segment)
  }

  peaks_all_segments <- matrix(0, nrow = S, ncol = 2)
  for (i in 1:S) {
    peaks_all_segments[i,] <- get_clonal_peaks(karyotype, purity)
  }

  return(data_all_segments)
}





#' simulate_mutations Function
#'
#' This function allows you to obtain sample mutations to use for CN Timing inference.
#' @param karyotype karyotype
#' @param time_interval time
#' @param tau tau generated
#' @param l length of the segment considered
#' @param mu mutation rate
#' @param w  cell division rate
#' @param segment_id segment id
#' @keywords mutations
#' @export
#' @examples
#' simulate_mutations()


simulate_mutations = function(karyotype, time_interval, tau, l, mu, w, segment_id = "segment_id") {
  available_karyotypes = c('2:0', '2:1', '2:2')
  # Check input
  if (!(time_interval > 0)) stop("time_interval must be a positive real number!")
  if (!(dplyr::between(tau, 0, 1))) stop("tau must be in [0,1]!")
  if (!(l > 0)) stop("l, which is the length of the segment considered, must be a positive real number!")
  if (!(mu > 0)) stop("mu, the mutation rate, must be a positive real number!")
  if (!(w > 0)) stop("w, the cell division rate, must be a positive real number!")
  if (!(karyotype %in% available_karyotypes)) stop("input karyotype is not available!")

  # Create time intervals
  t1 = time_interval * tau
  t2 = time_interval - t1

  # Extract nA and nB
  multiplicities = strsplit(karyotype, ":") %>% unlist() %>% as.numeric()
  nA = multiplicities[1]
  nB = multiplicities[2]

  # Init vector of used mutations
  used = c()
  all_mutations_df <- data.frame(mutation = c(), allele = c(), type = c())

  generate_mutations <- function(allele, l, mu, w, dt, type) {

    rate <- l * mu * w * dt
    n_mutations <- rpois(1, lambda = rate)

    if (n_mutations >= 1) {
      mutations <- lapply(c(1:n_mutations), function(i) {
        repeat{
          id = sample(LETTERS, 8, replace = TRUE) %>% paste(collapse = "")
          if(!id %in% used ) break
        }
        used <<- c(used, id) # <<- changes global value
        id
      })

      mutations_df <- data.frame(
        mutation = mutations %>% unlist(),
        allele = allele,
        type = type
      )
      return(mutations_df)
    } else {
      return(data.frame(mutation = c(), allele = c(), type = c()))
    }
  }

  # Simulate mutations on A
  if (nA == 2) {
    m1 <- generate_mutations("A1A2", l, mu, w, t1, "shared")
    m2 <- generate_mutations("A1", l, mu, w, t2, "private")
    m3 <- generate_mutations("A2", l, mu, w, t2, "private")

    all_mutations_df <- bind_rows(all_mutations_df, m1, m2, m3)
  }

  # Simulate mutations on B
  if (nB == 2) {
    m1 <- generate_mutations("B1B2", l, mu, w, t1, "shared")
    m2 <- generate_mutations("B1", l, mu, w, t2, "private")
    m3 <- generate_mutations("B2", l, mu, w, t2, "private")

    all_mutations_df <- bind_rows(all_mutations_df, m1, m2, m3)
  } else if (nB == 1) {
    m1 <- generate_mutations("B1", l, mu, w, time_interval, "private")
    all_mutations_df <- bind_rows(all_mutations_df, m1)
  }

  if (!(nrow(all_mutations_df) > 0)) stop("no mutations were simulated with the given parameters!")
  all_mutations_df <- all_mutations_df %>% mutate(karyotype = karyotype, segment_id = segment_id, , segment_name_real = segment_id)  # understand how to handle it in a better way
  all_mutations_df
}



#' add_DP_and_NV Function
#'
#' This function allows you to simulate the NV and DP for the CN events.
#' @param karyotype karyotype
#' @param simulated_mutations mutations
#' @param coverage coverage
#' @param purity purity
#' @keywords NV
#' @export
#' @examples
#' simulate_mutations()


add_DP_and_NV = function(karyotype, simulated_mutations, coverage, purity) {
  multiplicities = strsplit(karyotype, ":") %>% unlist() %>% as.numeric()
  nA = multiplicities[1]
  nB = multiplicities[2]

  peak <- function(m, n, purity) {
    num <- m * purity
    den <- (1 - purity) * 2 + n * purity
    return(num/den)
  }
  # Simulate DP and NV values for shared mutation
  p <- peak(nA, nA + nB, purity)
  shared_mutations <- simulated_mutations %>% filter(type == "shared")

  n <- nrow(shared_mutations)
  shared_mutations$DP <- rpois(n, coverage)
  shared_mutations$NV <- rbinom(n, size = shared_mutations$DP, prob = p)

  # Simulate DP and NV values for private mutation
  p <- peak(1, nA + nB, purity)
  private_mutations <- simulated_mutations %>% filter(type == "private")

  n <- nrow(private_mutations)
  private_mutations$DP <- rpois(n, coverage)
  private_mutations$NV <- rbinom(n, size = private_mutations$DP, prob = p)

  return(bind_rows(private_mutations, shared_mutations))
}





get_major = function (k, get_minor = FALSE){
  Major = 0
  minor = 0
  
  k1 = stringr::str_split_1(k, pattern = ":")[1]
  k2 = stringr::str_split_1(k, pattern = ":")[2]
  
  if (k1 > k2){
    Major = k1
    minor = k2} else{
      Major = k2
      minor = k1}
  major_minor = list(Major=Major, minor = minor)
  
  if (get_minor == TRUE){
    return(minor)
  }else{
    return(Major)}
}



get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)
  
  # get only Major and 1
  multiplicities <- c(1, major)
  
  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
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




convert_to_vcf <- function(mutations, cna, metadata) {
  # Validate input
  if(missing(mutations) || missing(cna) || missing(metadata)) {
    stop("mutations, cna, and metadata are required")
  }
  
  # Create VRanges object
  v <- VRanges(
    seqnames = unlist(lapply(mutations$chr, function(c) {stringr::str_replace_all(c, "chr", "")})),
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
