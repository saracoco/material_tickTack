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
