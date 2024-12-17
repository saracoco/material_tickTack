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