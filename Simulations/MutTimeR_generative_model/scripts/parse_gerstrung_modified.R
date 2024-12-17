parse_gerstrung_modified = function(vcf, cna) {
  # Mutations
  NV <- vcf@info$t_alt_count
  DP <- vcf@info$t_alt_count + vcf@info$t_ref_count
  from <- vcf@rowRanges@ranges@start
  to <- vcf@rowRanges@ranges@start + vcf@rowRanges@ranges@width - 1
  chr <- paste0("chr", vcf@rowRanges@seqnames)
  
  mutations <- dplyr::tibble(chr = chr, from = from, to = to, NV = NV, DP = DP)
  
  # Segments
  chr <- paste0("chr", cna@seqnames %>% as.numeric())
  from <- cna@ranges@start
  to <- cna@ranges@start + cna@ranges@width - 1
  
  Major <- cna$major_cn
  minor <- cna$minor_cn
  purity <- cna$clonal_frequency
  
  segments <- dplyr::tibble(chr = chr, from = from, to = to, Major = Major, minor = minor, purity = purity)
  
  return(list(mutations = mutations, cna = segments, metadata = tibble(purity = unique(purity)) ))
}
