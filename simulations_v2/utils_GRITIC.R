library(ggplot2)

plot_CNAqc = function (data, ref = 'hg19'){
  input_data_CNAqc <- CNAqc::init(mutations=data$mutations, 
                                  cna = data$cna %>% mutate(Major_2 = NA, minor_2 = NA), 
                                  purity = data$metadata$purity, 
                                  ref = ref)
  chrom = c(paste0('chr', 1:22),'chrX', 'chrY')
  shifted_dataset_input = data.frame()
  for (c in chrom){
    chr_from = CNAqc:::get_reference(ref) %>% filter(chr==c) %>% pull(from)
    my_df = data$cna %>% filter(chr==c) %>% 
      mutate(from_abs = from+chr_from, to_abs=to+chr_from)
    shifted_dataset_input = rbind(shifted_dataset_input, my_df)
  }
  
  p <- CNAqc:::blank_genome() +
       geom_segment(data = shifted_dataset_input, 
                 aes(x = from_abs, xend=to_abs, y=Major, yend=Major))

  return(p)}



# this function write the tsv files that are needed by the GRITIC function
convert_CNAqc_to_GRITIC = function(data){
  mutation_tsv = data$mutations %>% dplyr::mutate(Chromosome=chr, Position = as.integer(from), Tumor_Ref_Count=as.integer(DP), Tumor_Alt_Count=as.integer(NV)) 
  mutation_tsv = mutation_tsv %>% dplyr::select(Tumor_Ref_Count,Tumor_Alt_Count,Position,Chromosome)
  # write_tsv(mutation_tsv, "mutation_tsv.tsv")
  
  cna_tsv = data$cna %>% dplyr::mutate(Chromosome=chr, Segment_Start = as.integer(from), Segment_End = as.integer(to), Major_CN=as.integer(Major), Minor_CN=as.integer(minor)) 
  cna_tsv = cna_tsv %>% dplyr::select(Segment_Start,Segment_End,Major_CN,Minor_CN,Chromosome)
  cna_tsv <- cna_tsv %>%
    mutate(
      Segment_Start = as.integer(Segment_Start),
      Segment_End = as.integer(Segment_End),
      Major_CN = as.integer(Major_CN),
      Minor_CN = as.integer(Minor_CN)
    )   %>% dplyr::filter(!is.na(Major_CN), !is.na(Minor_CN))
  # write_tsv(cna_tsv, "cna_tsv.tsv")
  return(list(mut=mutation_tsv, cn=cna_tsv))
}



# this function mimic what happens in gritic preprocess step in the data loading procedure
merge_cn_segments <- function(cn_table) {
  cn_table <- cn_table %>% mutate(Gain_Type = paste0(Major_CN, "_", Minor_CN))
  setDT(cn_table)
  cn_table[, Row_Index := .I]
  
  indexes_to_delete <- c()
  for (chromosome in unique(cn_table$Chromosome)) {
    chr_data <- cn_table[Chromosome == chromosome]
    
    for (i in 1:(nrow(chr_data) - 1)) {
      index <- chr_data[i, Row_Index]
      forward_index <- chr_data[i + 1, Row_Index]
      
      if (cn_table[index, Gain_Type] == cn_table[forward_index, Gain_Type]) {
        indexes_to_delete <- c(indexes_to_delete, index)
        cn_table[forward_index, Segment_Start := cn_table[index, Segment_Start]]
      }
    }
  }
  merged_table <- cn_table[!Row_Index %in% indexes_to_delete]
  merged_table[, Row_Index := NULL]
  merged_table <- as.data.frame(merged_table)
  rownames(merged_table) <- NULL
  
  return(merged_table)
}

fit_gritic <- function(mutation_path="./mutation_tsv.tsv",copy_number_path = "./cna_tsv.tsv",purity = pi,sample_id = "S123",
                       output_dir = "results_gritic_example_PCAWG/",
                       wgd_status = TRUE,
                       penalty = FALSE,
                       subclone_path = NULL,
                       sex = "XY",
                       merge_cn=FALSE,
                       plot_trees = TRUE){
  example_gritic <- reticulate::import_from_path("run_gritic_module", path = ".")
  example_gritic$run_gritic(
    mutation_path = paste0(sub_dir,"/mutation_tsv.tsv"),
    copy_number_path = paste0(sub_dir,"/cna_tsv.tsv"),
    purity = pi,
    sample_id = "S123",
    output_dir = paste0(sub_dir,"/results_gritic/"),
    wgd_status = TRUE,
    penalty = FALSE,
    subclone_path = NULL,
    sex = "XY",
    merge_cn=FALSE,
    plot_trees = TRUE,
    reads_correction=FALSE
  )
  
}
