results_races <- function(phylo_forest){
  
  library(rRACES)
  library(dplyr)
  
  # phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")
  
  sample_name <- phylo_forest$get_samples_info()[["name"]]


  # Sequencing results from Sample A
  Sample_A <-phylo_forest$get_bulk_allelic_fragmentation("Sample A")



  CNA_Clone1_chr_5 <- Sample_A %>% filter(begin==107707518)
  print(CNA_Clone1_chr_5)
  CNA_Clone1_chr_12 <- Sample_A %>% filter(begin==25205246)
  print(CNA_Clone1_chr_12)


  CNA_Clone2_chr_6 <- Sample_A %>% filter(begin==25100000)
  print(CNA_Clone2_chr_6)
  CNA_Clone2_chr_13 <- Sample_A %>% filter(begin==39500001)
  print(CNA_Clone2_chr_13)


  CNA_Clone3_chr_10 <- Sample_A %>% filter(begin==87862638)
  print(CNA_Clone3_chr_10)

  Sample_A_results <- bind_rows(CNA_Clone1_chr_5, CNA_Clone1_chr_12, CNA_Clone2_chr_6, CNA_Clone3_chr_10)

  # # Sequencing results from Sample B
  # Sample_B <-phylo_forest$get_bulk_allelic_fragmentation("Sample B")
  # 
  # 
  # CNA_Clone1_chr_5 <- Sample_B %>% filter(begin==107707518)
  # print(CNA_Clone1_chr_5)
  # CNA_Clone1_chr_12 <- Sample_B %>% filter(begin==25205246)
  # print(CNA_Clone1_chr_12)
  # 
  # 
  # CNA_Clone2_chr_6 <- Sample_B %>% filter(begin==25100000)
  # print(CNA_Clone2_chr_6)
  # CNA_Clone2_chr_13 <- Sample_B %>% filter(begin==39500001)
  # print(CNA_Clone2_chr_13)
  # 
  # 
  # CNA_Clone3_chr_10 <- Sample_B %>% filter(begin==87862638)
  # print(CNA_Clone3_chr_10)
  # 
  # CNA_Clone4_chr_22 <- Sample_B %>% filter(begin==20303470)
  # print(CNA_Clone4_chr_22)
  # 
  # Sample_B_results <- bind_rows(CNA_Clone1_chr_5, CNA_Clone1_chr_12, CNA_Clone2_chr_6, CNA_Clone2_chr_13, CNA_Clone3_chr_10, CNA_Clone4_chr_22)
  # 
  return(Sample_A_results)
}