.libPaths("~/R/rstudio_v3/")
library(AmplificationTimeR)
library(dpclust3p)





# our simulation for AmplificationTimer
sim <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/simulations_v2/results/sim_1_5_0.1_5_10/1/sim.rds")



muts <- sim$muts
cn <- sim$cn 
purity <- 0.1
ploidy = 2



if(file.exists("./pseudo_dpclust_files")){
} else {
  dir.create("./pseudo_dpclust_files")
}
if(file.exists("./pseudo_dpclust_files/pseudo_allele_frequency/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_allele_frequency/")
}
if(file.exists("./pseudo_dpclust_files/pseudo_cellularity/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_cellularity/")
}
if(file.exists("./pseudo_dpclust_files/pseudo_loci_file/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_loci_file/")
}
if(file.exists("./pseudo_dpclust_files/pseudo_subclones/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_subclones/")
}
if(file.exists("./dpclust_info")){
} else {
  dir.create("./dpclust_info")
}

write.table(cn, "./pseudo_dpclust_files/pseudo_subclones/sample_1_subclones.txt", col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")


# First, initialize counts to 0
muts$Count_A <- 0
muts$Count_C <- 0
muts$Count_G <- 0
muts$Count_T <- 0

# Assign reference allele counts
muts$Count_A[muts$ref == "A"] <- muts$DP[muts$ref == "A"] - muts$NV[muts$ref == "A"]
muts$Count_C[muts$ref == "C"] <- muts$DP[muts$ref == "C"] - muts$NV[muts$ref == "C"]
muts$Count_G[muts$ref == "G"] <- muts$DP[muts$ref == "G"] - muts$NV[muts$ref == "G"]
muts$Count_T[muts$ref == "T"] <- muts$DP[muts$ref == "T"] - muts$NV[muts$ref == "T"]

# Assign alt allele counts
muts$Count_A[muts$alt == "A"] <- muts$NV[muts$alt == "A"]
muts$Count_C[muts$alt == "C"] <- muts$NV[muts$alt == "C"]
muts$Count_G[muts$alt == "G"] <- muts$NV[muts$alt == "G"]
muts$Count_T[muts$alt == "T"] <- muts$NV[muts$alt == "T"]

# Total depth (sanity check)
muts$Total_depth <- muts$Count_A + muts$Count_C + muts$Count_G + muts$Count_T

muts_alleles <- (muts[,c("chr","start","Count_A","Count_C","Count_G","Count_T","Total_depth")])
muts_alleles[is.na(muts_alleles)] <- 0
head(muts_alleles)


colnames(muts_alleles) <- c("CHR","START","Count_A","Count_C","Count_G","Count_T","Total_depth")
write.table(muts_alleles, "./pseudo_dpclust_files/pseudo_allele_frequency/sample_1_pseudo_allele_frequency.tsv", col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")

muts_loci <- muts[,c("chr","start","ref","alt")]
colnames(muts_loci) <- c("CHR","START","REF","ALT")
write.table(muts_loci, "./pseudo_dpclust_files/pseudo_loci_file/sample_1_pseudo_loci.txt", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")



sample_1_pseudo_cellularity <- matrix(nrow = 3, ncol = 4, dimnames = list(c("ASCAT","FRAC_GENOME","REF_SEG"),c("rho","psi","distance","is.best")),
                                      data = c(NA, NA, NA, NA, purity, ploidy, NA, "TRUE", NA, NA, 0, "FALSE"), byrow = TRUE)

write.table(sample_1_pseudo_cellularity, "./pseudo_dpclust_files/pseudo_cellularity/sample_1_pseudo_cellularity.tsv", col.names = TRUE, quote = FALSE, row.names = TRUE, sep = "\t")



runGetDirichletProcessInfo(loci_file = "./pseudo_dpclust_files/pseudo_loci_file/sample_1_pseudo_loci.txt",
                           allele_frequencies_file = "./pseudo_dpclust_files/pseudo_allele_frequency/sample_1_pseudo_allele_frequency.tsv",
                           cellularity_file = "./pseudo_dpclust_files/pseudo_cellularity/sample_1_pseudo_cellularity.tsv",
                           subclone_file = "./pseudo_dpclust_files/pseudo_subclones/sample_1_subclones.txt",
                           gender = "male",
                           SNP.phase.file = NA, mut.phase.file = NA,
                           output_file = "./dpclust_info/sample_1_dpclust_info")

sample_1_dpclust_info <- read.delim("./dpclust_info/sample_1_dpclust_info", sep = "\t", header = TRUE)
head(sample_1_dpclust_info)


