res <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/segments_analysis/data/info_segments_1e+06.rds")
unique(res$type)
res_BRCA <- res %>% filter(type=="BRCA")
saveRDS(res_BRCA, "./data/info_segments_BRCA_1e6.rds")
