# setwd(~/Documents/GitHub/material_tickTack/PCAWG)
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("utils.R")

# sshfs Orfeo:/orfeo/scratch/cdslab dati_Orfeo
results_path <- "~/dati_Orfeo/material_tickTack/PCAWG/results_whole/"
ttypes <- read.delim("data/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv", sep = "\t") %>% 
  #select(sample_id, ttype) %>% 
  dplyr::distinct()

IDs <- list.files(results_path)
RES = readRDS("results/summary_all_samples.rds")
# load("data/gene_coordinates_hg19.rda")
chromosomes_table = CNAqc:::get_reference('hg19') %>%
  mutate(p_arm_length = centromerStart - from, q_arm_length = to - centromerEnd)
# Classes: whole chromosome, arm-event (p-q), focal event
classify_segment = function(s, data){
  #print(segment)
  segment = data[s,]
  print(segment)
  if (segment$chr %in% chromosomes_table$chr){
    targeted_chr =  chromosomes_table %>% filter(chr == segment$chr) 
    segment %>% mutate(len= to-from) %>% mutate(
      class = case_when(
        ## Amplifications 
        # whole chromosome : the CNA involves more than 80% of the total length of the chromosome
        len >= targeted_chr$length * .8 ~ 'whole_chromosome',
        # p-arm : The CNA is on the p-arm and involves more than 80% of its length
        (from + targeted_chr$from < targeted_chr$centromerStart ) & (len >= targeted_chr$p_arm_length * .8) ~ 'p_arm',
        # q-arm : The CNA is on the q-arm and involves more than 80% of its length
        (from + targeted_chr$from > targeted_chr$centromerEnd)  & (len >= targeted_chr$q_arm_length * .8) ~ 'q_arm',
        # focal : The CNA is on one arm and involves less than 80% of its length
        .default = 'focal'
      )
    )
  }else{
    segment$class = NA
  }
  
}
annotate_cna = lapply(1:nrow(RES), classify_segment, data= RES)
annotate_cna =Reduce(rbind, annotate_cna)
annotate_cna = annotate_cna %>% filter(class != 'focal') %>% 
  mutate(class=paste0(chr, ':', class))

write.csv(annotate_cna, 'data/annotated_cnas.csv')
# Retrieve ttypes 
# retrieve_ttypes = function(i){
#   print(i)
#   cnaqc_fits = '~/dati_Orfeo/data/clonal_analysis_PCAWG/'
#   fit = readRDS(paste0(cnaqc_fits, i, '/fit.rds'))
#   ttype = strsplit(fit$snvs$project_code,'-')[[1]][1]
#   data.frame('sample_id'=i, 'code'=ttype)
# }
# with_cna_ttypes = lapply(strsplit(annotate_cna$sample_id %>% unique(), '.rds') %>% unlist(), retrieve_ttypes)  
# with_cna_ttypes = Reduce(rbind, with_cna_ttypes)

# annotate_cna = left_join(annotate_cna %>% rowwise() %>% 
#                            mutate(sample_id = strsplit(sample_id, '.rds')[[1]]), 
#                          with_cna_ttypes, by= 'sample_id')

# annotate_cna = annotate_cna %>% filter(class != 'focal') %>% group_by(code, class) %>% mutate(n= n()) %>% 
#   filter(n > 3)
annotate_cna %>% group_by(ttype, class) %>% mutate(n= n()) %>% 
  filter(n > 3) %>% ggplot(aes(x = class)) + 
  geom_bar() + facet_wrap(~ttype, scales = 'free')

repeated_events = annotate_cna %>% pull(class) %>% unique() ##%>% group_by(code, class) %>% mutate(n= n()) %>% pull(class) #%>% filter(n > 3) %>% pull(class)

# score_df = lapply(annotate_cna$code %>% unique(), function(t){
#   print(t)
#   df = annotate_cna %>% filter(code == t)
#   combinations = expand.grid(df$class %>% unique(), df$class %>% unique()) %>% filter(Var1 != Var2)
#   #distinct(combinations)
#   seen = c()
#   tt_df = data.frame('code'=c(), 'event1'=c(), 'event2'=c(), 'cooccurrence'=c(), 'score'=c())
#   # For each combination
#   for (c in 1:nrow(combinations)){
#     combo = c(combinations[c,1], combinations[c,2])
#     
#     if ( !(paste0(combo[1],',',combo[2]) %in% seen)){
#       co_occurrence = 0
#       score = 0
#       seen = c(paste0(combo[1],',',combo[2]), paste0(combo[2],',',combo[1]))
#       # For each sample
#       for (s in (df$sample_id %>% unique())){
#         filtered_df = df %>% filter(sample_id == s, class %in% c(combo))
#         # Se i due eventi sono entrambi presenti 
#         if (length(filtered_df$class %>% unique())==2){
#           co_occurrence = co_occurrence+1
#           clock_1 = filtered_df %>% filter(class == combo[1]) %>% pull(clock_rank)
#           clock_2 = filtered_df %>% filter(class == combo[2]) %>% pull(clock_rank)
#           if (clock_1 > clock_2){score=score+1}
#           if (clock_2 > clock_1){score=score-1}
#         }
#       }
#       tt_df = rbind(
#         tt_df,
#         data.frame('code'=t, 'event1'=combo[1], 'event2'=combo[2], 'cooccurrence'=co_occurrence, 'score'=score)
#       )
#     }
#   
#   }
#   tt_df
# })
#   
# score_df = Reduce(rbind, score_df)
# 
# library(ggsci)
# 
# write.csv(score_df, file=paste0(getwd(), '/annotated_CNAs_PCAWG.csv'))
# score_df = read.csv(paste0(getwd(), '/annotated_CNAs_PCAWG.csv'))
# score_df %>% filter(cooccurrence >= 2, abs(score) > 1) %>% ggplot() +
#   geom_tile(aes(x = event1, y = event2, fill= score))+
#   facet_wrap(~code, scales = 'free') +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle=45, hjust=1)
#   ) +
#   scale_fill_gsea()

# Create all pairs of drivers
driver_pairs <- expand.grid(first_driver = repeated_events, second_driver = repeated_events) %>%
  dplyr::mutate(first_driver = as.character(first_driver), second_driver=as.character(second_driver)) %>% 
  dplyr::filter(first_driver <= second_driver) # Avoid duplicate/reverse pairs

# Precompute necessary data from res_w_drivers
res_processed <- annotate_cna %>%
  group_by(sample_id) %>%
  mutate(driver_in_pair = class %in% repeated_events) %>%
  filter(driver_in_pair) %>%
  ungroup()

# Initialize results as a list for better performance
results <- vector("list", nrow(driver_pairs))

# Calculate scores for each pair
#k = 69
all_scores = data.frame()
for (t in res_processed$ttype %>% unique()){
  res_processed_t = res_processed %>% filter(ttype==t)
  scores_df = lapply(1:nrow(driver_pairs), function(k) {
    print(k)
    pair <- driver_pairs[k, ]
    
    co_occurr_df = res_processed_t %>%
      filter(class %in% c(pair$first_driver, pair$second_driver)) %>%
      group_by(sample_id) %>%
      mutate(n = n()) %>%
      filter(n != 1) %>% 
      dplyr::select(sample_id, clock_rank, class) %>% 
      tidyr::pivot_wider(values_from = clock_rank, names_from = class)
    sample_ids = paste(co_occurr_df$sample_id %>% unique(), collapse=' ')
    if (nrow(co_occurr_df)) {
      n_pre = sum(co_occurr_df[,pair$first_driver] < co_occurr_df[,pair$second_driver])
      n_coo = sum(co_occurr_df[,pair$first_driver] == co_occurr_df[,pair$second_driver])
      n_post = sum(co_occurr_df[,pair$first_driver] > co_occurr_df[,pair$second_driver])
      
      # Simulated data: Replace these counts with your actual observed data
      O_counts <- c(n_pre, n_coo, n_post)  # Observed frequencies for (-1, 0, 1)
      #DescTools::MultinomCI(O_counts)
      
      # is the mean different from zero?
      t_test <- tryCatch({
        t.test(c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post)), mu = 0)
      }, error = function(e) {
        #cat("Error at iteration", i, "- setting default value\n")
        return(list('p.value'=NA))  # Set your desired default value here
      })
      #t_test <- t.test(c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post)), mu = 0)
      t_test$p.value
      
      score = mean(c(rep(-1, n_pre), rep(0, n_coo), rep(1, n_post)))
      n_samples = sum(O_counts)
      
      # Store result
      return(tibble(first_driver = pair$first_driver, 
                    second_driver = pair$second_driver, score = score, n_samples=n_samples, 
                    p.value = t_test$p.value, type = t, samples = sample_ids))
    }
  }) %>% do.call("bind_rows", .)
  all_scores = rbind(all_scores,scores_df)
}

#cores_df

all_scores = all_scores %>% rowwise() %>% 
  mutate(
    first_driver = paste0(strsplit(strsplit(first_driver, 'chr')[[1]][2], ':')[[1]][1],strsplit(strsplit(first_driver, ':')[[1]][2],'_')[[1]][1]),
    second_driver = paste0(strsplit(strsplit(second_driver, 'chr')[[1]][2], ':')[[1]][1],strsplit(strsplit(second_driver, ':')[[1]][2],'_')[[1]][1])
  ) 
write.csv(all_scores, 'data/arm_scores.csv')

all_scores %>%  dplyr::filter(p.value <= .1) %>% 
  na.omit() %>% 
  ggplot(mapping = aes(x=first_driver, y=second_driver, fill=score)) +
  geom_tile() +
  geom_text(aes(label=n_samples)) +
  scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white") +
  theme_bw() +
  labs(x = "CNA 1", y = 'CNA 2') +
  facet_wrap(~type, scales = 'free')+
  theme(
    #axis.text.x = element_text(angle= 45, hjust=1)
  )
#ggtitle(tumour_type)


ggsave(paste0("plot/matrices/arm_events.png"), width = 10, height =10, units="in", dpi=300)



########## Circus plot
library(circlize)
library(dplyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
source("utils.R")
all_scores = read.csv('data/arm_scores.csv') %>% filter(p.value < .1)
annotate_cna = read.csv('data/annotated_cnas.csv') %>% rowwise() %>%
  mutate(class = paste0(strsplit(strsplit(class, 'chr')[[1]][2], ':')[[1]][1],strsplit(strsplit(class, ':')[[1]][2],'_')[[1]][1]))
all_events = c(all_scores$first_driver %>% unique(), all_scores$second_driver %>% unique()) %>% unique() %>% sort()

# Annotate the karyotypes
karyo_annotations = lapply(all_events, function(e){
  samples = c(all_scores %>% filter(first_driver == e) %>% pull(samples),
              all_scores %>% filter(second_driver == e) %>% pull(samples))
  samples = strsplit(samples, ' ') %>% unlist()
  k_df = annotate_cna %>% filter(class == e, sample_id %in% samples) %>% group_by(karyotype) %>% summarise(n=n())
  k_df$event = e
  k_df
})
karyo_annotations = Reduce(rbind, karyo_annotations)

# Annotate the events with drivers
load("data/gene_coordinates_hg19.rda")
drivers = read.csv('./IntOGen_2024/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep='\t') %>%
  filter(CANCER_TYPE %in% c('LIHB','PANCREAS', 'PAAD', 'MEL', 'BRCA', 'COAD','COADREAD')) %>%
  filter(SAMPLES > 10)
gene_coordinates_hg19 = gene_coordinates_hg19 %>% filter(gene %in% drivers$SYMBOL)

genes = lapply(all_events, function(c){
  if (grepl('p',c)){
    chromosome = strsplit(c, 'p')[[1]][1]
    from_e = CNAqc:::get_reference('hg19') %>% 
      filter(chr==paste0('chr',chromosome)) %>% pull(from)
    to_e = CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% 
      pull(centromerStart)
    genes = gene_coordinates_hg19 %>% filter(chr==paste0('chr',chromosome), from >= 0, to <= (to_e-from_e)) 
  }
  if (grepl('q',c)){
    chromosome = strsplit(c, 'q')[[1]][1]
    from_0 = CNAqc:::get_reference('hg19') %>% 
      filter(chr==paste0('chr',chromosome)) %>% pull(from)
    
    from_e = CNAqc:::get_reference('hg19') %>% 
      filter(chr==paste0('chr',chromosome)) %>% pull(centromerEnd)
    to_e = CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% 
      pull(to)
    genes = gene_coordinates_hg19 %>% filter(chr==paste0('chr',chromosome), from >= from_e-from_0, to <= (to_e-from_0)) 
  }
  if (grepl('whole',c)){
    chromosome = strsplit(c, 'whole')[[1]][1]
    from_e = CNAqc:::get_reference('hg19') %>% 
      filter(chr==paste0('chr',chromosome)) %>% pull(from)
    to_e = CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% 
      pull(to)
    genes = gene_coordinates_hg19 %>% filter(chr==paste0('chr',chromosome), from >= from_e, to <= (to_e-from_e)) 
  }
  genes$arm = c
  genes
})
genes = Reduce(rbind, genes)

# By ttype
# for (tt in (all_scores$type %>% unique())){
#   pdf(paste0("plot/circos_plot",tt,".pdf"), width = 8, height = 8)
#   all_scores_tt = all_scores %>% filter(type == tt) %>% filter(p.value < .1)
#   sectors <- c(all_scores_tt$first_driver %>% unique(), all_scores_tt$second_driver %>% unique()) %>% unique() %>% sort()
#   values = lapply(sectors, function(c){
#     if (grepl('p',c)){
#       chromosome = strsplit(c, 'p')[[1]][1]
#       arm='p'
#       len = (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(centromerStart)) - (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(from))
#     }
#     if (grepl('q',c)){
#         chromosome = strsplit(c, 'q')[[1]][1]
#         arm='q'
#         len = (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(to)) - (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(centromerEnd))
#     }
#     if (grepl('whole',c)){
#       chromosome = strsplit(c, 'whole')[[1]][1]
#       arm='whole'
#       len = (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(to)) - (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(from))
#     }
#     data.frame('sector'= c,'value'=len)
#     #CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) 
#   }) #%>% Reduce(rbind)
#   values = Reduce(rbind,values)
#   
#   values$value_scaled <- values$value / sum(values$value)
#   df = values 
#   
#   df$start <- c(0, head(cumsum(df$value_scaled), -1))
#   df$end <- cumsum(df$value_scaled)
#   
#   # Initialize the circular layout with proportional widths
#   circos.initialize(factors = df$sector, xlim = df[, c("start", "end")],
#                     circos.par(cell.padding = c(0.02, 0, 0.02, 0)))
#   
#   # Create sector tracks with labels
#   circos.trackPlotRegion(factors = df$sector, ylim = c(0, 1), 
#                          panel.fun = function(x, y) {
#                            sector.name <- get.cell.meta.data("sector.index")
#                            xcenter <- mean(get.cell.meta.data("xlim"))
#                            circos.text(xcenter, 0.5, sector.name, facing = "inside", niceFacing = TRUE)
#                          },
#                          bg.border = rep("white", nrow(df)),
#                          bg.col = rep("gainsboro", nrow(df)))
#   for (i in 1:nrow(all_scores_tt)){
#     first = all_scores_tt$first_driver[i]
#     second = all_scores_tt$second_driver[i]
#     if (all_scores_tt$score[i] < 0){
#       transparent_violet <- rgb(148/255, 131/255, 204/255, alpha = 0.5)
#       circos.link(first, 
#                   c(df %>% filter(sector == first) %>% pull(start), df %>% filter(sector == first) %>% pull(end)), 
#                   second, 
#                   c(df %>% filter(sector == second) %>% pull(start), df %>% filter(sector == second) %>% pull(end)), 
#                   col = transparent_violet, border = transparent_violet#, directional = 1, arr.length = 0.4
#       )
#     }else{
#       transparent_orange <- rgb(255/255, 165/255, 0, alpha = 0.5)
#       circos.link(first, 
#                   c(df %>% filter(sector == first) %>% pull(start), df %>% filter(sector == first) %>% pull(end)), 
#                   second, 
#                   c(df %>% filter(sector == second) %>% pull(start), df %>% filter(sector == second) %>% pull(end)), 
#                   col = transparent_orange, border = transparent_orange#, directional = 1, arr.length = 0.4
#       )
#     }
#     
#     }
#   
#   # Clear the Circos plot
#   circos.clear()
#   dev.off()
# }

# All types
pdf("plot/circos_plot.pdf", width = 8, height = 8)
sectors <- c(all_scores$first_driver %>% unique(), all_scores$second_driver %>% unique()) %>% unique() %>% sort()
values = lapply(sectors, function(c){
  if (grepl('p',c)){
    chromosome = strsplit(c, 'p')[[1]][1]
    arm='p'
    len = (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(centromerStart)) - (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(from))
  }
  if (grepl('q',c)){
    chromosome = strsplit(c, 'q')[[1]][1]
    arm='q'
    len = (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(to)) - (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(centromerEnd))
  }
  if (grepl('whole',c)){
    chromosome = strsplit(c, 'whole')[[1]][1]
    arm='whole'
    len = (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(to)) - (CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) %>% pull(from))
  }
  data.frame('sector'= c,'value'=len)
  #CNAqc:::get_reference('hg19') %>% filter(chr==paste0('chr',chromosome)) 
}) #%>% Reduce(rbind)
values = Reduce(rbind,values)

values$value_scaled <- values$value / sum(values$value)
df = values 
df$start <- c(0, head(cumsum(df$value_scaled), -1))
df$end <- cumsum(df$value_scaled)
# Initialize the circular layout with proportional widths
# circos.initialize(factors = df$sector, xlim = df[, c("start", "end")],
#                     circos.par(cell.padding = c(0.02, 0, 0.02, 0)))
circos.initialize(factors = df$sector, xlim = df[, c("start", "end")],
                  circos.par(start.degree = 90, gap.degree = 15))

genes_scaled = lapply(genes$gene, function(g){
  gene_info <- genes %>% filter(gene==g)
  sector_name <- gene_info$arm
  sector_data <- df %>% filter(sector == sector_name)
  gene_scaled <- sector_data$start + 
    (gene_info$from / max(values$value)) * (sector_data$end - sector_data$start)
  gene_info$scaled = gene_scaled
  gene_info
})
genes_scaled = Reduce(rbind, genes_scaled)
circos.labels(genes_scaled$arm, 
              x = genes_scaled$scaled, 
              labels = genes_scaled$gene,
              cex = 0.7,
              side = "outside")

# Create sector tracks with labels
circos.trackPlotRegion(factors = df$sector, ylim = c(0, 1), 
                       panel.fun = function(x, y) {
                         sector.name <- get.cell.meta.data("sector.index")
                         xcenter <- mean(get.cell.meta.data("xlim"))
                         circos.text(xcenter, 0.3, sector.name, facing = "inside", niceFacing = TRUE)
                       },
                       bg.border = rep("grey94", nrow(df)),
                       bg.col = rep("grey94", nrow(df)))

# for (e in all_events){
#   
#   circos.track(factors = e, ylim = c(0, 1), panel.fun = function(x, y) {
#     value = karyo_annotations %>% filter(event == e) %>% pull(n)
#     position1 = df %>% filter(sector == e) %>% pull(start)
#     position2 = mean(df %>% filter(sector == e) %>% pull(start),
#                      df %>% filter(sector == e) %>% pull(end))
#     position3 = df %>% filter(sector == e) %>% pull(end)
#     circos.barplot(value, pos = c(position1,position2,position3), 
#                    col =karyo_annotations$karyotype,
#                    sector.index = karyo_annotations$event)
#   })
#   
# }
#circos.track(ylim = c(0, 1))

# for (i in 1:nrow(genes)) {
#   gene_info <- genes[i, ]
#   sector_name <- gene_info$arm  # Match arm names with Circos sectors
#   
#   # Get the start position of the sector
#   sector_data <- df %>% filter(sector == sector_name)
#   
#   if (nrow(sector_data) > 0) {
#     # Scale gene position relative to the sector
#     gene_scaled <- sector_data$start + 
#       (gene_info$from / max(values$value)) * (sector_data$end - sector_data$start)
#     
#     # Define y-positions
#     y_start <- 1.05  # Start of annotation line
#     y_end <- 1.2     # Where gene label appears
#     
#     # Draw outward annotation line
#     circos.lines(x = rep(gene_scaled, 2), y = c(1.5, y_start), sector.index = sector_name)
#     
#     # Add gene name, facing outward
#     circos.text(x = gene_scaled, y = y_end, labels = gene_info$gene, 
#                 sector.index = sector_name, facing = "clockwise", 
#                 cex = .7, adj = c(-.55, 0), niceFacing = T)
#   }
# }


fill = 
  c("LIRI"=rgb(148/255, 131/255, 204/255, alpha = 0.5),
    "PACA"=rgb(247/255, 216/255, 133/255, alpha = .5),
    "ESAD"=rgb(136/255, 181/255, 215/255, alpha = .5),
    "MALY"=rgb(144/255, 252/255, 206/255, alpha = .5),
    "MELA"=rgb(255/255, 255/255, 167/255, alpha = 0.5),
    "BRCA"=rgb(246/255, 196/255, 205/255, alpha = 0.5),
    "BOCA"=rgb(203/255, 204/255, 250/255, alpha = 0.5))
color = 
  c("LIRI"=rgb(148/255, 131/255, 204/255, alpha = 1),
    "PACA"=rgb(247/255, 216/255, 133/255, alpha = 1),
    "ESAD"=rgb(136/255, 181/255, 215/255, alpha = 1),
    "MALY"=rgb(144/255, 252/255, 206/255, alpha = 1),
    "MELA"=rgb(255/255, 255/255, 167/255, alpha = 1),
    "BRCA"=rgb(246/255, 196/255, 205/255, alpha = 1),
    "BOCA"=rgb(203/255, 204/255, 250/255, alpha = 1))

legend(1,1, legend=names(color), fill=color)

for (i in 1:nrow(all_scores)){
  first = all_scores$first_driver[i]
  second = all_scores$second_driver[i]
  t = all_scores$type[i]
  if (all_scores$score[i] < 0){
    # Ribbon
    circos.link(first,
                c(df %>% filter(sector == first) %>% pull(start), df %>% filter(sector == first) %>% pull(end)),
                second,
                c(df %>% filter(sector == second) %>% pull(start), df %>% filter(sector == second) %>% pull(end)),
                col = fill[t], border = color[t], directional = 1, arr.length = 0.1,arr.type = 'big.arrow')
    # Arrow
    # circos.link(first, #median(df %>% filter(sector == first) %>% pull(start), df %>% filter(sector == first) %>% pull(end)), 
    #             second, 
    #             lwd = 2,h=0.5,l=0.5,#median(df %>% filter(sector == second) %>% pull(start), df %>% filter(sector == second) %>% pull(end)), 
    #             col = color[t], border = color[t], directional = 1, arr.length = 0.5)
  }else{
    # Ribbon
    circos.link(first,
                c(df %>% filter(sector == first) %>% pull(start), df %>% filter(sector == first) %>% pull(end)),
                second,
                c(df %>% filter(sector == second) %>% pull(start), df %>% filter(sector == second) %>% pull(end)),
                col = fill[t], border = color[t], directional = -1, arr.length = 0.1,arr.type = 'big.arrow')
    # Arrow
    # circos.link(first, #median(df %>% filter(sector == first) %>% pull(start), df %>% filter(sector == first) %>% pull(end)), 
    #             second, 
    #             lwd = 2,h=0.5,l=0.5, #median(df %>% filter(sector == second) %>% pull(start), df %>% filter(sector == second) %>% pull(end)), 
    #             col = color[t], border = color[t], directional = -1, arr.length = 0.5)
  }
  
  
  
}

# Clear the Circos plot
circos.clear()
dev.off()




