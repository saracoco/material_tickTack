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
# Retrieve ttypes 
retrieve_ttypes = function(i){
  print(i)
  cnaqc_fits = '~/dati_Orfeo/data/clonal_analysis_PCAWG/'
  fit = readRDS(paste0(cnaqc_fits, i, '/fit.rds'))
  ttype = strsplit(fit$snvs$project_code,'-')[[1]][1]
  data.frame('sample_id'=i, 'code'=ttype)
}
with_cna_ttypes = lapply(strsplit(annotate_cna$sample_id %>% unique(), '.rds') %>% unlist(), retrieve_ttypes)  
with_cna_ttypes = Reduce(rbind, with_cna_ttypes)

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
                    p.value = t_test$p.value, type = t))
    }
  }) %>% do.call("bind_rows", .)
  all_scores = rbind(all_scores,scores_df)
}

#cores_df

all_scores %>% 
  dplyr::filter(p.value <= .1) %>% 
  na.omit() %>% 
  ggplot(mapping = aes(x=first_driver, y=second_driver, fill=score)) +
  geom_tile() +
  geom_text(aes(label=n_samples)) +
  scale_fill_gradient2(low = "#998ec3", high = "#f1a340", mid = "white") +
  theme_bw() +
  labs(x = "CNA 1", y = 'CNA 2') +
  facet_wrap(~type, scales = 'free')+
  theme(
    axis.text.x = element_text(angle= 45, hjust=1)
  )
  #ggtitle(tumour_type)


ggsave(paste0("plot/matrices/arm_events.png"), width = 10, height =10, units="in", dpi=300)










  