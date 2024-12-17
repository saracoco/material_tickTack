source("../scripts/simulate_functions.R") # to sample tau and karyo

get_simulation_MutationTimeR = function(number_clocks, number_events, purity, coverage, epsilon, seed = 123, vector_karyo = c("2:0", "2:1", "2:2")) {
  
  
  # prepare simulation
  set.seed(seed=seed)
  
  
  weights_karyo <- rep(1/length(vector_karyo), length(vector_karyo))
  
  vector_tau = runif(number_clocks, 0, 1)
  while(min(diff(vector_tau)) <= epsilon) {
    vector_tau = runif(number_clocks, 0, 1)
  }
  
  weights_tau <- rep(1/number_clocks, number_clocks)
  data_simulation <- get_taus_karyo_gerstrung(number_events, vector_tau, vector_karyo, weights_tau, weights_karyo)
  
  
  
  # simulation from mutationtimer
  cn <- MutationTimeR:::refLengths[c(1:23)]
  
  clusters <- data.frame(cluster=1:1, proportion=c(purity), n_ssms=c(100))
  clusters <- data.frame(cluster=1:1, proportion=c(purity), n_ssms=c(100))
  
  time2pi <- function(N,n,t1,t2){
    if(N==2 & n ==1)
      pi <-  c(3-2*t1, t1)
    else if(N==2 & n %in% c(0,2))
      pi <- c(2 -2*t1, t1)
    else pi <- 1
    pi <- pi/sum(pi)
  }
  
  
  ######## my modifications to simulation #############
  major_cn <- lapply(data_simulation$karyo, function(x) as.numeric(strsplit(x,split=":")[[1]][1]))%>%unlist
  minor_cn <- lapply(data_simulation$karyo, function(x) as.numeric(strsplit(x,split=":")[[1]][2]))%>%unlist
  
  df <- bind_cols(major_cn=major_cn,minor_cn=minor_cn, chr = data_simulation$chr, taus = data_simulation$taus)
  
  all_chr <- tibble(chr = 1:23)
  
  df <- all_chr %>%
    left_join(df, by = "chr") %>%
    mutate(
      major_cn = if_else(is.na(major_cn), 1, major_cn),
      minor_cn = if_else(is.na(minor_cn), 1, minor_cn)
    )
  
  cn$major_cn <- df$major_cn
  cn$minor_cn <- df$minor_cn
  
  cn$clonal_frequency <- purity
  ##############
  
  
  cn$timing_param <- MutationTimeR:::defineMcnStates(cn, purity=purity, clusters=clusters, deltaFreq = 0)
  
  
  for(i in seq_along(cn)){
    t1 <- df$taus[i]
    pi <- time2pi(cn$major_cn[i], cn$minor_cn[i], t1, t2)
    pi_sub <- clusters$n_ssms
    pi_sub[1] <- pi_sub[1] * (t1*2 + (1-t1)*(cn$major_cn[i]+ cn$minor_cn[i])) / 2
    pi_sub[-1] <- pi_sub[-1] * (cn$major_cn[i]+ cn$minor_cn[i]) / 2
    pi_sub <- pi_sub/sum(pi_sub)
    pi <- c(pi * pi_sub[1], pi_sub[-1])
    cn$timing_param[[i]][,"P.m.sX"] <- pi
    cn$timing_param[[i]][,"power.m.s"] <- rep(1, length(pi))
  }
  
  
  
  # capire parametri
  cn$n.snv_mnv <- width(MutationTimeR:::refLengths[1:23])/1e6 * 10
  
  vcf <- MutationTimeR:::simulateMutations(cn, rho=0.01, n=40)
  cn_timed <- cn[,c("major_cn","minor_cn","clonal_frequency")]
  
  res = list( df = df, vcf=vcf, cn = cn, cn_timed=cn_timed, clusters = clusters, purity = purity)
  return(res)
}

