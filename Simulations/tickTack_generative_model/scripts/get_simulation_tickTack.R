source("../scripts/simulate_functions.R") 

get_simulation_tickTack = function(number_clocks, number_events, purity, coverage, epsilon, seed = 123, vector_karyo = c("2:0", "2:1", "2:2")) {
  
  # prepare simulation
  set.seed(seed=seed)
  
  weights_karyo <- rep(1/length(vector_karyo), length(vector_karyo))
  
  vector_tau = runif(number_clocks, 0, 1)
  while(min(diff(vector_tau)) <= epsilon) {
    vector_tau = runif(number_clocks, 0, 1)
  }
  
  weights_tau <- rep(1/number_clocks, number_clocks)
  data_simulation <- get_taus_karyo(number_events, vector_tau, vector_karyo, weights_tau, weights_karyo)
  
  taus <- data_simulation$taus
  karyo = data_simulation$karyo
  
  
  
  # get simulation parametes
  coverage = 100 # average number of reads that align to a reference base
  cat("coverage: ", coverage)
  
  mu = 1e-4 # mutation rate
  cat("mutation rate: ", mu)
  
  w = 1e-2 # cell division rate
  cat("cell division rate: ", w)
  
  l = 2e7 # length of the segment
  cat("length of the segment: ", l)
  
  time_interval = 7
  cat("time_interval: ", time_interval)
  
  
  options(bitmapType='cairo')
  
  
  simulation_data_all_segments = get_simulation(data_simulation$taus, data_simulation$karyo, purity, time_interval, l, mu, w, coverage) # the other parameters have default value assigned if none is specified
  data_simulation_mutations <- simulation_data_all_segments  
  
  
  x <- list( mutations = tibble(chr = data_simulation_mutations$segment_id, from = 2, to = l-1, ref = "", alt = "", DP = data_simulation_mutations$DP, NV = data_simulation_mutations$NV, VAF = data_simulation_mutations$NV/data_simulation_mutations$DP, sample = 1), 
             cna = tibble(chr = unique(data_simulation_mutations$segment_id), from = 1, to = l, Major = as.numeric(lapply(data_simulation$karyo, get_major) %>% unlist()) , minor = as.numeric(lapply(data_simulation$karyo, get_major, get_minor = TRUE) %>% unlist()),   CCF = 0, total_cn = Major + minor), 
             metadata = tibble(purity = purity, sample = "sample 1") 
  )
  
  res = list (data_simulation = data_simulation, x = x)
  
  return(res)
}

