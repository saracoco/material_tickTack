simulate_tissue <- function( name_sim = "SPN01", width_sim=1e3, height_sim=1e3, seed = 123, death_activation_level = 10, history_delta = 1,
                             clones_number = 3, K_clocks, n_events, taus, simulation_time = 500){
  
  print(taus)
  clones_number = K_clocks
  
  
  # Prep simulation ####
  sim <- SpatialSimulation(name_sim, width=width_sim, height=height_sim, seed = seed)
  sim$history_delta <- history_delta
  sim$death_activation_level <- death_activation_level
  #sim <- new(Simulation, seed = seed, save_snapshot = F)
  #sim$duplicate_internal_cells <- T
  
  clocks = list()
  print(K_clocks)
  print(taus)
  
  i = 0
  print(i)
  print(sim)
  sim$add_mutant(name = paste0("Clone ", i), growth_rates = .1, death_rates = .01)
  sim$place_cell(paste0("Clone ", i), 500, 500)
  
  for (i in 1:K_clocks){
    print(i)
    print(sim)
    sim$add_mutant(name = paste0("Clone ", i), growth_rates = .1, death_rates = .01)
    
    sim$mutate_progeny(sim$choose_cell_in(paste0("Clone ", i-1)), paste0("Clone ", i))
    sim$update_rates(paste0("Clone ", i - 1 ), rates = c(growth = 0, death = .9))
    clocks<-append ( clocks, sim$get_clock()) #add
    print(paste0("clone ", i, " start: ", clocks[i]))
    print(paste0("clone ", i, " start: ",  sim$get_clock()))

    sim$run_up_to_time(taus[i])
  }
  
  sim$run_up_to_time(simulation_time)
  
  
  sim$add_mutant(name = paste0("Final clone"), growth_rates = 0.1, death_rates = .01)
  clocks<-append (clocks, sim$get_clock()) #add
  print(paste0("clone", K_clocks + 1, "start: ", clocks[K_clocks + 1]))
  
  sim$mutate_progeny(sim$choose_cell_in(paste0("Clone ", K_clocks)), paste0("Final clone"))
  
  sim$update_rates(paste0("Clone ", K_clocks ), rates = c(growth = 0, death = .9))
  
  sim$run_up_to_time(simulation_time + 200)
  
  print(paste0("End simulation :", sim$get_clock()))
  
  
  #plot_tissue(sim)
  #ggsave("tissue/tissue_00.pdf", dpi=300, width = 8, height = 8)
  #plot_muller(sim)
  #ggsave("tissue/muller_00.pdf", dpi=300, width = 8, height = 8)
  
  
  
  
  print(sim)
  
  
  ####### SAMPLE A #######
  n_w <- n_h <- 20
  ncells <- 0.99*n_w*n_h
  name = paste0("Clone ", K_clocks + 1)
  bbox <- sim$search_sample(c("Final clone" = ncells), n_w, n_h)
  sim$sample_cells("Sample A", bbox$lower_corner, bbox$upper_corner)
  
  t1 <- plot_tissue(sim, num_of_bins = 300)
  time_sample_A <- sim$get_clock() #add
  print(paste0("time_sample_A: ", time_sample_A)) #add
  
  
  
  
  
  
  
  # Save
  muller <- plot_muller(sim)
  forest <- sim$get_samples_forest()
  # forest$save("data/samples_forest.sff")
  
  plt_forest <- plot_forest(forest) %>%
    annotate_forest(forest)
  
  piechart <- plot_state(sim)
  timeseries <- plot_timeseries(sim)
  
  pl <- t1 + piechart + timeseries + muller + plt_forest + plot_layout(design = 'AABB\nEFGG\nHHHH\nHHHH\nHHHH')
  ggsave('plots/SPN04_tissue.png', plot = pl, width = 210, height = 297, units = "mm", dpi = 300)
  ggsave('plots/SPN04_tissue.pdf', plot = pl, width = 210, height = 297, units = "mm", dpi = 300)
  
  
  results <- list( samples_forest = forest, 
                   pl = pl)
  
  return(results)

}