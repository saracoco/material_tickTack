library(ggplot2)
library(dplyr)
library(tidyr)
#justify why I am using the median of tau MAP 



#' plotting_fit Function
#'
#' For non validation setting, this function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param res results
#' @param input_data input_data
#' @param all_sim all_sim
#' @param K K
#' @export
#' @examples
#' plotting_fit()

plotting_fit <- function(res, accepted_mutations, all_sim, K){
 
  S <- length(unique(accepted_mutations$segment_id))


  draws <- res$draws(format = "matrix")

  # check for NA

  names <- paste("tau[", 1:K, "]", sep = "")

  areas_tau <- mcmc_areas(
  draws,
  pars = names,
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 99%
  point_est = "median"
  )+
  labs(
  title = "Approximate Posterior distributions",
  subtitle = "With median and 80% and 95% intervals"
  )+
  xlim(0, 1) # + scale_x_continuous(breaks = c(1:5), labels = c("A", "B", "C", "D", "E"))

  color_scheme_set("blue")
  intervals_weigths_per_tau <- list()
  for (k in 1:K){
    names_weights <- paste("w[",1:S,",", k, "]", sep = "") 
    intervals_weigths_per_tau[[k]] <- mcmc_intervals(draws, pars = names_weights, point_est = "median", prob = 0.8, prob_outer = 0.95)+
      labs(
      title =  str_wrap( paste0("Posterior distributions of the weigths for tau ",k), width = 30 + K + sqrt(S)),
      subtitle = "With median and 80% and 95% intervals"
      )
  }
  intervals_weigths_per_tau <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol=K) #add global title

  #not using it, just need to add the general title to avoid repetitions
  intervals <- mcmc_intervals(draws, regex_pars = c("w"), point_est = "median", prob = 0.8, prob_outer = 0.95)+
    labs(
    title = "Posterior distributions of the weigths for each tau inferred",
    subtitle = "With median and 80% and 95% intervals"
    )
    
  Subtitle <- vector("list", (length(unique(accepted_mutations$segment_id))+1))
  
  
  Subtitle[[1]]  <- paste0("Number of mutations per segment: ")
  num_mutations_all <- c()
  for (i in seq_along(unique(accepted_mutations$segment_id))) {
    segment <- unique(accepted_mutations$segment_id)[i]
    num_mutations_single <- nrow(accepted_mutations %>% filter(segment_id == segment))
    num_mutations_all <- c(num_mutations_all, num_mutations_single)
    Subtitle[[i+1]] <- paste0(segment, "=", num_mutations_single," ")
  }
  Subtitle <- paste(Subtitle, collapse = "   ")
  cat(Subtitle)

  mean_mut <- mean(num_mutations_all)
  max_mut <- max(num_mutations_all)
  min_mut <- min(num_mutations_all)

  Subtitle_short <- paste0("Mutations per segment: average =", mean_mut, ",  min = ", min_mut, ", max = ", max_mut )
  
  print(paste0("accepted mutation plot filtered data", dim(accepted_mutations)))


  plot_filtered_data <- accepted_mutations %>%
  ggplot(mapping = aes(x = NV / DP, fill = segment_id)) +
  geom_histogram(alpha = .5, position = "identity", bins=60) +
  labs(x = "VAF")+
  labs(
  title = paste0( str_wrap("Histogram of the VAF spectrum, per segment, resulting from the simulation (only the data used in the inference after the filtering step are plotted here)", width = 90 + K + (S) ) ),
  subtitle = paste0( str_wrap(Subtitle_short, width = 90 + K + (S) ) )
  )+
  facet_wrap(vars(karyotype))



  final_plot <- plot_filtered_data /  (areas_tau) / intervals_weigths_per_tau +
    plot_layout(widths = c( 6, 6, 8), heights = c(15 + S + (S/1.3) , 15 + S + (S/1.3) + K, 15 + S + (S/1.3) + K*2)) +
    plot_annotation(
    title = paste0("Number of events: ", S),
    subtitle = " ",
    caption = ""
    ) & theme(text = element_text(size = 12+sqrt(S)), plot.title = element_text(size = 15+sqrt(S)), plot.subtitle = element_text(size = 12+sqrt(S)), axis.text = element_text(size = 12 + sqrt(S)), plot.caption = element_text(size = 10 + sqrt(S)))
    

  return(final_plot)
}







#' plotting_fit_all_k Function
#'
#' For non validation setting, this function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param res results
#' @param input_data input_data
#' @param all_sim all_sim
#' @param K K
#' @export
#' @examples
#' plotting_fit()

plotting_fit_all_k <- function(results, accepted_mutations, data){

  karyo <- data %>%
  group_by(segment_id) %>%
  summarise(karyotype = first(karyotype)) %>%
  pull(karyotype)

  if (length(karyo) <= 2){
    k_max = (length(karyo)) 
  } else if (length(karyo) <= 7){
    k_max = (length(karyo)-1) 
  } else if (length(karyo) <= 15) { 
    k_max = ((floor(length(karyo)/2))-1)
  } else{
    k_max = ceiling(sqrt(length(karyo))) 
  }

  model_selection_tibble <- dplyr::tibble()
  S <- length(unique(accepted_mutations$segment_id))

  for (K in 1:k_max) {
    
    res <- results[[K]]
    
    # Plot ELBO values over iterations
    output_files <- res$latent_dynamics_files()
    print(paste0("output_files ", output_files,"\n"))
    elbo_data <- read.csv(output_files, header = FALSE, comment.char = "#")
    colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
    iterations <- elbo_data$iter  # iteration column
    elbo_values <- elbo_data$ELBO  # ELBO column
    elbo_df <- data.frame(iteration = iterations, elbo = elbo_values)
        
    p <- ggplot(elbo_df, aes(x = iteration, y = elbo)) +
    geom_line(color = "blue") +
    labs(title = "ELBO Values over Iterations",
        x = "Iteration",
        y = "ELBO")
    saveRDS(p, paste0("./elbo_vs_iterations_",K,".rds"))  

    p <- plotting_fit(res, accepted_mutations, data, K)
        ggsave(paste0("./plots/plot_inference_",K,".png"), width = (12 + (S/2)), height = (16 + (S/2)), limitsize = FALSE, device = png, plot=p)




}
  # migliorare
    p_elbo_iter <- plotting_elbo(k_max)
    ggsave(paste0("./elbo_vs_iterations_.png"),width = (20 + (S/2)), height = (10 + (S/2)), plot = p_elbo_iter)

}




#' plotting_cluster_partition Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param res result
#' @param K number of mixture components
#' @param VALIDATION boolean
#' @export
#' @examples
#' plotting_cluster_partition()
 
plotting_cluster_partition <- function(res, K, input, VALIDATION = FALSE){

  accepted_mutations = input$accepted_mutations
    

    # Prepare data as for the MAE computation 
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    tau_inferred <- res$draws(names_tau, format = "matrix")
    tau_inferred_map <- lapply(1:ncol(tau_inferred), function(i) {median(tau_inferred[,i])} ) %>% unlist() 

    identity_matrix_RI = c()

    # Loop through each segment and determine its tau assignment
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
        segment <- unique(accepted_mutations$segment_id)[i]
        seg_modified <- str_split(segment, " ")
        segment_number <- seg_modified[[1]][2]

        names_weights <- paste("w[",segment_number,",", 1:K, "]", sep = "")  
        weights_inferred <- res$draws(names_weights, format = "matrix")
        weights_inferred_map <- lapply(1:ncol(weights_inferred), function(i) {median(weights_inferred[,i])} ) %>% unlist()

        tau_index_assigned <- which.max(weights_inferred_map)
        print(paste0("tau_index_assigned: ",tau_index_assigned))

        identity_matrix_RI <- c(identity_matrix_RI, tau_index_assigned)
    }
    print(paste0("identity_matrix_RI: ",identity_matrix_RI))
    print(paste0("tau_inferred_map: ",tau_inferred_map))

    # Prepare the data for ggplot
    segments <- unique(accepted_mutations$segment_id)
    plot_data <- data.frame(
        segment = segments,
        tau_index_assigned = identity_matrix_RI,
        tau_inferred_map = tau_inferred_map[identity_matrix_RI]
    )

    # Prepare identity matrix data
    identity_matrix_df <- data.frame(
        segment = rep(segments, each = K),
        tau_index = rep(1:K, times = length(segments)),
        value = as.numeric(rep(identity_matrix_RI, each = K) == rep(1:K, times = length(segments))),
        tau_inferred_map = rep(tau_inferred_map, times = length(segments))  # Correct the length mismatch
    )

    # Plot identity matrix with RColorBrewer palette
    library(RColorBrewer)

    # Plot identity matrix with RColorBrewer palette
    p <- ggplot(identity_matrix_df, aes(x = tau_index, y = segment)) +  # Change x to tau_index instead of tau_inferred_map
        geom_tile(aes(fill = factor(value)), color = "white") +
        scale_fill_brewer(palette = "Set1") + # Use a predefined palette
        labs(x = "Tau Index", y = "Segment", fill = "Membership") +
        theme_minimal(base_size = 14) +
        theme(
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            legend.background = element_rect(fill = "white"),
            legend.box.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 12, angle = 0, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.ticks = element_line(color = "grey70"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(color = "black"),
            plot.subtitle = element_text(color = "black")
        ) +
        scale_x_continuous(breaks = 1:K)  # Set x-axis to tau indices

    return(p) # Return the ggplot object
}





#' plotting_model_selection Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param results result
#' @export
#' @examples
#' plotting_model_selection()

# plotting_model_selection <- function(model_selection_tibble, best_K, entropy_per_segment_matrix_norm, entropy_per_segment_matrix, accepted_mutations){
plotting_model_selection <- function(model_selection_tibble, best_K, accepted_mutations, entropy){

  k_max = nrow(model_selection_tibble)
  model_selection = model_selection_tibble
# BIC plot
bic_plot <- ggplot(data = model_selection, aes(x = K, y = BIC)) + 
    geom_line(aes(colour = "BIC Line"), size = 1) + 
    geom_point(aes(colour = "BIC Point"), size = 3) +
    geom_point(aes(x = best_K, y = BIC[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[BIC == min(BIC)], y = min(BIC), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("BIC Line" = "steelblue", 
                                   "BIC Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  # Tema minimalista
    theme(legend.position = "top",  # Sposta la legenda sopra il grafico
          legend.title = element_blank(),  # Rimuove il titolo della legenda
          panel.grid.minor = element_blank(),  # Rimuove le griglie minori
          panel.grid.major.x = element_blank(),  # Rimuove le griglie verticali
          axis.title.x = element_text(margin = margin(t = 10)),  # Spazio per l'asse X
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))  # Spazio per l'asse Y


# AIC Plot
aic_plot <- ggplot(data = model_selection, aes(x = K, y = AIC)) + 
    geom_line(aes(colour = "AIC Line"), size = 1) + 
    geom_point(aes(colour = "AIC Point"), size = 3) +
    geom_point(aes(x = best_K, y = AIC[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[AIC == min(AIC)], y = min(AIC), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("AIC Line" = "steelblue", 
                                   "AIC Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  
    theme(legend.position = "top",  
          legend.title = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.grid.major.x = element_blank(),  
          axis.title.x = element_text(margin = margin(t = 10)),  
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))

# LOO Plot
loo_plot <- ggplot(data = model_selection, aes(x = K, y = LOO)) + 
    geom_line(aes(colour = "LOO Line"), size = 1) + 
    geom_point(aes(colour = "LOO Point"), size = 3) +
    geom_point(aes(x = best_K, y = LOO[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[LOO == min(LOO)], y = min(LOO), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("LOO Line" = "steelblue", 
                                   "LOO Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  
    theme(legend.position = "top",  
          legend.title = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.grid.major.x = element_blank(),  
          axis.title.x = element_text(margin = margin(t = 10)),  
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))

# Log Likelihood Plot
log_lik_plot <- ggplot(data = model_selection, aes(x = K, y = Log_lik)) + 
    geom_line(aes(colour = "Log likelihood Line"), size = 1) + 
    geom_point(aes(colour = "Log likelihood Point"), size = 3) +
    geom_point(aes(x = best_K, y = Log_lik[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[Log_lik == max(Log_lik)], y = max(Log_lik), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("Log likelihood Line" = "steelblue", 
                                   "Log likelihood Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  
    theme(legend.position = "top",  
          legend.title = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.grid.major.x = element_blank(),  
          axis.title.x = element_text(margin = margin(t = 10)),  
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))



# ICL plot
ICL_plot <- ggplot(data = model_selection, aes(x = K, y = ICL)) + 
    geom_line(aes(colour = "ICL Line"), size = 1) + 
    geom_point(aes(colour = "ICL Point"), size = 3) +
    geom_point(aes(x = best_K, y = ICL[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[ICL == min(ICL)], y = min(ICL), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("ICL Line" = "steelblue", 
                                   "ICL Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  # Tema minimalista
    theme(legend.position = "top",  # Sposta la legenda sopra il grafico
          legend.title = element_blank(),  # Rimuove il titolo della legenda
          panel.grid.minor = element_blank(),  # Rimuove le griglie minori
          panel.grid.major.x = element_blank(),  # Rimuove le griglie verticali
          axis.title.x = element_text(margin = margin(t = 10)),  # Spazio per l'asse X
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))  # Spazio per l'asse Y


# Update the BIC plot
bic_plot <- bic_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(BIC, 2)), nudge_y = -10, size = 3, fill = "white", color = "black") 
# Update the AIC plot
aic_plot <- aic_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(AIC, 2)), nudge_y = -2, size = 3, fill = "white", color = "black")
# Update the LOO plot
loo_plot <- loo_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(LOO, 2)), nudge_y = 50, size = 3, fill = "white", color = "black") 
# Update the Log Likelihood plot
log_lik_plot <- log_lik_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(Log_lik, 2)), nudge_y = -0.01, size = 3, fill = "white", color = "black") 
# Update the ICL plot
ICL_plot <- ICL_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(ICL, 2)), nudge_y = -10, size = 3, fill = "white", color = "black") 


  
  S <- length(unique(results$accepted_mutations$segment_id))
  Subtitle <- rep(NA, times=k_max)
  
  for (k in 1:k_max) {
    n_parameters <- K+((K-1)*S)+2 
    Subtitle[[k]] <- paste0(" ", k, ": ", n_parameters," parameters")
  }
  
  Subtitle <- paste(Subtitle, collapse = "\n")
  
  

# Combinazione dei grafici con layout e annotazioni
model_selection_plot <- (bic_plot | aic_plot) / (loo_plot | log_lik_plot) /(ICL_plot) +
  plot_annotation(
    title = "Model Selection Graphs: Score vs Number of Clusters",
    subtitle = paste0("Correspondence between number of clusters and number of parameters: \n", Subtitle),
    caption = "Source: Your Data",
    theme = theme_minimal(base_size = 12) +  # Tema minimalista con base 12
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
        plot.caption = element_text(size = 8, hjust = 0.5, margin = margin(t = 10)),
        plot.background = element_rect(fill = "white", color = NA),  # Sfondo bianco senza bordo
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray80", size = 0.2),  # Griglie sottili e discrete
        panel.grid.minor = element_blank()
      )
  ) & 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 10, face = "italic", margin = margin(b = 10)),
    axis.text = element_text(size = 8),
    plot.caption = element_text(size = 8)
  )

entropy_per_segment_matrix <- entropy$entropy_per_segment_matrix
entropy_per_segment_matrix_norm <- entropy$entropy_per_segment_matrix_norm

# plot the entropy behaviour 
  if (k_max>2){
  print(paste0("k_max: ", k_max))
  plot_entropy <- plotting_entropy(entropy_per_segment_matrix, entropy_per_segment_matrix_norm, k_max)
  ggsave(paste0("./plot_entropy.png"),width = (6 + (k_max*2)), height = (5), plot = plot_entropy)
  }


  return(model_selection_plot)
}





#' plotting_elbo Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param K number of clusters
#' @export
#' @examples
#' plotting_elbo()

plotting_elbo <- function(K){

  elbo_vs_iter <- list()
  for (k in 1:K){
        elbo_vs_iter[[k]] <- readRDS(paste0("./elbo_vs_iterations_",k,".rds"))+
                    labs(
                      title =  str_wrap( paste0("K = ",k))  # , width = 30 + K*2
                    )
  }
  elbo_vs_iter_plot <- gridExtra::grid.arrange(grobs = elbo_vs_iter, ncol=K) #add global title
  return(elbo_vs_iter_plot)

}



#' plotting_entropy Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param entropy_matrix number of clusters
#' @param entropy_matrix_norm number of clusters
#' @param K number of clusters
#' @export
#' @examples
#' plotting_entropy()

plotting_entropy <- function(entropy_per_segment_matrix, entropy_per_segment_matrix_norm, K){

  print("entropy_per_segment_matrix: ")
  print(entropy_per_segment_matrix)
  print("entropy_per_segment_matrix_norm: ")
  print(entropy_per_segment_matrix_norm)

  entropy_plot_list <- list()

  entropy_per_segment_matrix_trimmed <- entropy_per_segment_matrix[-1, ]
    df <- data.frame(
      k = 2:nrow(entropy_per_segment_matrix),  # Valori di k che vanno da 2 a k_max
      mean_per_column = rowMeans(entropy_per_segment_matrix_trimmed),  # Media delle colonne
      sum_per_column = rowSums(entropy_per_segment_matrix_trimmed)  # Somma delle colonne
    )

    entropy_plot_list[[1]] <- ggplot(df, aes(x = k, y = mean_per_column)) +
      geom_line() +
      geom_point() +
      labs(x = "k", y = "Mean across segment", title = paste0("Mean of entropy across segment for k from 2 to k_max = ", K), 
      subtitle = "The entropy of the distribution of w is multiplied by the number of mutations per segment") +
      theme_bw() #theme_minimal()




    entropy_per_segment_matrix_norm_trimmed <- entropy_per_segment_matrix_norm[-1, ]
    df <- data.frame(
      k = 2:nrow(entropy_per_segment_matrix_norm),  # Valori di k che vanno da 2 a k_max
      mean_per_column_norm = rowMeans(entropy_per_segment_matrix_norm_trimmed),  # Media delle colonne
      sum_per_column_norm = rowSums(entropy_per_segment_matrix_norm_trimmed)  # Somma delle colonne
    )

    entropy_plot_list[[2]] <- ggplot(df, aes(x = k, y = mean_per_column_norm)) +
      geom_line() +
      geom_point() +
      labs(x = "k", y = "Mean across segment", title = paste0("Mean of normalized entropy across segment for k from 2 to k_max = ", K)) +
      theme_bw() #theme_minimal()


  entropy_plot <- gridExtra::grid.arrange(grobs = entropy_plot_list, ncol=2) #add global title
  return(entropy_plot)

}




plot_simulate <- function(all_sim){

  segment_counts <- all_sim %>%
    group_by(tau) %>%
    summarise(num_segments = n_distinct(segment_id)) %>%
    ungroup()
  
  # plot for validation tau parameters of the inference ##############################################################
  tau_segments_plot <- ggplot(segment_counts, aes(x = factor(round(tau, 2)), y = num_segments, fill = factor(tau))) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(x = "Tau", y = "Number of Segments") +
    labs(
      title = paste0( str_wrap("Number of Unique Segments Associated with Each Tau Value", width = 40 ) ),
      subtitle = " "
    ) +
    theme_minimal()
      
  segment_counts <- all_sim %>%
    group_by(karyotype) %>%
    summarise(num_segments = n_distinct(segment_id)) %>%
    ungroup()
  
  karyo_segments_plot <- ggplot(segment_counts, aes(x = factor(karyotype), y = num_segments, fill = factor(karyotype))) +
    scale_fill_manual(values = c("skyblue", "orange", "lightgreen")) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(x = "Karyotype", y = "Number of Segments") +
    labs(
      title = paste0( str_wrap("Number of Unique Segments Associated with Each Karyotype configuration", width = 40 ) ),
      subtitle = " "
    ) +
    theme_minimal()


  simulation_plot <- (tau_segments_plot | karyo_segments_plot) +
      plot_layout(widths = 25, heights = 10 ) +
      plot_annotation(
        title = " ",
        subtitle = " ",
        caption = ""
      )

  return(simulation_plot)


}