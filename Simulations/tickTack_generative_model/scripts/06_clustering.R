
require(transport)
require(tidy)
library(cluster)

res_SingleTT$inference_results %>% 
  dplyr::filter(tau <= 1) %>% 
  ggplot(mapping = aes(x=tau, fill=as.factor(segment))) +
  geom_histogram(bins = 100)

library(ggplot2)
library(reshape2) # For melting the matrix

plot_heatmap <- function(matrix_data, 
                         x_label = "Columns", 
                         y_label = "Rows", 
                         fill_label = "Value", 
                         low_color = "blue", 
                         high_color = "red", 
                         title = "Heatmap") {
  # Check if input is a matrix
  if (!is.matrix(matrix_data)) {
    stop("Input must be a matrix.")
  }
  
  # Convert matrix to a long-format data frame
  df <- melt(matrix_data)
  colnames(df) <- c("Row", "Column", "Value")
  
  # Create the heatmap
  ggplot(df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = low_color, high = high_color, name = fill_label) +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}


compute_metrics <- function(posterior_draws) {
  # Input: 
  # posterior_draws: A list where each element contains posterior draws for a sample
  
  # Output:
  # A matrix with Wasserstein distance, KL divergence, and overlap-based similarity
  
  n <- length(posterior_draws)
  results <- matrix(NA, nrow = n, ncol = n, 
                    dimnames = list(paste0("Sample", 1:n), paste0("Sample", 1:n)))
  
  wasserstein <- results
  kl_divergence <- results
  overlap_similarity <- results
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Extract draws for sample i and j
      draws_i <- posterior_draws[[i]]
      draws_j <- posterior_draws[[j]]
      
      # Compute Wasserstein distance
      wasserstein[i, j] <- wasserstein1d(draws_i, draws_j)
      wasserstein[j, i] <- wasserstein[i, j]
      
      # Compute KL divergence using density estimates
      density_i <- density(draws_i)
      density_j <- density(draws_j)
      
      # Interpolate densities to match the support
      common_x <- sort(unique(c(density_i$x, density_j$x)))
      p <- approx(density_i$x, density_i$y, xout = common_x, rule = 2)$y
      q <- approx(density_j$x, density_j$y, xout = common_x, rule = 2)$y
      
      # Normalize to avoid NaN issues
      p <- p / sum(p)
      q <- q / sum(q)
      
      kl_divergence[i, j] <- sum(ifelse(p > 0 & q > 0, p * log(p / q), 0))
      kl_divergence[j, i] <- sum(ifelse(q > 0 & p > 0, q * log(q / p), 0))
      
      # Compute overlap-based similarity
      overlap_similarity[i, j] <- sum(pmin(p, q))
      overlap_similarity[j, i] <- overlap_similarity[i, j]
    }
  }
  
  list(
    Wasserstein = wasserstein,
    KL_Divergence = kl_divergence,
    Overlap_Similarity = overlap_similarity
  )
}


draws = lapply(res_SingleTT$inference_results$segment %>% unique(), function(idx) {
  res_SingleTT$inference_result %>% 
    dplyr::filter(segment == idx) %>% 
    dplyr::filter(tau <= 1, tau >= 0) %>% 
    pull(tau)
})

hierarchical_optimal_clusters <- function(metrics_matrix, max_clusters = 10, method = "ward.D2") {
  # Check if the input is a valid matrix
  if (!is.matrix(metrics_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Convert the metrics matrix to a distance object
  dist_matrix <- as.dist(metrics_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = method)
  
  # Initialize variables to store evaluation results
  silhouette_scores <- numeric(max_clusters - 1)
  
  # Loop over the number of clusters (from 2 to max_clusters)
  for (k in 2:max_clusters) {
    # Cut the dendrogram into k clusters
    cluster_assignments <- cutree(hc, k = k)
    
    # Calculate silhouette scores for the current k
    silhouette_result <- silhouette(cluster_assignments, dist_matrix)
    silhouette_scores[k - 1] <- mean(silhouette_result[, 3]) # Mean silhouette width
  }
  
  # Find the optimal number of clusters (max silhouette score)
  optimal_clusters <- which.max(silhouette_scores) + 1
  
  # Return results
  return(list(
    optimal_clusters = optimal_clusters,
    silhouette_scores = silhouette_scores,
    hc = hc # Return the hierarchical clustering object for plotting or further analysis
  ))
}

cluster_hierarchical <- function(metrics_matrix, num_clusters) {
  # Check if the input is a valid matrix
  if (!is.matrix(metrics_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Convert the metrics matrix to a distance object
  dist_matrix <- as.dist(metrics_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  # Cut the dendrogram into clusters
  clusters <- cutree(hc, k = num_clusters)
  
  # Return the cluster assignments
  return(clusters)
}

metrics = compute_metrics(draws)
plot_heatmap(metrics$Wasserstein)
plot_heatmap(metrics$Overlap_Similarity)
plot_heatmap(metrics$KL_Divergence)

hierarchical_optimal_clusters(metrics$Overlap_Similarity, max_clusters = 9)
clustering = cluster_hierarchical(metrics$Overlap_Similarity, 3)
true_clust = compare_assignment$real_clocks %>% 
  as.factor() %>% 
  as.numeric()

fossil::rand.index(clustering, true_clust)











































