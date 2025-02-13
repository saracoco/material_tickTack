# Load required packages
library(dbscan)
library(ggplot2)
library(dplyr)

# Compute Z-scores for CNA events and Clusters assigned
data <- df %>%
  dplyr::rename(CNA_events=n_events, Clusters_assigned=n_clusters)

# Prepare data for DBSCAN (excluding Sample_ID)
dbscan_data <- data[, c("CNA_events", "Clusters_assigned")]

# Perform DBSCAN clustering
dbscan_result <- dbscan(dbscan_data, eps = 5, minPts = 20)  # Adjust eps & minPts if needed

# Add DBSCAN cluster labels to the dataset
data$DBSCAN_Cluster <- dbscan_result$cluster

# Identify hopeful monsters: DBSCAN labels outliers as "0"
data <- data %>%
  mutate(Hopeful_Monster = DBSCAN_Cluster == 0)

# Scatter plot to visualize DBSCAN clusters and hopeful monsters
ggplot(data, aes(x = CNA_events, y = Clusters_assigned, color = as.factor(DBSCAN_Cluster))) +
  geom_point(size = 3) +
  theme_minimal() +
  facet_wrap(~ploidy) +
  labs(title = "Hopeful Monster Detection Using DBSCAN",
       x = "Number of CNA Events",
       y = "Number of Clusters Assigned",
       color = "DBSCAN Cluster")

# Display hopeful monster samples
data %>% filter(Hopeful_Monster == TRUE)
