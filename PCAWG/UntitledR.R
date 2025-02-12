# Load required packages
library(ggplot2)
library(dplyr)

# Compute Z-scores for CNA events and Clusters assigned
data <- df %>%
  mutate(
    Z_CNA = (n_events - mean(n_events)) / sd(n_events),
    Z_Clusters = (n_clusters - mean(n_clusters)) / sd(n_clusters)
  )

# Define hopeful monsters: High CNA Z-score (≥ 2) & Low Cluster Z-score (≤ -2)
z_threshold_high <- 1  # High threshold for CNA events
z_threshold_low <- -1  # Low threshold for clusters assigned

data <- data %>%
  mutate(Hopeful_Monster = (Z_CNA >= z_threshold_high) & (Z_Clusters <= z_threshold_low))

# Scatter plot to visualize hopeful monsters
ggplot(data, aes(x = n_events, y = n_clusters, color = Hopeful_Monster)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Hopeful Monster Detection (Z-Score Based)",
       x = "Number of CNA Events",
       y = "Number of Clusters Assigned",
       color = "Hopeful Monster")

# Display hopeful monster samples
data %>% filter(Hopeful_Monster == TRUE)
