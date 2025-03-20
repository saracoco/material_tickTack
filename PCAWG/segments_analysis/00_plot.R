.libPaths(new="~/R/rstudio_v3/") 

rm(list= ls())
library(ggplot2)
library(dplyr)
smoothing = "5+e05"

info_segments <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/segments_analysis/data/info_segments_",smoothing,".rds"))
info_segments

fittable_samples <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/data/fittable_samples.RDS")) #fittable_samples_5Mb_smoothed_15mm_0.4pi
info_segments <- info_segments %>% filter(sample %in% basename(fittable_samples))
#plot 1
# Count the number of rows per type
type_counts <- info_segments %>%
  group_by(type) %>%
  summarise(n = n()) %>%
  mutate(facet_label = paste0(type, " (n=", n, ")"))

# Merge with original data to use updated facet labels
info_segments <- info_segments %>%
  left_join(type_counts, by = "type")

# Calculate the median of median_length per type
median_per_type <- info_segments %>%
  group_by(facet_label) %>%
  summarise(median_length = median(median_length, na.rm = TRUE), .groups = 'drop') 


# Create the plot
ggplot(info_segments, aes(x = n_cna_simple, y = median_length)) +
  geom_point(size = 0.4, alpha = 0.7, color = "blue") +  # Scatter points
  facet_wrap(~ facet_label) +  # Facet by type with row counts
  geom_hline(data = median_per_type, aes(yintercept = median_length), 
             color = "red", linetype = "dashed", linewidth = 0.6) +  # Add median line
  geom_text(data = median_per_type, aes(x = Inf, y = median_length, 
                                        label = paste("Median:", round((median_length/10**6),2))),
            vjust = -1, hjust = 1, size = 3, color = "black") +  # Add median value text
  labs(
    title = "Scatter Plot of n_cna_simple vs. median_length per Type",
    x = "Number of CNA Simple",
    y = "Median Length",
    caption = paste0("Max_distance applied in the CNAqc smoothing function:",smoothing)
  ) +
  theme_minimal() +  # Clean theme
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Style facet labels
    panel.background = element_rect(fill = "white"),  # White background for the panels
    plot.background = element_rect(fill = "white")  # White background for the whole plot
  )

ggsave(paste0("./plot/scatter_plot_with_counts_",smoothing,".pdf"), width = 10, height = 8, dpi = 300)






# plot 2
l <- 10
info_segments_binned <- info_segments %>%
  mutate(bin = cut_number(median_length, n = l)) %>%
  group_by(bin) %>%
  summarise(mean_mutation_number = mean(meadian_mutation_number, na.rm = TRUE), .groups = "drop")

ggplot(info_segments_binned, aes(x = bin, y = mean_mutation_number)) +
  geom_col(fill = "steelblue", color = "black") +
  labs(
    x = "Median Length (Binned into l Groups)",
    y = "Mean of Median Mutation Number",
    title = paste("Histogram of Median Mutation Number (", l, "Bins)", sep = ""),
    caption = paste0("Max_distance applied in the CNAqc smoothing function:",smoothing)
  ) +
  geom_text(aes(label = round(mean_mutation_number,2)), vjust = -0.3, size = 3) +  # Add labels on top of bars
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),  # White background for the panels
        plot.background = element_rect(fill = "white") )  # Rotate x-axis labels for readability
ggsave(paste0("./plot/mutations_per_segment_all_",smoothing,".png"), width = 10, height = 8, dpi = 300)


# faceting by type
unique_types <- unique(info_segments$type)
split_idx <- ceiling(length(unique_types) / 2)  # Split into two roughly equal parts

types_group1 <- unique_types[1:split_idx]
types_group2 <- unique_types[(split_idx + 1):length(unique_types)]

# Function to generate the binned data
bin_data <- function(df, types) {
  df %>%
    filter(type %in% types) %>%  # Select only relevant types
    group_by(type) %>%
    mutate(bin = cut_number(median_length, n = l)) %>%
    group_by(type, bin) %>%
    summarise(mean_mutation_number = mean(meadian_mutation_number, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      bin = substr(as.character(bin), 1, 10),  # Truncate bin labels
      mean_mutation_number = round(mean_mutation_number, 2)  # Round values for better readability
    )
}

# Generate data for each group
info_segments_binned_1 <- bin_data(info_segments, types_group1)
types_group2 <- types_group2[(types_group2!= "DLBC") & (types_group2!= "KICH") &  (types_group2!="LIHC") & (types_group2!="READ")]
info_segments_binned_2 <- bin_data(info_segments, types_group2)

# Function to plot the data with value labels
plot_histogram <- function(data, title_suffix) {
  ggplot(data, aes(x = bin, y = mean_mutation_number)) +
    geom_col(fill = "steelblue", color = "black") +
    geom_text(aes(label = mean_mutation_number), vjust = -0.3, size = 3) +  # Add labels on top of bars
    facet_wrap(~type, scales = "free_x") +  # Facet by type
    labs(
      x = "Median Length (Binned into l Groups)",
      y = "Mean of Median Mutation Number per segment",
      title = paste("Histogram of Median Mutation Number -", title_suffix),
      caption = paste0("Max_distance applied in the CNAqc smoothing function:",smoothing)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      strip.text = element_text(size = 10, face = "bold"),
      panel.background = element_rect(fill = "white"),  # White background for the panels
      plot.background = element_rect(fill = "white") 
    )
}

# Create plots
plot1 <- plot_histogram(info_segments_binned_1, "Group 1")
plot2 <- plot_histogram(info_segments_binned_2, "Group 2")

# Display both plots
plot1
plot2

ggsave(paste0("./plot/mutations_per_segment_per_type_",smoothing,"_1.png"), width = 20, height = 12, dpi = 300, plot=plot1)
ggsave(paste0("./plot/mutations_per_segment_per_type_",smoothing,"_2.png"), width = 20, height = 12, dpi = 300, plot=plot2)




# karyotype distribution 
library(ggplot2)
library(dplyr)
library(forcats) 

info_karyo <- readRDS("./data/info_complex_segments_5e+06.rds")
df <- info_karyo

percentage_threshold <- 0.01 
df %>%
  filter(karyotype != "1:1") %>% 
  count(karyotype) %>%
  mutate(percent = n / sum(n)) %>%  
  filter(percent >= percentage_threshold) %>%  
  ggplot(aes(y = fct_reorder(karyotype, n), x = n, fill = karyotype)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 14) + 
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Karyotype Distribution (Filtered by Percentage (> 0.01), Excluding 1:1)",
       y = "Karyotype",
       x = "Count") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(title = "Karyotype"))
ggsave(paste0("./plot/Overall_distribution_of_karyotypes_",smoothing,".png"), width = 20, height = 12, dpi = 300)

# percent 
df %>%
  filter(karyotype != "1:1") %>% 
  count(karyotype) %>%
  mutate(percent = n / sum(n)) %>%  
  filter(percent >= percentage_threshold) %>%  
  ggplot(aes(y = fct_reorder(karyotype, percent), x = percent, fill = karyotype)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 14) + 
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Karyotype Distribution (Filtered by Percentage (> 0.01), Excluding 1:1)",
       y = "Karyotype",
       x = "Count") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(title = "Karyotype"))
ggsave(paste0("./plot/Overall_distribution_of_karyotypes_percent_",smoothing,".png"), width = 20, height = 12, dpi = 300)



# per tumour Type
# count
percentage_threshold <- 0.06
df %>%
  filter(karyotype != "1:1") %>%  
  count(type, karyotype) %>%
  group_by(type) %>%
  mutate(percent = n / sum(n)) %>% 
  filter(percent >= percentage_threshold) %>% 
  ungroup() %>%
  ggplot(aes(x = fct_reorder(karyotype, n), y = n, fill = karyotype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type, scales = "free_y") +
  coord_flip() +
  theme_classic(base_size = 14) +  # White background
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = paste0("Karyotype Distribution Across Tumour Types (Filtered, Excluding 1:1, percentage > ", percentage_threshold, " )"),
       x = "Karyotype",
       y = "Count") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(title = "Karyotype"))  
ggsave(paste0("./plot/Per_tumour_distribution_of_karyotypes_",smoothing,".png"), width = 20, height = 12, dpi = 300)


# percent
percentage_threshold <- 0.06
df %>%
  filter(karyotype != "1:1") %>%  
  count(type, karyotype) %>%
  group_by(type) %>%
  mutate(percent = n / sum(n)) %>% 
  filter(percent >= percentage_threshold) %>% 
  ungroup() %>%
  ggplot(aes(x = fct_reorder(karyotype, percent), y = percent, fill = karyotype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type, scales = "free_y") +
  coord_flip() +
  theme_classic(base_size = 14) +  # White background
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = paste0("Karyotype Distribution Across Tumour Types (Filtered, Excluding 1:1, percentage > ", percentage_threshold, " )"),
       x = "Karyotype",
       y = "Count") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(title = "Karyotype")) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  # Add the total number of samples per tumor type
  geom_text(data = df %>%
    filter(karyotype != "1:1") %>%
    group_by(type) %>%
    summarise(total_samples = n_distinct(sample)),  # Count distinct samples
    aes(x = 1, y = 0, label = paste("N:", total_samples)), 
    inherit.aes = FALSE, 
    size = 3, 
    hjust = -3, 
    vjust = 0, 
    color = "black")  # Add total sample count
ggsave(paste0("./plot/Per_tumour_distribution_of_karyotypes_percent_",smoothing,".png"), width = 20, height = 12, dpi = 300)









# distribution of 3 most common karyotypes per samples
percentage_threshold <- 0.05

# Step 1: Count occurrences of each karyotype per sample
top_karyotypes_per_sample <- df %>%
  filter(karyotype != "1:1") %>%  
  group_by(sample, karyotype, type) %>%
  summarise(count = n(), .groups = "drop")  # Count the occurrences of each karyotype per sample

# Step 2: Calculate frequency for each karyotype within each sample
top_karyotypes_per_sample <- top_karyotypes_per_sample %>%
  group_by(sample) %>%
  mutate(frequency = count / sum(count)) %>%  # Calculate frequency for each karyotype
  ungroup()

# Step 3: Get the top 3 most frequent karyotypes for each sample
top_3_karyotypes <- top_karyotypes_per_sample %>%
  group_by(sample) %>%
  slice_max(order_by = count, n = 3, with_ties = TRUE) %>%  # Select top 3 by count (or frequency)
  ungroup()



percentage_threshold <- 0.01 
top_3_karyotypes %>%
  filter(karyotype != "1:1") %>% 
  count(karyotype) %>%
  mutate(percent = n / sum(n)) %>%  
  filter(percent >= percentage_threshold) %>%  
  ggplot(aes(y = fct_reorder(karyotype, n), x = n, fill = karyotype)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 14) + 
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Karyotype Distribution (Filtered by Percentage (> 0.01), Excluding 1:1)",
       y = "Karyotype",
       x = "Count") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(title = "Karyotype"))
ggsave(paste0("./plot/3_Most_common_distribution_of_karyotypes_",smoothing,".png"), width = 20, height = 12, dpi = 300)



percentage_threshold <- 0.05
df %>%
  filter(karyotype != "1:1") %>%  
  count(type, karyotype) %>%
  group_by(type) %>%
  mutate(percent = n / sum(n)) %>% 
  filter(percent >= percentage_threshold) %>% 
  ungroup() %>%
  ggplot(aes(x = fct_reorder(karyotype, percent), y = percent, fill = karyotype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type, scales = "free_y") +
  coord_flip() +
  theme_classic(base_size = 14) +  # White background
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = paste0("Karyotype Distribution Across Tumour Types (Filtered, Excluding 1:1, percentage > ", percentage_threshold, " )"),
       x = "Karyotype",
       y = "Count") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(title = "Karyotype")) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  # Add the total number of samples per tumor type
  geom_text(data = df %>%
              filter(karyotype != "1:1") %>%
              group_by(type) %>%
              summarise(total_samples = n_distinct(sample)),  # Count distinct samples
            aes(x = 1, y = 0, label = paste("N:", total_samples)), 
            inherit.aes = FALSE, 
            size = 3, 
            hjust = -3, 
            vjust = 0, 
            color = "black")  # Add total sample count
ggsave(paste0("./plot/3_Most_common_distribution_Per_tumour_",smoothing,".png"), width = 20, height = 12, dpi = 300)





# 
# 
# outliers <- info_segments%>%filter(n_cna_simple>100)
# outliers_first_per_type <- info_segments %>%
#   filter(n_cna_simple > 100) %>%
#   group_by(facet_label) %>%
#   slice(1) %>%  # Select the first row per type
#   ungroup() 
# 
# 
# outliers_first_per_type <- outliers
# 
# type_counts <- outliers_first_per_type %>%
#   group_by(type) %>%
#   summarise(n = n()) %>%
#   mutate(facet_label = paste0(type, " (n=", n, ")"))
# 
# # Merge with original data to use updated facet labels
# outliers_first_per_type <- outliers_first_per_type %>%
#   left_join(type_counts, by = "type")
# 
# # Calculate the median of median_length per type
# median_per_type <- outliers_first_per_type %>%
#   group_by(facet_label) %>%
#   summarise(median_length = median(median_length, na.rm = TRUE), .groups = 'drop')
# 
# # Create the plot
# ggplot(outliers_first_per_type, aes(x = n_cna_simple, y = median_length)) +
#   geom_point(size = 1, alpha = 0.7, color = "blue") +  # Scatter points
#   facet_wrap(~ facet_label) +  # Facet by type with row counts
#   geom_hline(data = median_per_type, aes(yintercept = median_length), 
#              color = "red", linetype = "dashed", linewidth = 0.6) +  # Add median line
#   geom_text(data = median_per_type, aes(x = Inf, y = median_length, 
#                                         label = paste("Median:", round((median_length/10**6),2))),
#             vjust = -1, hjust = 1, size = 3, color = "black") +  # Add median value text
#   labs(
#     title = "Scatter Plot of n_cna_simple vs. median_length per Type",
#     x = "Number of CNA Simple",
#     y = "Median Length"
#   ) +
#   theme_minimal() +  # Clean theme
#   theme(
#     strip.text = element_text(size = 12, face = "bold"),  # Style facet labels
#     panel.background = element_rect(fill = "white"),  # White background for the panels
#     plot.background = element_rect(fill = "white")  # White background for the whole plot
#   )
# 



