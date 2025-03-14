# Load required libraries
library(ggplot2)
library(viridis)
library(dplyr)

setwd("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states")

dir.exists(file.path("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states"))
#Both <-read.table(file = "pomca_chr24_g10477.t1.simap_vs_omega.tsv", sep = '\t', header = TRUE)
Loss <-read.table(file = "./simapvsomega_allgenes/chr11_g33711.t1.simap_vs_omega.tsv",sep = '\t',  header = TRUE)
Gain <-read.table(file = "./simapvsomega_allgenes/chr02_g52043.t1.simap_vs_omega.tsv",sep = '\t', header = TRUE)



# Filter rows with highest weight within each group
#filtered_data <- Gain %>%
#  group_by(node_name) %>%
#  filter(weights == max(weights)) %>%
#  ungroup()
#print(filtered_data)

# Calculate the number of data points in each category
#category_counts <- filtered_data %>%
#  count(transition)
#print(category_counts)

# Calculate the average contrast for each category
#avg_contrasts <- filtered_data %>%
#  group_by(transition) %>%
#  summarise(avg_contrast = mean(contrast))
#print(avg_contrasts)

# Normalize the average contrast based on the number of data points in each category
#normalized_avg_contrasts <- avg_contrasts %>%
#  left_join(category_counts, by = "transition") %>%
#  mutate(normalized_avg_contrast = avg_contrast / n)
#print(normalized_avg_contrasts)

# Now you can use normalized_avg_contrasts for further analysis or visualization
#plot using a scatterplot of the weights 10.21.24

# Filter for 'gain' and 'loss' categories
#filtered_data_gain_loss <- filtered_data %>%
#  filter(transition %in% c("gain", "loss"))

# Scatterplot with labels for 'gain' and 'loss'
#ggplot(filtered_data, aes(x = weights, y = contrast, color = transition)) +
#  geom_point(size = 2) +  # Points for all categories
#  geom_text(data = filtered_data_gain_loss, aes(label = node_name), vjust = -1, size = 3) + # Labels for 'gain' and 'loss' categories
# labs(title = "gch2",
#       x = "Weights",
#       y = expression("Δω")) +
#  theme_classic() + 
#  theme(legend.position = "right") +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  theme(plot.title = element_text(face = "bold", size = 15)) +
#  theme(legend.key.size = unit(1, 'cm'))
#######################################################
# Step 1: Calculate Weight Difference
weight_summary <- Loss %>%
  mutate(weights = case_when(
    transition == "gain" ~ weights,
    transition == "same" ~ 0,
    transition == "loss" ~ -weights
  )) %>%
  group_by(node_name) %>%
  summarise(
    weight_difference = sum(weights, na.rm = TRUE),  # Total weight difference
    mean_contrast = mean(contrast, na.rm = TRUE),     # Mean contrast
    .groups = 'drop'
  )

# Display the summary data to check the results
print(weight_summary)

# Step 2: Specify the nodes for labeling
label_nodes <- c("Ane19_PBI", "GA023_LAT", "node_12", "GA058_TRI", 
                 "node_24", "node_8", "node_13", "GA056_MCC", 
                 "004_NIG", "GA057_EPH")

# Step 3: Create Scatter Plot with Gradient Color and Labels
ggplot(weight_summary, aes(x = weight_difference, y = mean_contrast, color = weight_difference)) +
  geom_point(size = 3) +  # Points for each node
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Line at y = 0
  scale_color_gradient(low = "blue", high = "red") +  # Gradient color from blue (close to -1) to red (close to 1)
  scale_x_continuous(limits = c(-1.1, 1.1), breaks = c(-1, 0, 1)) +  # Set x-axis limits and breaks
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +  # Adjust y-axis to allow space for labels
  theme_classic() + 
  labs(title = "ADAM17",
       x = "Weights",
       y = expression("Δω")) +
  theme(legend.position = "none") +  # Hide the legend
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold.italic", size = 15)) +
  geom_text(data = weight_summary %>% filter(node_name %in% label_nodes), 
            aes(label = node_name), 
            vjust = -1.5,  # Adjust vertical position of text
            size = 3)   # Adjust size of text

