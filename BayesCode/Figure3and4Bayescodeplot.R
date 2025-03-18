# Load required libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(fishualize)
#fish_palettes()

setwd("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states")

# Define function to load and process each gene's data
process_gene_data <- function(file_path, gene_name) {
  data <- read.table(file = file_path, sep = '\t', header = TRUE)
  data$gene_name <- gene_name  
  return(data)
}

# Load data for Loss and Gain genes
Loss_files <- c(
  "./simapvsomega_allgenes/chr15_g33521.t1.simap_vs_omega.tsv", # vps11
  "./simapvsomega_allgenes/chr04_g50543.t1.simap_vs_omega.tsv", # th
  "./simapvsomega_allgenes/chr09_g57882.t1.simap_vs_omega.tsv", # drd2a
  "./simapvsomega_allgenes/chr18_g29860.t1.simap_vs_omega.tsv", # gchfr
  "./simapvsomega_allgenes/chr16_g55920.t1.simap_vs_omega.tsv", # spr(a)
  "./simapvsomega_allgenes/chr11_g33711.t1.simap_vs_omega.tsv" # ADAM17
)

Gain_files <- c(
  "./simapvsomega_allgenes/chr19_g19445.t1.simap_vs_omega.tsv", # gch2
  "./simapvsomega_allgenes/chr01_g56646.t1.simap_vs_omega.tsv", # leo1
  "./simapvsomega_allgenes/chr04_g12514.t1.simap_vs_omega.tsv", # ANKRD27
  "./simapvsomega_allgenes/chr03_g37341.t1.simap_vs_omega.tsv", # slc2a15b
  "./simapvsomega_allgenes/chr02_g52043.t1.simap_vs_omega.tsv" # OCA2
)

# Gene names for Loss and Gain genes
Loss_gene_names <- c("vps11", "th", "drd2a", "gchfr", "spra", "adam17b")
Gain_gene_names <- c("gch2", "leo1", "ankrd27", "slc2a15b", "oca2")

# Load and process Loss and Gain gene data
Loss_data <- lapply(1:length(Loss_files), function(i) {
  process_gene_data(Loss_files[i], Loss_gene_names[i])  
})

Gain_data <- lapply(1:length(Gain_files), function(i) {
  process_gene_data(Gain_files[i], Gain_gene_names[i])  
})


# Combine the lists into data frames
Loss_data <- bind_rows(Loss_data)
Gain_data <- bind_rows(Gain_data)

# Step 1: Calculate Weight Difference for each gene
calculate_weight_summary <- function(data) {
  data %>%
    mutate(weights = case_when(
      transition == "gain" ~ weights,
      transition == "same" ~ 0,
      transition == "loss" ~ -weights
    )) %>%
    group_by(node_name, gene_name) %>%
    summarise(
      weight_difference = sum(weights, na.rm = TRUE),
      mean_contrast = mean(contrast, na.rm = TRUE),
      .groups = 'drop'
    )
}

# Apply the function to Loss and Gain data
weight_summary_loss <- calculate_weight_summary(Loss_data)
weight_summary_gain <- calculate_weight_summary(Gain_data)

# Display the summary data for loss and gain genes
print(weight_summary_loss)
print(weight_summary_gain)

# Specify the nodes for labeling (adjust as necessary)
label_nodes <- c("Ane19_PBI", "GA023_LAT", "node_12", "GA058_TRI", 
                 "node_24", "node_8", "node_13", "GA056_MCC", 
                 "004_NIG", "GA057_EPH")





#############################################################################
#########################   For the loss genes      #########################
#############################################################################

# Step 1: Get the node names that are below the threshold (-0.5) for the Loss plot
colored_nodes_loss <- weight_summary_loss %>%
  filter(weight_difference < -0.5) %>%
  pull(node_name) %>%
  unique()

# Step 2: Generate contrasting colors for these nodes 
node_colors_loss <- setNames(fish(n = length(colored_nodes_loss), option = "Clepticus_parrae"), colored_nodes_loss)
#print(node_colors_loss)

# Step 3: Add the "below_threshold" category with gray for non-colored nodes
# This will ensure the rest of the nodes are gray if they don't meet the condition
non_colored_nodes_loss <- setdiff(label_nodes, colored_nodes_loss)
node_colors_loss <- c(node_colors_loss, setNames("lightgrey", "below_threshold"))

plot_loss <- ggplot(weight_summary_loss, aes(x = weight_difference, y = mean_contrast)) +
  geom_point(aes(color = ifelse(weight_difference < -0.5, node_name, "below_threshold")), size = 3) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +  
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +  
  
  # Apply the generated node colors (ensure that 'below_threshold' nodes are gray)
  scale_color_manual(values = node_colors_loss, name = "Node Label") +  
  
  scale_x_continuous(limits = c(-1.1, 1.1), breaks = c(-1, 0, 1)) +  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +  
  theme_classic() +  
  labs(title = "Significant Loss Genes", x = "Weight Difference", y = expression("Δω")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15)) +  
  
  # Facet by gene
  facet_grid(gene_name ~ ., scales = "free_y", switch = "both") +  
  
  theme(
    strip.text.y = element_text(face = "bold.italic", size = 10, angle = 0),  
    axis.title.x = element_text(margin = margin(t = 10)),  
    axis.title.y = element_text(size = 12),  
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  
    axis.line.x = element_line(),  
    axis.ticks.x = element_line(),  
    panel.spacing = unit(1, "lines"),  
    strip.placement = "outside"
  )

print(plot_loss)

ggsave("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/plot_loss.png", plot = plot_loss, width = 8, height = 6, dpi = 300)

#############################################################################
#########################   For the gain genes      #########################
#############################################################################
# Step 1: Get the node names that are above the threshold (0.5) for the Gain plot
colored_nodes_gain <- weight_summary_gain %>% 
  filter(weight_difference > 0.5) %>%  # Change the condition to > 0.5
  pull(node_name) %>% 
  unique()

# Step 2: Generate contrasting colors for these nodes 
node_colors_gain <- setNames(fish(n = length(colored_nodes_gain), option = "Chlorurus_spilurus"), colored_nodes_gain)

# Step 3: Add the "below_threshold" category with gray for non-colored nodes
non_colored_nodes_gain <- setdiff(label_nodes, colored_nodes_gain)
# Add 'below_threshold' explicitly to the color mapping
node_colors_gain <- c(node_colors_gain, setNames("lightgrey", "below_threshold"))

# Create the scatter plot for Gain genes with conditional coloring
plot_gain <- ggplot(weight_summary_gain, aes(x = weight_difference, y = mean_contrast)) +  
  geom_point(aes(color = ifelse(weight_difference > 0.5, node_name, "below_threshold")), size = 3) +  # Use the > 0.5 condition for coloring
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") + 
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") + 
  scale_color_manual(values = node_colors_gain, name = "Node Label") +  # Apply the node colors for the Gain plot
  scale_x_continuous(limits = c(-1.1, 1.1), breaks = c(-1, 0, 1)) +  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +  
  theme_classic() +  
  labs(title = "Significant Gain Genes", x = "Weight Difference", y = expression("Δω")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15)) + 
  facet_grid(gene_name ~ ., scales = "free_y", switch = "both") +  
  theme(
    strip.text.y = element_text(face = "bold.italic", size = 10, angle = 0),  
    axis.title.x = element_text(margin = margin(t = 10)),   
    axis.title.y = element_text(size = 12),  
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  
    axis.line.x = element_line(),  
    axis.ticks.x = element_line(),  
    panel.spacing = unit(1, "lines"),  
    strip.placement = "outside"
  )

# Display the plot
print(plot_gain)

ggsave("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/plot_gain.png", plot = plot_gain, width = 8, height = 6, dpi = 300) 

######################################################################
###########            Add tree to the plot                ###########
######################################################################
# Load necessary libraries

library(ape)
#############################################################################
#########################   For the loss genes      #########################
#############################################################################

# Read tree
tree <- read.tree("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Data/calibrated_tree.tre") # /ASR/Data/calibrated_tree.tre

# Ensure `colored_node_hex_loss` is a named vector, and the names correspond to branch IDs
colored_node_hex_loss <- c("2" = "#7C43A9FF", "14" = "#E28928FF", "23" = "#C464AAFF", "40" = "#0646CFFF", "41" = "#001847FF")
loss_nodes <- c(2, 14, 23, 40, 41)  # For loss branches

# Check for missing nodes in your branch-color mapping 
missing_nodes <- setdiff(loss_nodes, names(colored_node_hex_loss))
if(length(missing_nodes) > 0) {
  print(paste("Missing colors for nodes:", paste(missing_nodes, collapse = ", ")))
}

# Initialize edge color vector (darkgrey for all edges by default)
edge_colors <- rep("lightgrey", length(tree$edge[, 1]))

# Loop through the edges and color those that correspond to the blue nodes (loss branches)
for (i in 1:length(tree$edge[, 1])) {
  # Check if the target node (tree$edge[i, 2]) is in the list of blue_nodes (i.e., loss branches)
  if (tree$edge[i, 2] %in% loss_nodes) {
    # Use the custom color from `colored_node_hex_loss` for the corresponding branch
    if (as.character(tree$edge[i, 2]) %in% names(colored_node_hex_loss)) {
      edge_colors[i] <- colored_node_hex_loss[as.character(tree$edge[i, 2])]
    } else {
      edge_colors[i] <- "lightgray"  # Default color if no match found
    }
  }
}

#plot(tree, type = "phylogram", edge.col = edge_colors, cex = 1.2, edge.width = 10, show.tip.label = TRUE)
#nodelabels(node = 1:(length(tree$tip.label) + length(tree$edge) - 1), cex = 0.5, frame = "none", adj = c(1.5, 0))

# Plot the tree with the custom colors for the loss branches
png("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/losstree_plot.png", width = 2000, height = 2000, res = 300)
par(mar = c(2, 2, 2, 2)) 
plot(tree, type = "phylogram", edge.col = edge_colors, cex = 1.2, edge.width = 10, show.tip.label = TRUE)
dev.off() 

# Optionally, add node labels to visualize the internal node numbers
#nodelabels(node = 1:(length(tree$tip.label) + length(tree$edge) - 1), cex = 0.5, frame = "none", adj = c(1.5, 0))

library(grid)
library(gridExtra)

# Read the saved images
loss_genes <- rasterGrob(png::readPNG("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/plot_loss.png"), interpolate = TRUE)
loss_tree <- rasterGrob(png::readPNG("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/losstree_plot.png"), interpolate = TRUE)

# Define layout matrix
layout <- matrix(c(1, 2,
                   1, NA),  # NA leaves an empty space
                 nrow = 2, byrow = TRUE)

# Arrange side-by-side
pdf("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/Loss_Genes.pdf", width = 14, height = 8) #Figure 4 in paper
grid.arrange(loss_genes, loss_tree,layout_matrix = layout, widths = c(3, 1), heights = c(1, 1))
dev.off() 

#############################################################################
#########################   For the gain genes      #########################
#############################################################################

# Read tree
tree <- read.tree("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Data/calibrated_tree.tre")

# Ensure `colored_node_hex_loss` is a named vector, and the names correspond to branch IDs
colored_node_hex_gain <- c("6" = "#00A5FFFF", "25" = "#4A79BCFF", "28" = "#A790DBFF", "39" = "#04EAB8FF", "52" = "#D8D643FF")
red_nodes <- c(6, 25, 28, 39, 52)  # For gain branches

print(node_colors_gain)
# Check for missing nodes in your branch-color mapping (if any)
missing_nodes <- setdiff(red_nodes, names(colored_node_hex_gain))
if(length(missing_nodes) > 0) {
  print(paste("Missing colors for nodes:", paste(missing_nodes, collapse = ", ")))
}

# Initialize edge color vector (darkgrey for all edges by default)
edge_colors <- rep("lightgrey", length(tree$edge[, 1]))

# Loop through the edges and color those that correspond to the blue nodes (loss branches)
for (i in 1:length(tree$edge[, 1])) {
  # Check if the target node (tree$edge[i, 2]) is in the list of blue_nodes (i.e., loss branches)
  if (tree$edge[i, 2] %in% red_nodes) {
    # Use the custom color from `colored_node_hex_loss` for the corresponding branch
    if (as.character(tree$edge[i, 2]) %in% names(colored_node_hex_gain)) {
      edge_colors[i] <- colored_node_hex_gain[as.character(tree$edge[i, 2])]
    } else {
      edge_colors[i] <- "lightgray"  # Default color if no match found
    }
  }
}

# Plot the tree with the custom colors for the gain branches
png("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/gaintree_plot.png", width = 2000, height = 2000, res = 300)
par(mar = c(2, 2, 2, 2)) 
plot(tree, type = "phylogram", edge.col = edge_colors, cex = 1.2, edge.width = 10, show.tip.label = TRUE)
dev.off() 

# Optionally, add node labels to visualize the internal node numbers
#nodelabels(node = 1:(length(tree$tip.label) + length(tree$edge) - 1), cex = 0.5, frame = "none", adj = c(1.5, 0))
 
library(grid)
library(gridExtra)

# Read the saved images
gain_genes <- rasterGrob(png::readPNG("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/plot_gain.png"), interpolate = TRUE)
gain_tree <- rasterGrob(png::readPNG("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/gaintree_plot.png"), interpolate = TRUE)


# Define layout matrix
layout <- matrix(c(1, 2,
                   1, NA),  # NA leaves an empty space
                 nrow = 2, byrow = TRUE)


# Arrange plots and save as PDF
pdf("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Plots/Gain_Genes.pdf", width = 14, height = 8) #Figure 3 in paper
grid.arrange(gain_genes, gain_tree, layout_matrix = layout, widths = c(3, 1), heights = c(1, 1))
dev.off()

