setwd("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/")
library(png)       # For reading PNG images
library(phytools)  # For plotting the tree
library(fishualize)

fish_palette <- fish(n = 4, option = "Epinephelus_lanceolatus") 
fishualize(n = 4, option = "Epinephelus_lanceolatus")

stripes_bin <- read.csv("Data/adult_stripes_4states.csv", stringsAsFactors = F, row.names = 1)
clownfishstripes <- read.csv("Data/adult_stripes_binary1_4states.csv", stringsAsFactors = T, row.names = 1)

##  phylogenetic tree
clownfishtree <- read.tree("./Data/calibrated_tree.tre")
Stripes<-setNames(clownfishstripes[[1]], rownames(stripes_bin))

# Load the clownfish maps
clownfish.maps <- read.simmap("./Results/model0.6_clownfish4states.map", format="phylip", version=1)

# Define the correct legend order
legend_labels <- c("None", "One", "Two", "Three")

# Reorder the color mapping based on the desired order
ordered_levels <- c("None", "One", "Two", "Three")  # New order

# Define map colors using fishualize
map.cols <- setNames(fish_palette, levels(Stripes))

# Add space for the axis and legend
par(mar = c(6, 4, 2, 2), oma = c(2, 0, 0, 0))

# Plot the summary tree with adjusted layout
plot(summary(clownfish.maps), colors = map.cols, ftype = "i")
legend(x = "bottomleft",title= "Number of Vertical Bars", 
       legend = legend_labels, col = map.cols[ordered_levels], pt.cex = 3, 
       pch=16, cex=1.2, bty="n", ncol=2)
axisPhylo(backward = TRUE)


gain_nodes <- c(28, 25, 52, 39, 6)

# Extract x and y coordinates for selected nodes
gain_x <- tree_coords$xx[gain_nodes]
gain_y <- tree_coords$yy[gain_nodes]
gain_x <- gain_x + -0.20  # Shift left
gain_y <- gain_y + 0.58  # Shift up

# Add "Gain" labels to the tree
text(x = gain_x, y = gain_y, labels = "G", col = "red", cex = 0.8, font = 2)


loss_nodes <- c(23, 14, 40, 41, 2)

# Extract x and y coordinates for selected nodes
loss_x <- tree_coords$xx[loss_nodes]
loss_y <- tree_coords$yy[loss_nodes]
loss_x <- loss_x + -0.20  # Shift left
loss_y <- loss_y + 0.58  # Shift up

# Add "Loss" labels to the tree
text(x = loss_x, y = loss_y, labels = "L", col = "purple", cex = 0.8, font = 2)

#nodelabels()
#tiplabels()

same_nodes <- c(30, 31, 32, 33, 34, 35, 36, 37, 38,
                42, 45, 47, 50, 51, 49, 48, 46, 44, 43,
                53, 54, 55, 27, 26, 24, 22, 21, 20,
                19, 18, 17, 16, 15, 13, 12, 11, 10, 9, 
                8, 7, 5, 4, 3, 1)

# Extract x and y coordinates for selected nodes
same_x <- tree_coords$xx[same_nodes]
same_y <- tree_coords$yy[same_nodes]
same_x <- same_x + -0.20  # Shift left
same_y <- same_y + 0.58  # Shift up

# Add "Same" labels to the tree
text(x = same_x, y = same_y, labels = "S", col = "gray30", cex = 0.8, font = 2)

##########################################################
######  Add Tree2 internal node labels to Tree 1 #########
##########################################################
setwd("/Users/lfitzger/Test_BayesCode/")

# Load the ape package
library(ape)

# Read tree files
tree1 <- read.tree("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/Data/calibrated_tree.tre")
tree2 <- read.tree("/Users/lfitzger/Test_BayesCode/Clown.DatedPhylogeny.annotated.tree")

# Extract internal node labels from tree2
node_labels_tree2 <- tree2$node.label

# Assign those labels to tree1's internal nodes (tree1 must have the same number of internal nodes as tree2)
tree1$node.label <- node_labels_tree2

# Plot the tree without node labels and tip labels
plot(tree1, show.node.label = FALSE, show.tip.label = TRUE)

# Add node labels with a smaller font size
nodelabels(tree1$node.label, cex = 0.6)

