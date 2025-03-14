setwd("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/")
library(png)       # For reading PNG images
library(phytools)  # For plotting the tree
library(fishualize)

#Choose a fish species for the color palette
#option = "Epinephelus_lanceolatus"
fish_palette <- fish(n = 4, option = "Epinephelus_lanceolatus")  # You can change "Clownfish" to other species
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
#par(mar = c(6, 4, 2, 2), oma = c(2, 0, 0, 0))
par(mar = c(5, 4, 2, 2), oma = c(4, 0, 0, 0))  #try to add bottom space
# Plot the summary tree with adjusted layout
plot(summary(clownfish.maps), colors = map.cols, ftype = "i")
legend(x = "bottomleft",title= "Number of Vertical Bars", 
       legend = legend_labels, col = map.cols[ordered_levels], pt.cex = 3, 
       pch=16, cex=1.2, bty="n", ncol=2)
axisPhylo(backward = TRUE)

# Add a properly constrained x-axis title
mtext("MYA", side = 1, line = 3, cex = 1.2)

# Load tree structure to get tip labels and coordinates
image_folder <- "./nemo_vectors"
tip_images <- setNames(file.path(image_folder, paste0(clownfishtree$tip.label, ".png")),
                       clownfishtree$tip.label)

# Retrieve tree coordinates
tree_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Select five node positions (adjust these numbers based on print output)
gain_nodes <- c(28, 25, 52, 39, 6)

# Extract x and y coordinates for selected nodes
gain_x <- tree_coords$xx[gain_nodes]
gain_y <- tree_coords$yy[gain_nodes]
gain_x <- gain_x + -0.20  # Shift left
gain_y <- gain_y + 0.50  # Shift up

# Add "Gain" labels to the tree
text(x = gain_x, y = gain_y, labels = "G", col = "red", cex = 0.8, font = 2)

# Select five node positions (adjust these numbers based on print output)
loss_nodes <- c(23, 14, 40, 41, 2)

# Extract x and y coordinates for selected nodes
loss_x <- tree_coords$xx[loss_nodes]
loss_y <- tree_coords$yy[loss_nodes]
loss_x <- loss_x + -0.20  # Shift left
loss_y <- loss_y + 0.50  # Shift up

# Add "Loss" labels to the tree
text(x = loss_x, y = loss_y, labels = "L", col = "red", cex = 0.8, font = 2)

#nodelabels()
#tiplabels()

# Updated function to plot an image with correct aspect ratio
plot_img <- function(img, coords, scale) {
  # Calculate aspect ratio dynamically
  img_width <- ncol(img)
  img_height <- nrow(img)
  aspect_ratio <- img_width / img_height
  
  # Debugging: print the aspect ratio of each image
  cat(sprintf("Aspect ratio for image %s: %.2f (Width: %d, Height: %d)\n", tip_label, aspect_ratio, img_width, img_height))
  
  # Calculate scaled width and height
  scaled_height <- scale
  scaled_width <- scale * aspect_ratio

  
  # Overlay the image at specified coordinates
  rasterImage(img,
              xleft = coords[1] - scaled_width / 2,
              xright = coords[1] + scaled_width / 2,
              ybottom = coords[2] - scaled_width / 2,
              ytop = coords[2] + scaled_width / 2,
              interpolate = TRUE)
}

# Scaling factors for resizing images 
image_scale <- 0.4  # Reduce image size (adjust smaller if needed)

# Overlay images at tree tips
for (i in seq_along(clownfishtree$tip.label)) {
  tip_label <- clownfishtree$tip.label[i]
  
  # Check if the corresponding image exists
  if (file.exists(tip_images[tip_label])) {
    # Read the image
    img <- readPNG(tip_images[tip_label])
    
    # Tip coordinates
    x <- tree_coords$xx[i] + 1.98  # Adjusted x-coordinate for image placement
    y <- tree_coords$yy[i]
    
    # Plot the image at the tip
    plot_img(img, scale = image_scale, coords = c(x, y))
    
   } else {
    cat(sprintf("Image for %s not found: %s\n", tip_label, tip_images[tip_label]))
   }
}

# Final adjustments to the plot
axisPhylo(backward = TRUE)

