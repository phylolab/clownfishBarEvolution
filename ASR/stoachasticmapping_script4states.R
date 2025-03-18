#Script examining number of stripes with stochastic mapping and ancestral state
#reconstruction on the most updated tree

setwd("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/")
library(phytools)
#install.packages('corHMM')
library(corHMM)

stripes_bin <- read.csv("Data/adult_stripes_4states.csv", stringsAsFactors = F, row.names = 1)
stripes_bin <- data.frame(species = rownames(stripes_bin), stripes_bin)
clownfishstripes <- read.csv("Data/adult_stripes_binary1_4states.csv", stringsAsFactors = T, row.names = 1)
head (clownfishstripes)

##  phylogenetic tree
clownfishtree <- read.tree("./Data/calibrated_tree.tre")
nodelabels()

## Phytools tutorial for plotting stripes on the tree
Stripes<-setNames(clownfishstripes[[1]],
                  rownames(stripes_bin))

plotTree(clownfishtree,show.tip.label = TRUE)
nodelabels()
X<-strsplit(setNames(as.character(Stripes),names(Stripes)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(clownfishtree),4,
             dimnames=list(clownfishtree$tip.label,
                           c("None","One","Two", "Three")))
for(i in 1:Ntip(clownfishtree)) 
  pies[clownfishtree$tip.label[i],X[[clownfishtree$tip.label[i]]]]<-
        rep(1/length(X[[clownfishtree$tip.label[i]]]),
        length(X[[clownfishtree$tip.label[i]]]))
cols<-c("black","red","green","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("None","One","Two","Three"),
       pt.cex=2,pch=16,col=cols,  title = "Number of Stripes",bty = "n")


#Make own matrix for testing 
st_rates <- matrix(0,nrow=4, ncol=4)
dimnames(st_rates) <-list(levels(Stripes), levels(Stripes))

##model0.1 hidden states = 0 ARD, when every rate is different 
rates0.1 <- st_rates
rates0.1[,]  <- c(0,1:4,0,5:8,0,9:12,0) 
rates0.1[,]
model0.1 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.1, root.p=NULL)

##model0.2 hidden states = 0 SYM, where forward and reverse are the same
rates0.2 <- st_rates
rates0.2[,]
rates0.2[1,2]<-rates0.2[1,2]+1
rates0.2[2,1]<-rates0.2[2,1]+1
rates0.2[1,3]<-rates0.2[1,3]+2
rates0.2[3,1]<-rates0.2[3,1]+2
rates0.2[1,4]<-rates0.2[1,4]+3
rates0.2[4,1]<-rates0.2[4,1]+3
rates0.2[2,3]<-rates0.2[2,3]+4
rates0.2[3,2]<-rates0.2[3,2]+4
rates0.2[2,4]<-rates0.2[2,4]+5
rates0.2[4,2]<-rates0.2[4,2]+5
rates0.2[3,4]<-rates0.2[3,4]+6
rates0.2[4,3]<-rates0.2[4,3]+6
rates0.2[,]
model0.2 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.2, root.p=NULL)

##model0.3 hidden states = 0 ER, single parameter governs all (everything except diagonal, 1 rate)
rates0.3 <- st_rates
rates0.3[,]
rates0.3[1,2]<-rates0.3[1,2]+1
rates0.3[2,1]<-rates0.3[2,1]+1
rates0.3[1,3]<-rates0.3[1,3]+1
rates0.3[3,1]<-rates0.3[3,1]+1
rates0.3[1,4]<-rates0.3[1,4]+1
rates0.3[4,1]<-rates0.3[4,1]+1
rates0.3[2,3]<-rates0.3[2,3]+1
rates0.3[3,2]<-rates0.3[3,2]+1
rates0.3[2,4]<-rates0.3[2,4]+1
rates0.3[4,2]<-rates0.3[4,2]+1
rates0.3[3,4]<-rates0.3[3,4]+1
rates0.3[4,3]<-rates0.3[4,3]+1
rates0.3[,]
model0.3 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.3, root.p=NULL)

##model 0.4 biological model - bi directional, SYM 
rates0.4 <- st_rates
rates0.4[,]
rates0.4[1,2]<-rates0.4[1,2]+1
rates0.4[2,1]<-rates0.4[2,1]+1
rates0.4[2,4]<-rates0.4[2,4]+5
rates0.4[3,4]<-rates0.4[3,4]+6
rates0.4[4,2]<-rates0.4[4,2]+5
rates0.4[4,3]<-rates0.4[4,3]+6
rates0.4[,]
model0.4 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.4, root.p=NULL)

##model 0.5 biological model - ARD
rates0.5 <- st_rates
rates0.5[,]
rates0.5[1,2]<-rates0.5[1,2]+4
rates0.5[2,1]<-rates0.5[2,1]+1
rates0.5[2,4]<-rates0.5[2,4]+11
rates0.5[3,4]<-rates0.5[3,4]+12
rates0.5[4,2]<-rates0.5[4,2]+6
rates0.5[4,3]<-rates0.5[4,3]+9
rates0.5[,]
model0.5 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.5, root.p=NULL)

##model 0.6 biological model - ER 
rates0.6 <- st_rates
rates0.6[,]
rates0.6[1,2]<-rates0.6[1,2]+1
rates0.6[2,1]<-rates0.6[2,1]+1
rates0.6[2,4]<-rates0.6[2,4]+1
rates0.6[3,4]<-rates0.6[3,4]+1
rates0.6[4,2]<-rates0.6[4,2]+1
rates0.6[4,3]<-rates0.6[4,3]+1
rates0.6[,]
model0.6 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.6, root.p=NULL)

##model 0.8 biological model - ARD bi-directional
rates0.8 <- st_rates
rates0.8[,]
rates0.8[1,2]<-rates0.8[1,2]+1
rates0.8[2,1]<-rates0.8[2,1]+1
rates0.8[2,4]<-rates0.8[2,4]+2
rates0.8[3,4]<-rates0.8[3,4]+3
rates0.8[4,2]<-rates0.8[4,2]+2
rates0.8[4,3]<-rates0.8[4,3]+3
rates0.8[,]
model0.8 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.8, root.p=NULL)

##model 0.9 biological model - ARD uni-directional
rates0.9 <- st_rates
rates0.9[,]
rates0.9[1,2]<-rates0.9[1,2]+2
rates0.9[2,1]<-rates0.9[2,1]+1
rates0.9[2,4]<-rates0.9[2,4]+5
rates0.9[3,4]<-rates0.9[3,4]+6
rates0.9[4,2]<-rates0.9[4,2]+3
rates0.9[4,3]<-rates0.9[4,3]+4
rates0.9[,]
model0.9 <- corHMM(phy = clownfishtree, data = stripes_bin, rate.cat = 1, rate.mat = rates0.9, root.p=NULL)


#compare models
model0.1$AIC
model0.2$AIC
model0.3$AIC
model0.4$AIC
model0.5$AIC
model0.6$AIC
model0.8$AIC
model0.9$AIC


#FitMK models
rates0.1[is.na(rates0.1)] = 0
fit.ARD<- fitMk(clownfishtree, Stripes,rates0.1)
rates0.2[is.na(rates0.2)] = 0
fit.SYM<- fitMk(clownfishtree, Stripes,rates0.2)
rates0.3[is.na(rates0.3)] = 0
fit.ER<- fitMk(clownfishtree, Stripes,rates0.3)
rates0.4[is.na(rates0.4)] = 0
fit.BIOSYM<- fitMk(clownfishtree, Stripes,rates0.4)
rates0.5[is.na(rates0.5)] = 0
fit.BIOARD<- fitMk(clownfishtree, Stripes,rates0.5)
rates0.6[is.na(rates0.6)] = 0
fit.BIOER<- fitMk(clownfishtree, Stripes,rates0.6)
fit.BIOER.pi<- fitMk(clownfishtree, Stripes,rates0.6, pi = c(0,0,1,0))
rates0.8[is.na(rates0.8)] = 0
fit.BIOARD2<- fitMk(clownfishtree, Stripes,rates0.8)
rates0.9[is.na(rates0.9)] = 0
fit.BIOARD3<- fitMk(clownfishtree, Stripes,rates0.9)

AIC(fit.ARD,fit.SYM,fit.ER,fit.BIOSYM,fit.BIOARD,fit.BIOER,fit.BIOARD2,fit.BIOARD3)
########################################################################################
clownfish.ordered2<-fitMk(clownfishtree, Stripes,rates0.6)
Q<-as.Qmatrix(clownfish.ordered2)
dim(Q)
print(Q,digits=3)

X<-clownfish.ordered2$data
head(X)
df<-read.csv("Data/adult_stripes_4states.csv",stringsAsFactors = T, row.names = 1, check.names=FALSE)
mymatrix <- as.matrix(df)
dim(mymatrix)
head(mymatrix)
#pi = constrains the root node to 3
#clownfish.maps<-make.simmap(clownfishtree,x=mymatrix,Q=Q,nsim=100,pi = c(0,0,1,0)) 
#pi = constrains the root node to 2
#clownfish.maps<-make.simmap(clownfishtree,x=mymatrix,Q=Q,nsim=100,pi = c(0,0,0,1)) 
clownfish.maps<-make.simmap(clownfishtree,x=mymatrix,Q=Q,nsim=100) 
clownfish.maps
########################################################################################
#saveRDS(list(simmap = clownfish.maps), file = "Results/model0.6_adult_stripes_4states.rds")

#write.simmap(clownfish.maps, file="Results/model0.6_1clownfish4states.map", append=FALSE, map.order=NULL)
#to write all the 100 simmap files for stripe association
#lapply(clownfish.maps,write.simmap,file="Results/model0.6_clownfish4states.map",append=TRUE) 
#########################################################################################
# If you already have a saved simmap, use read.simmap
file.exists("./Results/model0.6_clownfish4states.map")
clownfish.maps <- read.simmap("./Results/model0.6_clownfish4states.map",format="phylip",version=1)

#Extracting the probability at the root 
# Check if the file exists
if (file.exists("./Results/model0.6_clownfish4states.map")) {
  
  # Read the stochastic mapping results
  clownfish.maps <- read.simmap("./Results/model0.6_clownfish4states.map", format = "phylip", version = 1)
  
  # Get summary statistics from the simmap object
  simmap_summary <- summary(clownfish.maps)
  
  # Extract root probabilities
  root_probs <- simmap_summary$ace[1,]  # The first row corresponds to the root node
  
  # Print root probabilities
  print(root_probs)
} else {
  cat("File not found: ./Results/model0.6_clownfish4states.map\n")
}

#########################################################################################
map.cols<-setNames(colorRampPalette(c("black","red","blue","green"))(4),
                   levels(Stripes))

# Add space for the axis and legend
par(mar = c(6, 4, 2, 2), oma = c(2, 0, 0, 0)) 

# Plot with adjusted layout
plot(summary(clownfish.maps), colors = map.cols, ftype = "i")
legend(x = "bottomleft", legend = levels(Stripes), col = map.cols, pt.cex = 1.2, pch=15, cex=0.8, bty="n", ncol=2)
axisPhylo(backward = TRUE)

####################################################################################
##plotting nemo vectors at the tips 11.19.24
# Load required library
# Load necessary libraries
library(png)       # For reading PNG images
library(phytools)  # For plotting the tree

# Define the folder where images are stored
image_folder <- "./nemo_vectors"

# Match each tip label with its corresponding image file
tip_images <- setNames(file.path(image_folder, paste0(clownfishtree$tip.label, ".png")),
                       clownfishtree$tip.label)

# Set label size for the plot
label_cex <- 2

# Plot the tree to get its coordinates
par(mar = c(8, 6, 4, 10), oma = c(2, 0.5, 0.5, 2))  # Adjust margins for space
plot(summary(clownfish.maps), colors = map.cols, ftype = "i", label.cex = label_cex)

# Retrieve tree coordinates after plotting
tree_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Add a legend if needed
legend(x = "bottomleft", legend = levels(Stripes), col = map.cols, pt.cex = 2, pch = 16, cex = 0.8, bty = "n", title = "Number of Vertical Bars", ncol = 2)

# Add the phylogenetic axis
axisPhylo(backward = TRUE)

# Define the fixed horizontal offset (3 inches to the right)
horizontal_offset_inch <- 3.34  # Horizontal offset in inches
x_range <- diff(par()$usr[1:2])  # Get the x-axis range
inch_to_plot_units <- horizontal_offset_inch / (x_range / 8)  # Convert inches to plot units (assuming 8 inches for the x-range)

# Overlay images at tips, aligned vertically in the same position (no variation due to label length)
for (i in seq_along(clownfishtree$tip.label)) {
  tip_label <- clownfishtree$tip.label[i]
  
  # Check if the corresponding image exists
  if (file.exists(tip_images[tip_label])) {
    # Read the image
    img <- readPNG(tip_images[tip_label])
    
    # Convert transparent areas to NA for correct rendering
    transparent <- img[,,4] == 0
    img <- as.raster(img[,,1:3])  # Remove alpha channel
    img[transparent] <- NA  # Set transparent pixels to NA
    
    # Tip coordinates (fixed x position with offset, keep y position as is)
    x <- tree_coords$xx[i] + inch_to_plot_units  # Add the horizontal offset (in plot units)
    y <- tree_coords$yy[i]  # Keep y position as it is (same as tip)
    
    # Get the image's original aspect ratio
    img_width_orig <- ncol(img)
    img_height_orig <- nrow(img)
    aspect_ratio <- img_width_orig / img_height_orig
    
    # Adjust image size dynamically based on tree coordinates
    y_range <- diff(par()$usr[3:4])  # Get y-axis range
    
    # Set the image width relative to the x-axis range
    img_width <- 0.05 * x_range  # Image width proportional to x-axis range
    img_height <- img_width / aspect_ratio  # Maintain aspect ratio
    
    # Adjust image height to prevent squishing
    img_height <- img_height * 1.8  # Allow more space vertically (scale by 1.8 times)
    
    # Ensure that the image height does not exceed the y-axis range
    if (img_height > y_range * 0.35) {
      img_height <- y_range * 0.35
      img_width <- img_height * aspect_ratio
    }
    
    # Overlay the image at the fixed x position and the original y position (no change in y)
    rasterImage(img, 
                xleft = x - img_width / 2, 
                xright = x + img_width / 2, 
                ybottom = y - img_height / 2, 
                ytop = y + img_height / 2, 
                interpolate = TRUE)
  } else {
    cat(sprintf("Image for %s not found: %s\n", tip_label, tip_images[tip_label]))
  }
}


####################################################################################

##one stochastic map
plot(clownfish.maps[[10]],map.cols,ftype="i",
     ylim=c(0,1.05*Ntip(clownfishtree)))
legend(x="bottomleft",legend=levels(Stripes), col=map.cols,pt.cex=1.2,pch=15,cex=0.8, bty="n",ncol=2)
obj<-summary(clownfish.maps)
par(fg="black")
nodelabels(pie=obj$ace[1:clownfishtree$Nnode,],piecol=map.cols,
           cex=0.3)

#####################################################################################
## From Theo's script
#Plot the nemo vector images on the tips of the tree
plot_img <- function(img, fh, coords, col.rep = NULL){
  transparent <- img[,,4] == 0
  img <- as.raster(img[,,1:3])
  img[transparent] <- NA
  # image dimensions
  img_dim <- dim(img)[1:2]
  fig.dim <- par()$fin#/c(diff(par()$fig[1:2]), diff(par()$fig[3:4]))
  if(!is.null(col.rep)){
    img <- as.raster(img)
    img[img != "#FFFFFF00"] <- col.rep
  }
  rasterImage(img, xleft = coords[1]-(fh*img_dim[2]/fig.dim[1])/2, ybottom = coords[2]-(fh*img_dim[1]/fig.dim[2])/2, xright = coords[1]+(fh*img_dim[2]/fig.dim[1])/2, ytop = coords[2]+(fh*img_dim[1]/fig.dim[2])/2, interpolate = F)
}
2:23
#read the image with readPNG()
2:25
md <- get('last_plot.phylo', envir = .PlotPhyloEnv)
#To get the coordinates of tips and nodes after plotting a tree

