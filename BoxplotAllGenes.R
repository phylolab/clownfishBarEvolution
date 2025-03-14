# Load required libraries
library(ggplot2)
library(viridis)
library(dplyr)

setwd("/Users/lfitzger/Test_BayesCode_desktop")

dir.exists(file.path("/Users/lfitzger/Test_BayesCode_desktop"))

Loss <-read.csv(file = "Sig_Loss_contrasts_allgenes.csv", header = TRUE) #by branch
#Loss <-read.csv(file = "Sig_Loss_avercontrasts_allgenes.csv", header = TRUE) #by gene
Gain <-read.csv(file = "Sig_Gain_contrasts_allgenes.csv", header = TRUE) #by branch
#Gain <-read.csv(file = "Sig_Gain_avercontrasts_allgenes.csv", header = TRUE) #by gene

# Reshape the data to long format only use when plotting by gene NOT by branch 
#data_long <- tidyr::pivot_longer(Gain, cols = c(gain, loss, same), names_to = "transition", values_to = "contrast")

ggplot(Loss, aes(x = transition, y = contrast, fill = transition)) +
  geom_boxplot(show.legend = FALSE) +  # Add boxplot without legend
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add jittered points for individual data points
  labs(title = expression(paste("Average ", Delta*omega, " for Loss")),
       x = "Type",
       y = expression("Δω")) +
  theme_classic() + 
  coord_cartesian(ylim = c(0, 5)) +  # Set the y-axis limits to focus on the range between 0 and 1
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))  # Center and bold the title
