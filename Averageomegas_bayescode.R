# Load necessary library
library(ggplot2)

# Define file path
avg_omega_file <- "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/Avg_Omega_per_Gene_with_Significance.csv"

# Load the dataset
avg_omega_df <- read.csv(avg_omega_file, header=TRUE)

# Load the dataset
avg_omega_df <- read.csv(avg_omega_file, header=TRUE)

# Ensure Significance is a factor with all four categories
avg_omega_df$Significance <- factor(avg_omega_df$Significance, 
                                    levels = c("Significant_Gain", "Significant_Loss", 
                                               "Significant_Gain_Loss", "Insignificant"))

# Compute summary statistics per category
summary_stats <- aggregate(Avg_Omega ~ Significance, data = avg_omega_df, 
                           FUN = function(x) c(mean=mean(x), min=min(x), max=max(x)))

# Convert list-columns into separate columns
summary_stats <- do.call(data.frame, summary_stats)

# Rename columns properly
colnames(summary_stats) <- c("Significance", "Mean_Omega", "Min_Omega", "Max_Omega")

# Print the summary statistics in the console
print(summary_stats)

# Print unique categories to confirm they are recognized
print(unique(avg_omega_df$Significance))

# Remove extreme outliers (if necessary)
avg_omega_df <- avg_omega_df[avg_omega_df$Avg_Omega < 3, ]  # Adjust threshold if needed

# Compute summary statistics per category
summary_stats <- aggregate(Avg_Omega ~ Significance, data = avg_omega_df, 
                           FUN = function(x) c(mean=mean(x), min=min(x), max=max(x)))

# Convert list-columns into separate columns
summary_stats <- do.call(data.frame, summary_stats)

# Rename columns properly
colnames(summary_stats) <- c("Significance", "Mean_Omega", "Min_Omega", "Max_Omega")

# Print the summary statistics in the console
print(summary_stats)
# Create the violin plot with Omega max = 3
ggplot(avg_omega_df, aes(x=Significance, y=Avg_Omega, fill=Significance)) +
  geom_violin(trim=FALSE, alpha=0.7) +  # Violin plot with full distribution
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  #geom_jitter(width=0.2, size=1, alpha=0.5)  # Adds individual points# Add a boxplot inside for summary stats
  theme_minimal() +
  labs(title="Distribution of Average Omega by Significance",
       x="Significance Category",
       y="Average Omega") +
  scale_fill_manual(values=c("Significant_Gain"="red", 
                             "Significant_Loss"="blue",
                             "Significant_Gain_Loss"="purple",
                             "Insignificant"="gray")) +
  coord_cartesian(ylim=c(0, max(avg_omega_df$Avg_Omega, na.rm=TRUE)))  # Proper y-axis scaling


# adjusts the violin width for the number of points 
ggplot(avg_omega_df, aes(x=Significance, y=Avg_Omega, fill=Significance)) +
  geom_violin(trim=FALSE, scale="count", alpha=0.7) +  # Adjust width based on sample size
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +  # Add boxplot for summary stats
  theme_minimal() +
  labs(title="Distribution of Average Omega by Significance",
       x="Significance Category",
       y="Average Omega") +
  scale_fill_manual(values=c("Significant_Gain"="red", 
                             "Significant_Loss"="blue",
                             "Significant_Gain_Loss"="purple",
                             "Insignificant"="gray")) +
  coord_cartesian(ylim=c(0, max(avg_omega_df$Avg_Omega, na.rm=TRUE)))  # Proper y-axis scaling

#Add points to show data 
ggplot(avg_omega_df, aes(x=Significance, y=Avg_Omega, fill=Significance)) +
  geom_violin(trim=FALSE, scale="width", alpha=0.7) +  # Keeps violin width equal
  geom_jitter(width=0.2, size=1, alpha=0.5) +  # Adds individual points
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +  
  theme_minimal() +
  labs(title="Distribution of Average Omega by Significance",
       x="Significance Category",
       y="Average Omega") +
  scale_fill_manual(values=c("Significant_Gain"="red", 
                             "Significant_Loss"="blue",
                             "Significant_Gain_Loss"="purple",
                             "Insignificant"="gray")) +
  coord_cartesian(ylim=c(0, max(avg_omega_df$Avg_Omega, na.rm=TRUE)))


#If "Insignificant" has too many points, you can randomly select 500–1000 of them:
set.seed(123)  # For reproducibility
insignificant_subset <- avg_omega_df[avg_omega_df$Significance == "Insignificant", ]
insignificant_sample <- insignificant_subset[sample(nrow(insignificant_subset), 1000), ]

# Keep all significant genes and only sampled insignificant ones
balanced_data <- rbind(avg_omega_df[avg_omega_df$Significance != "Insignificant", ], insignificant_sample)

# Re-run the plot with this smaller dataset
ggplot(balanced_data, aes(x=Significance, y=Avg_Omega, fill=Significance)) +
  geom_violin(trim=FALSE, scale="width", alpha=0.7) +  
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +  
  theme_minimal() +
  labs(title="Balanced Violin Plot of Average Omega",
       x="Significance Category",
       y="Average Omega") +
  scale_fill_manual(values=c("Significant_Gain"="red", 
                             "Significant_Loss"="blue",
                             "Significant_Gain_Loss"="purple",
                             "Insignificant"="gray")) +
  coord_cartesian(ylim=c(0, max(balanced_data$Avg_Omega, na.rm=TRUE)))

