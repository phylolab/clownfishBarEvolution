setwd("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/CodeML")
getwd()
library(stats)

# Read your data
df <- read.csv("/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/CodeML/LogL_Comparison_LossResults_with_LRT_pvalues_and_foreground_w_H1.csv", header=T)

# Inspect data to ensure it's correct
str(df)
summary(df)

# Remove any rows with NA p-values
df <- df[!is.na(df$p_value), ]

# Extract the chromosome and gene information
#chromosome <- df$gene  # Assuming your column is named 'chromosome'
genes <- df$gene

# Strip the unwanted suffix from gene names (e.g., ".Loss.Vert.Bar")
# Using sub to remove ".Loss.Vert.Bar" or any suffix after ".t1" or ".t2"
genes_cleaned <- sub("\\.Gain\\.Vert\\.Bar$", "", genes)
genes_cleaned <- sub("\\.t1.*$", ".t1", genes_cleaned)
genes_cleaned <- sub("\\.t2.*$", ".t2", genes_cleaned)

# Concatenate chromosome and cleaned gene into a single identifier
#chr_gene <- paste( genes_cleaned, sep="_")

# Extract the p-values
pvalues <- df$p_value

#extract the foreground_w
foregroundomega <- df$best_foreground_w_H1

# Apply BH adjustment (Benjamini-Hochberg)
pvalues_adj <- p.adjust(pvalues, method = "BH")

# Compare the original and adjusted p-values visually
hist(pvalues, main="Original P-values", xlab="P-values")
hist(pvalues_adj, main="Adjusted P-values (BH)", xlab="Adjusted P-values")

# Create a data frame with combined chromosome and gene, and corresponding adjusted p-values
df_adj <- data.frame(Chr_Gene = genes_cleaned, foregroundomega, Adjusted_P_Value = pvalues_adj)

# Sort the data frame by adjusted p-value
df_adj <- df_adj[order(df_adj$Adjusted_P_Value), ]

# Write the result to a CSV file
write.csv(df_adj, "/Users/lfitzger/Desktop/PhD/Chapter2_ClownfishStripes/CodeML/All_loss_genes_foregroundomega.csv", row.names = FALSE)

