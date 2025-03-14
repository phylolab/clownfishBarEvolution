setwd("/Users/lfitzger/Test_BayesCode_desktop/scripts/DiscreteTraits")
getwd()
library(stats)

## My weighted regression output from the TransitionGainorLoss.py script for discrete traits ##
df <- read.csv("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/Color_genes_weighted_regression_results_slope50k.tsv", header=T, sep="\t")
head(df)

##########################################################################################
### Nicolas's test script to show the calculation behind the p.adjust formula 03.20.24 ###

pvalues <- c(runif(900), runif(100, 0, 0.05))

len_pval <- length(pvalues)

pvalue_rank <- rank(pvalues, ties.method = "last") #rank of pvalues
pvalue_over_rank <- pvalues * len_pval / pvalue_rank #transform pvalues to pvalue/rank

pvalues_adj <- vector(length = len_pval)

#calculate adjusted pvalues
for (i in 1:len_pval) {
  tmp_rank <- pvalue_rank[i]
  #For each p-value at rank i, the adjusted p-value is the minimum between
  #.   1 and the smallest ratio found by dividing the total number of hypotheses (m)
  #.   by the rank (j) (ie pvalue_over_rank) for all j greater than or equal to i
  #see BH controlling procedures in https://en.wikipedia.org/wiki/False_discovery_rate
  pvalues_adj[i] <- min(1, min(pvalue_over_rank[pvalue_rank>=tmp_rank]))
}

hist(pvalues)
hist(pvalues_adj)

pvalues_fct <- p.adjust(pvalues, method = "BH") #to check that we indeed get the same as the p.adjust function
hist(pvalues_fct)
#########################################################################################

# Check with my data ##

# Extract the genes and p-values for gain and loss separately
genes_gain <- df$Gene
pvalues_gain <- df$P.Value_Gain
genes_loss <- df$Gene
pvalues_loss <- df$P.Value_Loss

len_pval_gain <- length(pvalues_gain)
len_pval_loss <- length(pvalues_loss)

# Rank the p-values for gain and loss separately
pvalue_rank_gain <- rank(pvalues_gain, ties.method = "last")
pvalue_rank_loss <- rank(pvalues_loss, ties.method = "last")

# Transform p-values to pvalue/rank for gain and loss separately
pvalue_over_rank_gain <- pvalues_gain * len_pval_gain / pvalue_rank_gain
pvalue_over_rank_loss <- pvalues_loss * len_pval_loss / pvalue_rank_loss

# Initialize vectors for adjusted p-values and genes for gain and loss
pvalues_adj_gain <- vector(length = len_pval_gain)
pvalues_adj_loss <- vector(length = len_pval_loss)

genes_adj_gain <- vector(length = len_pval_gain)
genes_adj_loss <- vector(length = len_pval_loss)

# Calculate adjusted p-values for gain
for (i in 1:len_pval_gain) {
  tmp_rank <- pvalue_rank_gain[i]
  pvalues_adj_gain[i] <- min(1, min(pvalue_over_rank_gain[pvalue_rank_gain >= tmp_rank]))
  genes_adj_gain[i] <- genes_gain[i] # Store corresponding gene
}

# Calculate adjusted p-values for loss
for (i in 1:len_pval_loss) {
  tmp_rank <- pvalue_rank_loss[i]
  pvalues_adj_loss[i] <- min(1, min(pvalue_over_rank_loss[pvalue_rank_loss >= tmp_rank]))
  genes_adj_loss[i] <- genes_loss[i] # Store corresponding gene
}

# You can also use p.adjust function to compare with your calculated adjusted p-values
pvalues_fct_gain <- p.adjust(pvalues_gain, method = "BH")
pvalues_fct_loss <- p.adjust(pvalues_loss, method = "BH")

# Now 'genes_adj_gain' and 'genes_adj_loss' contain the corresponding genes
# for the adjusted p-values in 'pvalues_adj_gain' and 'pvalues_adj_loss' respectively.
hist(pvalues_adj_gain)
hist(pvalues_fct_gain)
hist(pvalues_adj_loss)
hist(pvalues_fct_loss)

# Create a data frame with genes and corresponding adjusted p-values for gain
df_adj_gain <- data.frame(Gene = genes_adj_gain, Adjusted_P_Value = pvalues_adj_gain)

# Create a data frame with genes and corresponding adjusted p-values for loss
df_adj_loss <- data.frame(Gene = genes_adj_loss, Adjusted_P_Value = pvalues_adj_loss)

# Sort the data frames by adjusted p-value
df_adj_gain <- df_adj_gain[order(df_adj_gain$Adjusted_P_Value), ]
df_adj_loss <- df_adj_loss[order(df_adj_loss$Adjusted_P_Value), ]

# Write data frames to CSV files
write.csv(df_adj_gain, "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/Color_genes_adjusted_pvalues_gain.csv", row.names = FALSE)
write.csv(df_adj_loss, "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/Color_genes_adjusted_pvalues_loss.csv", row.names = FALSE)

