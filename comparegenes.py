import pandas as pd

# Define the paths to the input files
file_50k = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states/All_genes_weighted_regression_results_slope50k.tsv"
file_2k = "/Users/lfitzger/Test_BayesCode_desktop/Old_Results_DiscreteTraits/model0.6_maps4states/All_genes_weighted_regression_results_slope2k.tsv"
genes_file = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states/Sanity_check/gain_diffgenes.txt"

# Read the list of genes from the text file
with open(genes_file, 'r') as f:
    genes_to_compare = f.read().splitlines()

# Read the files into DataFrames
df_50k = pd.read_csv(file_50k, sep='\t')
df_2k = pd.read_csv(file_2k, sep='\t')

# Filter the DataFrames to include only the specified genes
df_50k_filtered = df_50k[df_50k['Gene'].isin(genes_to_compare)]
df_2k_filtered = df_2k[df_2k['Gene'].isin(genes_to_compare)]

# Merge the filtered DataFrames on the Gene column
merged_df = pd.merge(df_50k_filtered, df_2k_filtered, on='Gene', suffixes=('_50k', '_2k'))

# Calculate the differences between the slopes, t-statistics, and p-values
merged_df['Slope_Diff_Gain'] = merged_df['Slope_Gain_50k'] - merged_df['Slope_Gain_2k']
merged_df['Slope_Diff_Loss'] = merged_df['Slope_Loss_50k'] - merged_df['Slope_Loss_2k']
merged_df['T-Statistic_Diff_Gain'] = merged_df['T-Statistic_Gain_50k'] - merged_df['T-Statistic_Gain_2k']
merged_df['T-Statistic_Diff_Loss'] = merged_df['T-Statistic_Loss_50k'] - merged_df['T-Statistic_Loss_2k']
merged_df['P-Value_Diff_Gain'] = merged_df['P-Value_Gain_50k'] - merged_df['P-Value_Gain_2k']
merged_df['P-Value_Diff_Loss'] = merged_df['P-Value_Loss_50k'] - merged_df['P-Value_Loss_2k']

# Output the merged results with differences to a new TSV file
output_filename = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states/Sanity_check/Comparisonresults_newsiggenes_gain.tsv"
merged_df.to_csv(output_filename, sep='\t', index=False)
