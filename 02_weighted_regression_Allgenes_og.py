#!/usr/bin/env python3

import os
import pandas as pd
import statsmodels.api as sm

# Define the path to the input file
folder_path = "/Users/lfitzger/Test_BayesCode_desktop/Old_Results_DiscreteTraits/model0.6_maps4states/simapvsomega_colorgenes"

# Initialize an empty list to store results for all genes
all_results = []

# Loop through all files in the folder
for file_name in os.listdir(folder_path):
    # Check if the file is a TSV file
    if file_name.endswith(".tsv"):
        # Construct the full path to the file
        file_path = os.path.join(folder_path, file_name)
        
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep="\t")
        
        # Change the transitions into factors (for the regression)
        df['transition'] = pd.Categorical(df['transition'])
        
        # Put the level "same" as the basis
        df['transition'] = df['transition'].cat.reorder_categories(['same', 'gain', 'loss'])
        
        # Run the regression with weights (assuming you have weights in the dataframe)
        mdl_weight = sm.WLS.from_formula("contrast ~ transition", weights=df["weights"], data=df).fit()
        
        # Get the t-statistics and p-values for each effect (first gain, then loss) 
        tstat = mdl_weight.tvalues["transition[T.gain]"], mdl_weight.tvalues["transition[T.loss]"]
        pval = mdl_weight.pvalues["transition[T.gain]"], mdl_weight.pvalues["transition[T.loss]"]
        
        # Extract the gene name from the input file name
        gene_name = '.'.join(file_name.split('.')[:2])
        
        # Create a DataFrame to store results for the current gene
        gene_results_df = pd.DataFrame({
            'Gene': [gene_name],
            'T-Statistic_Gain': [tstat[0]],
            'P-Value_Gain': [pval[0]],
            'T-Statistic_Loss': [tstat[1]],
            'P-Value_Loss': [pval[1]]
        })
        
        # Append the results for the current gene to the list of all results
        all_results.append(gene_results_df)

# Concatenate results for all genes into a single DataFrame
final_results_df = pd.concat(all_results, ignore_index=True)

# Convert numeric columns to appropriate data types
final_results_df[['T-Statistic_Gain', 'P-Value_Gain', 'T-Statistic_Loss', 'P-Value_Loss']] = final_results_df[[
    'T-Statistic_Gain', 'P-Value_Gain', 'T-Statistic_Loss', 'P-Value_Loss']].astype(float)

# Output final results to a TSV file
output_filename = "/Users/lfitzger/Test_BayesCode_desktop/Old_Results_DiscreteTraits/model0.6_maps4states/Color_genes_weighted_regression_results.tsv"
final_results_df.to_csv(output_filename, sep='\t', index=False)
