#!/usr/bin/env python3

## with slopes ##

import os
import pandas as pd
import statsmodels.api as sm

# Define the path to the input file
folder_path = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/simapvsomega_allgenes"

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
        
        # Get the slope, t-statistics and p-values for each effect (first gain, then loss)
        slope_gain = mdl_weight.params["transition[T.gain]"]
        slope_loss = mdl_weight.params["transition[T.loss]"]
        tstat_gain = mdl_weight.tvalues["transition[T.gain]"]
        tstat_loss = mdl_weight.tvalues["transition[T.loss]"]
        pval_gain = mdl_weight.pvalues["transition[T.gain]"]
        pval_loss = mdl_weight.pvalues["transition[T.loss]"]
        
        # Extract the gene name from the input file name
        gene_name = '.'.join(file_name.split('.')[:2])
        
        # Create a DataFrame to store results for the current gene
        gene_results_df = pd.DataFrame({
            'Gene': [gene_name],
            'Slope_Gain': [slope_gain],
            'T-Statistic_Gain': [tstat_gain],
            'P-Value_Gain': [pval_gain],
            'Slope_Loss': [slope_loss],
            'T-Statistic_Loss': [tstat_loss],
            'P-Value_Loss': [pval_loss]
        })
        
        # Append the results for the current gene to the list of all results
        all_results.append(gene_results_df)

# Concatenate results for all genes into a single DataFrame
final_results_df = pd.concat(all_results, ignore_index=True)

# Convert numeric columns to appropriate data types
final_results_df[['Slope_Gain', 'T-Statistic_Gain', 'P-Value_Gain', 'Slope_Loss', 'T-Statistic_Loss', 'P-Value_Loss']] = final_results_df[[
    'Slope_Gain', 'T-Statistic_Gain', 'P-Value_Gain', 'Slope_Loss', 'T-Statistic_Loss', 'P-Value_Loss']].astype(float)

# Output final results to a TSV file
output_filename = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/All_genes_weighted_regression_results_slope50k.tsv"
final_results_df.to_csv(output_filename, sep='\t', index=False)
