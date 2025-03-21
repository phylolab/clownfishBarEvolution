##!/usr/bin/env python3

import os
import pandas as pd
import statsmodels.api as sm

# Define the path to the input file
folder_path = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states/simapvsomega_allgenes" #/Results/simapvsomega_allgenes.zip

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
        
        # Identify the transition category with the highest total weight
        max_weight_category = df.groupby('transition')['weights'].sum().idxmax()
        
        # Filter only rows that belong to this transition category
        df_filtered = df[df['transition'] == max_weight_category]
        
        # Compute the average omega for the filtered rows
        avg_omega = df_filtered['omega'].mean()
        
        # Change the transitions into factors (for the regression)
        df['transition'] = pd.Categorical(df['transition'])
        
        # Put the level "same" as the basis
        df['transition'] = df['transition'].cat.reorder_categories(['same', 'gain', 'loss'])
        
        # Run the regression with weights
        mdl_weight = sm.WLS.from_formula("contrast ~ transition", weights=df["weights"], data=df).fit()
        
        # Get the slope, t-statistics, and p-values for each effect
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
            'Max_Weight_Transition': [max_weight_category],  # Store the selected transition
            'Avg_Omega': [avg_omega],  # Average omega from selected transition
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
final_results_df[['Avg_Omega', 'Slope_Gain', 'T-Statistic_Gain', 'P-Value_Gain', 
                  'Slope_Loss', 'T-Statistic_Loss', 'P-Value_Loss']] = final_results_df[[
    'Avg_Omega', 'Slope_Gain', 'T-Statistic_Gain', 'P-Value_Gain', 
    'Slope_Loss', 'T-Statistic_Loss', 'P-Value_Loss']].astype(float)

# Output final results to a TSV file
output_filename = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states/All_genes_weighted_regression_results_slope50k_omega.tsv" #/Results/All_genes_weighted_regression_results_slope50k_omega.tsv
final_results_df.to_csv(output_filename, sep='\t', index=False)

print("Processing complete. Results saved to:", output_filename)
