import pandas as pd
import numpy as np
from scipy.stats import chi2

# Load the three input files
file1 = "All_loss_genes_results1.csv"
file2 = "All_loss_genes_results2.csv"
file3 = "All_loss_genes_results3.csv"

#file1 = "All_gain_genes_results1.csv"
#file2 = "All_gain_genes_results2.csv"
#file3 = "All_gain_genes_results3.csv"


df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
df3 = pd.read_csv(file3)

# Add a column for the run
df1['run'] = 1
df2['run'] = 2
df3['run'] = 3

# Extract required columns, including foreground_w_H1
df1 = df1[['gene', 'logL_H0', 'logL_H1', 'foreground_w_H1', 'run']]
df2 = df2[['gene', 'logL_H0', 'logL_H1', 'foreground_w_H1', 'run']]
df3 = df3[['gene', 'logL_H0', 'logL_H1', 'foreground_w_H1', 'run']]

# Combine the data from all three runs
combined_df = pd.concat([df1, df2, df3], ignore_index=True)

# Initialize a list to store the final results
final_results = []

# Process each gene
for gene, group in combined_df.groupby('gene'):
    # Get logL_H0, logL_H1, and foreground_w_H1 for each run
    logL_H0_values = group.set_index('run')['logL_H0'].to_dict()
    logL_H1_values = group.set_index('run')['logL_H1'].to_dict()
    foreground_w_H1_values = group.set_index('run')['foreground_w_H1'].to_dict()

    # Identify the worst (highest) logL_H0 and logL_H1 values
    worst_logL_H0_run = max(logL_H0_values, key=logL_H0_values.get)
    worst_logL_H1_run = max(logL_H1_values, key=logL_H1_values.get)

    worst_logL_H0 = logL_H0_values[worst_logL_H0_run]
    worst_logL_H1 = logL_H1_values[worst_logL_H1_run]

    # Get the corresponding foreground_w_H1 for the best run
    best_foreground_w_H1 = foreground_w_H1_values.get(worst_logL_H1_run, None)

    # Calculate the LRT (Likelihood Ratio Test)
    LRT = 2 * (worst_logL_H1 - worst_logL_H0)

    # Calculate the p-value
    p_value = chi2.sf(LRT, df=1)

    # Create a row for the gene with all required values
    row = {
        'gene': gene,
        'logL_H0_run1': logL_H0_values.get(1, None),
        'logL_H0_run2': logL_H0_values.get(2, None),
        'logL_H0_run3': logL_H0_values.get(3, None),
        'best_logL_H0': worst_logL_H0,
        'best_run_H0': worst_logL_H0_run,
        'logL_H1_run1': logL_H1_values.get(1, None),
        'logL_H1_run2': logL_H1_values.get(2, None),
        'logL_H1_run3': logL_H1_values.get(3, None),  
        'best_logL_H1': worst_logL_H1,
        'best_run_H1': worst_logL_H1_run,
        'foreground_w_H1_run1': foreground_w_H1_values.get(1, None),
        'foreground_w_H1_run2': foreground_w_H1_values.get(2, None),
        'foreground_w_H1_run3': foreground_w_H1_values.get(3, None),
        'best_foreground_w_H1': best_foreground_w_H1,
        'LRT': LRT,
        'p_value': p_value
    }

    final_results.append(row)

# Convert the final results to a DataFrame
final_df = pd.DataFrame(final_results)

# Save the results to a new CSV file
output_file = "LogL_Comparison_LossResults_with_LRT_pvalues_and_foreground_w_H1.csv"
#output_file = "LogL_Comparison_GainResults_with_LRT_pvalues_and_foreground_w_H1.csv"
final_df.to_csv(output_file, index=False)

print(f"Comparison results saved to {output_file}")
