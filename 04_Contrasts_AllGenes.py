import pandas as pd
import os

# Read gene names from the text file
with open("/Users/lfitzger/Test_BayesCode_desktop/Sig_Loss_Discrete.txt", "r") as file:
    gene_names = file.read().splitlines()

# Initialize an empty list to store dataframes for each gene
all_gene_data = []

# Loop through each gene name
for gene_name in gene_names:
    # Read gene data
    gene_data = pd.read_csv(os.path.join("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states_root3/simapvsomega_colorgenes", f"{gene_name}.simap_vs_omega.tsv"), sep="\t")

    # Group by 'node_name' and find the maximum weight within each group
    max_weights = gene_data.groupby('node_name')['weights'].transform('max')

    # Filter the original DataFrame to keep rows with maximum weights
    filtered_gene_data = gene_data[gene_data['weights'].isin(max_weights)].copy()

    # Modify the line where you add the gene_name column
    filtered_gene_data.loc[:, "gene_name"] = gene_name

    # Append 
    all_gene_data.append(filtered_gene_data)

# Concatenate all the dataframes in the list
combined_data = pd.concat(all_gene_data, ignore_index=True)

# Print pivoted data (optional)
print(combined_data)

# Export data to a single CSV file
combined_data.to_csv("/Users/lfitzger/Test_BayesCode_desktop/Sig_Loss_contrasts_colorgenes.csv", index=False)
