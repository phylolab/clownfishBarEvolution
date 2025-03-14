import pandas as pd
import os

# Read gene names from the text file
with open("/Users/lfitzger/Test_BayesCode_desktop/Sig_Gain_Discrete.txt", "r") as file:
    gene_names = file.read().splitlines()

# Initialize an empty list to store dataframes for each gene
all_gene_data = []

# Loop through each gene name
for gene_name in gene_names:
    # Read gene data
    gene_data = pd.read_csv(os.path.join("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states_root3/simapvsomega_colorgenes", f"{gene_name}.simap_vs_omega.tsv"), sep="\t")

    # Define a function to filter rows with highest weight within each group
    def filter_highest_weight(df):
        # Find the index of the row with the highest weight
        max_weight_idx = df["weights"].idxmax()
        
        # Keep the row with the highest weight
        return df.loc[max_weight_idx]

    # Apply the function to each group in the DataFrame
    filtered_gene_data = gene_data.loc[gene_data.groupby("node_name")["weights"].idxmax()].reset_index(drop=True)

    avg_contrasts = filtered_gene_data.groupby("transition").agg({"contrast": "mean"}).reset_index()

    # Add gene_name column to avg_contrasts
    avg_contrasts["gene_name"] = gene_name

    # Append avg_contrasts dataframe to the list
    all_gene_data.append(avg_contrasts)

# Concatenate all the dataframes in the list
combined_data = pd.concat(all_gene_data, ignore_index=True)

# Pivot the combined data
pivoted_data = combined_data.pivot_table(index='gene_name', columns='transition', values='contrast').reset_index()
pivoted_data.columns.name = None  # To remove the column name 'transition' from the header

# Print pivoted data (optional)
print(pivoted_data)

# Export pivoted data to a single CSV file
pivoted_data.to_csv("/Users/lfitzger/Test_BayesCode_desktop/Sig_Gain_avercontrasts_colorgenes.csv", index=False)
