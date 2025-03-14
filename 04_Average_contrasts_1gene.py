import pandas as pd
import os

# Read gene names from the text file
#with open("/Users/lfitzger/Test_BayesCode_desktop/Sig_Loss_Discrete.txt", "r") as file:
#    gene_names = file.read().splitlines()
gene_name = "chr06_g3597.t1"
gene_data = pd.read_csv(os.path.join("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/simapvsomega_allgenes", f"{gene_name}.simap_vs_omega.tsv"), sep="\t")

# Define a function to filter rows with highest weight within each group
def filter_highest_weight(df):
    # Find the index of the row with the highest weight
    max_weight_idx = df["weights"].idxmax()
    
    # Keep the row with the highest weight
    return df.loc[max_weight_idx]

# Apply the function to each group in the DataFrame
filtered_gene_data = gene_data.loc[gene_data.groupby("node_name")["weights"].idxmax()].reset_index(drop=True)

avg_contrasts = filtered_gene_data.groupby("transition").agg({"contrast": "mean"}).reset_index()

# Create a DataFrame with the gene name
gene_names_df = pd.DataFrame({'gene_name': [gene_name] * len(avg_contrasts)})

# Add gene_name column to avg_contrasts
avg_contrasts["gene_name"] = gene_name

# Merge gene_names_df with avg_contrasts
merged_data = pd.merge(gene_names_df, avg_contrasts, on='gene_name', how='inner')

pivoted_data = avg_contrasts.pivot_table(index='gene_name', columns='transition', values='contrast').reset_index()
pivoted_data.columns.name = None  # To remove the column name 'transition' from the header
print(pivoted_data)

pivoted_data.to_csv("/Users/lfitzger/Test_BayesCode_desktop/average_contrasts_per_gene.csv", index=False)
