import os
import pandas as pd

# Define the gene name
gene_name = "chr06_g3597.t1"

# Read the TSV file
file_path = os.path.join("/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/simapvsomega_allgenes", f"{gene_name}.simap_vs_omega.tsv")
gene_data = pd.read_csv(file_path, sep="\t")

# Keep rows with the highest weights for each node_name
gene_data = gene_data.loc[gene_data.groupby("node_name")["weights"].idxmax()]

# Display the filtered DataFrame
print(gene_data)

# Save the filtered data
gene_data.to_csv("filtered_output.csv", index=False, sep="\t")
