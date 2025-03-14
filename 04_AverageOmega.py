import pandas as pd
import os

# Define file paths
gain_gene_list_file = "/Users/lfitzger/Test_BayesCode_desktop/Sig_bayescode_gaincolorgenes.txt"
loss_gene_list_file = "/Users/lfitzger/Test_BayesCode_desktop/Sig_bayescode_losscolorgenes.txt"
input_folder = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/simapvsomega_allgenes"
output_file = "/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/Final_model0.6_maps4states/Avg_Omega_Significant_ColorGenes.csv"

# Define branches for gain and loss
gain_branches = {"A._biaculeatus", "A._latezonatus", "A._tricinctus", "node_12", "node_24"}
loss_branches = {"A._ephippium", "A._mccullochi", "A._nigripes", "node_13", "node_8"}

# Read gain significant genes
with open(gain_gene_list_file, "r") as file:
    significant_gain_genes = {gene.strip() for gene in file.read().splitlines()}  

# Read loss significant genes
with open(loss_gene_list_file, "r") as file:
    significant_loss_genes = {gene.strip() for gene in file.read().splitlines()}  

# Initialize a list to store results
avg_omega_list = []

# Loop through all TSV files in the input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".tsv"):  
        gene_name = file_name.replace(".simap_vs_omega.tsv", "")  
        gene_file = os.path.join(input_folder, file_name)
        
        # Read gene data
        gene_data = pd.read_csv(gene_file, sep="\t")

        # Group by 'node_name' and find the maximum weight within each group
        max_weights = gene_data.groupby('node_name')['weights'].transform('max')

        # Filter to keep only rows with the highest weights
        filtered_gene_data = gene_data[gene_data['weights'] == max_weights].copy()

        # Compute overall average omega
        avg_omega = filtered_gene_data['omega'].mean()

        # Compute gain omega (if the gene is in gain-related categories)
        gain_data = filtered_gene_data[filtered_gene_data['node_name'].isin(gain_branches)]
        avg_omega_gain = gain_data['omega'].mean() if not gain_data.empty else None

        # Compute loss omega (if the gene is in loss-related categories)
        loss_data = filtered_gene_data[filtered_gene_data['node_name'].isin(loss_branches)]
        avg_omega_loss = loss_data['omega'].mean() if not loss_data.empty else None

        # Determine significance category
        if gene_name in significant_gain_genes and gene_name in significant_loss_genes:
            significance = "Significant_Gain_Loss"
        elif gene_name in significant_gain_genes:
            significance = "Significant_Gain"
        elif gene_name in significant_loss_genes:
            significance = "Significant_Loss"
        else:
            significance = "Insignificant"

        # **Only store significant genes**
        if significance != "Insignificant":
            avg_omega_list.append({
                "Gene": gene_name,
                "Avg_Omega": avg_omega,
                "Average_Omega_Gain": avg_omega_gain if significance in {"Significant_Gain", "Significant_Gain_Loss"} else None,
                "Average_Omega_Loss": avg_omega_loss if significance in {"Significant_Loss", "Significant_Gain_Loss"} else None,
                "Significance": significance
            })

# Convert the results into a DataFrame
avg_omega_df = pd.DataFrame(avg_omega_list)

# Export to CSV (only significant genes)
avg_omega_df.to_csv(output_file, index=False)

print(f"Processing complete. Average omega values for significant genes saved to: {output_file}")
