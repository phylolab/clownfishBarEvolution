#!/usr/bin/env python3
from collections import defaultdict
import os
import numpy as np
from ete3 import Tree
import pandas as pd
from scipy.stats import mannwhitneyu
import scipy.stats as stats


# Specify the paths for input and output files
workdir = './Users/lfitzger/Test_BayesCode_desktop/'
input_directory = '/Users/lfitzger/Test_BayesCode_desktop/50k_cluster/1_results_NoTraits50k'
output_directory = '/Users/lfitzger/Test_BayesCode_desktop/Results_DiscreteTraits/model0.6_maps4states/simapvsomega_allgenes' #/Results/simapvsomega_allgenes.zip

# Define the directory where your files are located
path_simap = '/Users/lfitzger/Test_BayesCode_desktop/Models/model0.6_clownfish4states.map' #same file as generated from ASR "Results"

# Read gene names from the file
with open('/Users/lfitzger/Test_BayesCode_desktop/AllGenes.txt', 'r') as f:
    gene_names = [name.strip() for name in f.read().splitlines()]

# Loop through each folder and process it
for gene_name in gene_names:
    # Construct paths
    path_omega = os.path.join(input_directory,gene_name,f'{gene_name}_nodeomega.Omega.nhx')
    path_output = os.path.join(output_directory,f'{gene_name}.simap_vs_omega.tsv')

    print(f"Processing gene: {gene_name}")

    # Read the list of names and their corresponding species
    d = {}
    with open("/Users/lfitzger/Test_BayesCode_desktop/Names_List.txt") as f:
        for line in f:
            (key, val) = line.split()
            d[key] = val.replace('._', '')

    # Read the tree with omega values
    tree_omega = Tree(path_omega, format=1)
    tree_omega.name = 'root'
    omega_dico = {n.name: n for n in tree_omega.traverse()}
    omega_set_leaves = {n.name: set(n.get_leaf_names()) for n in tree_omega.traverse()}
    omega_leaves = set(tree_omega.get_leaf_names())


    # Get the omega value for each node in the simap tree
    def get_omega_node(n_simap):
        if n_simap.name in omega_leaves:
            return omega_dico[n_simap.name]

        else:
            simap_leaves = set(n_simap.get_leaf_names())
            # Find in node_set_leaves the node with the same set of leaves
            for node_name, node_set in omega_set_leaves.items():
                if simap_leaves == node_set:
                    return omega_dico[node_name]


    # Parse a simap tree
    def parse_simap(str_simap):
        str_simap = str_simap.replace('{', '[&&NHX:simap=').replace('}', ']').replace('._', "")
        for state in ['One', 'Two', 'Three', 'None']:
            str_simap = str_simap.replace(f"{state},", f"{state}_").replace(f":{state}", f"__{state}")
        for key, value in d.items():
            node = tree_omega.search_nodes(name=key)[0]
            str_simap = str_simap.replace(f"{value}:", f"{key}:{node.dist}")
        return str_simap.replace(':[', '[')


    # Read the simap trees
    list_simap = []
    with open(path_simap, 'r') as file:
        for str_newick in file:
            str_newick = parse_simap(str_newick)
            tree_simap = Tree(str_newick, format=1)
            list_simap.append(tree_simap)

    # conversion from states to numbers
    state_to_num = {
        "None": 0,
        "One": 1,
        "One+Two": 2,
        "Two": 2,
        "Three": 3,
        "Two+Three": 3,
    }

    # Get the transition for each node in the simap tree
    simap_dict = defaultdict(list)
    for i, tree_simap in enumerate(list_simap):
        for node_simap in tree_simap.traverse():
            omega_node = get_omega_node(node_simap)
            node_simap.name = omega_node.name
            if node_simap.name == 'root':
                continue
            simap_split = node_simap.simap.split('__')
            assert len(simap_split) >= 1, simap_split
            # The key are the states and the values are the transition time
            transition_simap = {i.split('_')[0]: i.split('_')[1] for i in simap_split}
            states = list(transition_simap.keys())
            if len(states) == 1:
                simap_dict[node_simap.name].append("same")
            else:
                if state_to_num[states[-1]] < state_to_num[states[0]]:
                    simap_dict[node_simap.name].append("gain")
                elif state_to_num[states[-1]] > state_to_num[states[0]]:
                    simap_dict[node_simap.name].append("loss")
                else:
                    simap_dict[node_simap.name].append("same")
    # Write the output file as a dataframe
    output_dict = defaultdict(list)
    for node_name, list_transitions in simap_dict.items():

        omega = float(omega_dico[node_name].Omega)
        omega_parent = float(omega_dico[node_name].up.Omega)
        time = omega_dico[node_name].dist
        contrast = (omega - omega_parent) / np.sqrt(time)

    # Count occurrences of each transition
        transition_counts = defaultdict(int)
        for transition in list_transitions:
            transition_counts[transition] += 1
        
        # Calculate weights for each transition
        total_transitions = sum(transition_counts.values())
        for transition, count in transition_counts.items():
            output_dict['node_name'].append(node_name)
            output_dict['omega'].append(omega)
            output_dict['contrast'].append(contrast)
            output_dict['transition'].append(transition)
            output_dict['weights'].append(count / total_transitions)

    output_df = pd.DataFrame(output_dict)
    output_df.to_csv(f"{path_output}", sep="\t", index=False)
