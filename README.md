# clownfishBarEvolution
Updated 3.18.25 Lucy Fitzgerald

Scripts are used to generate the results and make the figures for the paper **Genetic basis of the evolution of vertical bars in clownfishes**(https://www.biorxiv.org/content/10.1101/2025.03.14.643259v1)

**ASR Folder:**
- *stoachasticmapping_script4states.R*: this script was used to test different models of ancestral state reconstruction with the "Data" folder containing the data input for the script and "Results" folder containing the simmap results
- *Figure1ASRplot.R*: contains a script to plot Figure 1, a summary of the stochastic maps which was updated from the stoachasticmapping_script4states.R

**BayesCode Folder:**
- *01_Transition_GainorLoss_Allgenes.py*: takes the simmaps generated from the ASR and the Bayescode results (omegas and contrasts) and assigns categories of gain/loss/same to the different nodes and calculates the weights of each transition for every gene, uses the "Data", 1_results_NoTraits50k.zip and generates "Results", simapvsomega_allgenes.zip
    - Names_List.txt: text file used in script to make sure species names match between Bayescode and ASR results 
- *02_weighted_regression_AllGenes.py*: summarizes the results for all the genes in a .tsv file, "Results", All_genes_weighted_regression_results_slope50k_omega.tsv
- *03_Correct_Pvals_BH.R*: Correcting the p-values for all of the genes, "Results" All_genes_adjusted_pvalues_gain.csv and All_genes_adjusted_pvalues_loss.csv
- *04_AverageOmega.py*: calculates the average omega for all of the significant genes, "Results" Avg_Omega_Significant_Genes.csv
    - Sig_bayescode_gaingenes.txt: text file, input for script, list of significant BayesCode gain genes
    - Sig_bayescode_lossgenes.txt: text file, input for script, list of significant BayesCode loss genes 
- *Figure3and4Bayescodeplot.R*: Plots figures 3 and 4 from the paper for significant gain and loss color genes

**CodeML Folder:**
- *01_Create_CTLFile_CodeML.py*: Creates the codeml ctl files used to run codeml, raw alignment data is published in Gaboriau et al., 2024 (https://www.biorxiv.org/content/10.1101/2024.07.08.602550v1)
  -  ClownTree_OnlyNeededSpecies_Rooted.Gain.tree: input tree file with foreground branches for gain genes
  -  ClownTree_OnlyNeededSpecies_Rooted.Loss.tree: input tree file with foreground branches for loss genes
  -  "Results": results from running codeml three times for gain and loss analysis; All_gain_genes_results1.csv, All_gain_genes_results2.csv, All_gain_genes_results3.csv,
      All_loss_genes_results1.csv, All_loss_genes_results2.csv, All_loss_genes_results3.csv
- *02_ChecklogL_foregroundw.py*: summarizes the results for the three runs of codeml, computes log-likelihood and p-value for best run for each gene
  - "Results": LogL_Comparison_GainResults_with_LRT_pvalues_and_foreground_w_H1.csv and LogL_Comparison_LossResults_with_LRT_pvalues_and_foreground_w_H1.csv
- *03_Correct_Pvals_BH.R*: Corrects p-values for all of the genes
  - "Results": All_gain_genes_foregroundomega.csv and All_loss_genes_foregroundomega.csv
