# clownfishBarEvolution
Updated 3.18.25 Lucy Fitzgerald

Scripts are used to generate the results and make the figures for the paper:
**Genetic basis of the evolution of vertical bars in clownfishes** 
(https://www.biorxiv.org/content/10.1101/2025.03.14.643259v1)

**ASR:**
- *stoachasticmapping_script4states.R*: to test different models of ancestral state reconstruction 
- *Figure1ASRplot.R*: to plot Figure 1 from the paper

**BayesCode:**
- *01_Transition_GainorLoss_Allgenes.py*: takes the simmaps generated from the ASR and the Bayescode results (omegas and contrasts) and assigns categories of gain/loss/same to the different nodes and calculates the weights of each transition for every gene
    - Names_List.txt: text file used in script to make sure species names match between Bayescode and ASR results 
- *02_weighted_regression_AllGenes.py*: summarizes the results for all the genes in a .tsv file
- *03_Correct_Pvals_BH.R*: correcting the p-values for all of the genes
- *04_AverageOmega.py*: calculates the average omega for all of the significant genes
- *Figure3and4Bayescodeplot.R*: plots Figures 3 and 4 from the paper for significant gain and loss color genes

**CodeML:**
- *01_Create_CTLFile_CodeML.py*: creates the codeml ctl files used to run codeml
  -  ClownTree_OnlyNeededSpecies_Rooted.Gain.tree: input tree file with foreground branches for gain genes
  -  ClownTree_OnlyNeededSpecies_Rooted.Loss.tree: input tree file with foreground branches for loss genes
- *02_ChecklogL_foregroundw.py*: summarizes the results for the three runs of codeml, computes log-likelihood and p-value for best run for each gene
- *03_Correct_Pvals_BH.R*: corrects p-values for all of the genes

*Categorizingallbranches(SM).R*: used to plot supplementary figure 2 
