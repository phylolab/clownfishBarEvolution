import os
import sys

# Changes to Gain of Vertical bar
#path_alignments = "/scratch/lfitzger/CodeML_analyses/alidata/" # PATH where alignment files (PHYLIP format) are found
#path_to_ctl_files = "/scratch/lfitzger/CodeML_analyses/Gain_verticalbar/ctl_files/" # PATH were output control files are saved
#path_to_input_align_for_ctl = "/scratch/lfitzger/CodeML_analyses/alidata/"  # PATH for codeML output file (PATH to be used on control file)
#tree_file_name = "/scratch/lfitzger/CodeML_analyses/Gain_verticalbar/Clown.DatedPhylogeny.CodeML_Gain.tree" # PATH to labelled tree (PATH to be used on control file)
#path_to_output_for_ctl = "/scratch/lfitzger/CodeML_analyses/Gain_verticalbar/results/" # PATH to CodeML output file (PATH to be used on control file)
#out_suff = ".Gain.Vert.Bar."                  # Control file suffix

# Changes to Loss of Vertical bar
path_alignments = "/scratch/lfitzger/CodeML_analyses/alidata/" # PATH where alignment files (PHYLIP format) are found
path_to_ctl_files = "/scratch/lfitzger/CodeML_analyses/Loss_verticalbar/ctl_files" # PATH were output control files are saved
path_to_input_align_for_ctl = "/scratch/lfitzger/CodeML_analyses/alidata/"  # PATH for codeML output file (PATH to be used on control file)
tree_file_name = "/scratch/lfitzger/CodeML_analyses/Loss_verticalbar/Clown.DatedPhylogeny.CodeML_Loss.tree" # PATH to labelled tree (PATH to be used on control file)
path_to_output_for_ctl = "/scratch/lfitzger/CodeML_analyses/Loss_verticalbar/results/" # PATH to CodeML output file (PATH to be used on control file)
out_suff = ".Loss.Vert.Bar."                  # Control file suffix

Dictionary_Specification_models = {"H0":{"model":"2", "NSsites":"2", "fix_omega":"1", "omega":"1"}, "H1":{"model":"2", "NSsites":"2", "fix_omega":"0", "omega":"1"}, "M1a":{"model":"0", "NSsites":"1", "fix_omega":"0", "omega":"1"}}


with open('/scratch/lfitzger/CodeML_analyses/AllGenes.txt', 'r') as file:
    gene_names = [line.strip() for line in file.readlines()]

list_of_files = [f"{gene}_OnlyNeededSpecies.fa.phy" for gene in gene_names if os.path.isfile(os.path.join(path_alignments, f"{gene}_OnlyNeededSpecies.fa.phy"))]

for each_file in list_of_files:

    file_name = each_file.split("_")
    outfilemame = file_name[0]+"_"+file_name[1]
    
    align_name = path_to_input_align_for_ctl + each_file
    
    for each_model in Dictionary_Specification_models.keys():
        out_name_results = path_to_output_for_ctl + outfilemame + out_suff + each_model + ".out"
                        
        ctlfile = """
     seqfile = %s * sequence data file name
    treefile = %s * tree structure file name
     outfile = %s * main result file name

       noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
     verbose = 1   * 1: detailed output, 0: concise output
     runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

     seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
   CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       ndata = 1   * specifies the number of separate data sets in the file
       clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree

      aaDist = 0   * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
                   * 7:AAClasses

       model = %s   * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                        * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)

     NSsites = %s   * 0:one w; 1:neutral; 2:positive selection; 3:discrete; 4:freqs;
                   * 5:gamma; 6:2gamma; 7:beta; 8:beta&w; 9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

       icode = 0   * 0:universal code; 1:mammalian mt; 2-11:see below
       Mgene = 0   * 0:rates, 1:separate;

   fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
       kappa = 2.05154 * initial or fixed kappa

   fix_omega = %s   * 1: omega or omega_1 fixed, 0: estimate
       omega = %s

       getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .5e-6 * small value used in the difference approximation of derivatives

   cleandata = 1   * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed
      method = 0   * 0: simultaneous; 1: one branch at a time
     """ % (align_name, tree_file_name, out_name_results, 
            Dictionary_Specification_models[each_model]["model"], 
            Dictionary_Specification_models[each_model]["NSsites"], 
            Dictionary_Specification_models[each_model]["fix_omega"], 
            Dictionary_Specification_models[each_model]["omega"])
        
        ctl_file_path = os.path.join(path_to_ctl_files, outfilemame + out_suff + each_model + ".ctl")
        with open(ctl_file_path, "w") as outf:
            outf.write(ctlfile)
