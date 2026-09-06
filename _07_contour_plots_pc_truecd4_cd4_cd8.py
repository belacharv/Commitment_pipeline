################################################################################
### 07 COUNTOUR PLOTS/HEATMAPS ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

# Current script focuses on the 2nd layer (sequences present in exactly two samples)
# The goal is to draw a function estimating the correlation between Pc and true_cd4


# importing libraries
import numpy as np
import matplotlib.pyplot as plt   
from pathlib import Path
import os

script_dir = Path(__file__).parent
print(script_dir)
datatables_dir = script_dir.parent / "datatables"
plots_dir = script_dir.parent / "plots"
for directory in ["datatables","plots"]:
    #if os.path.exists(directory)==False:
    current_dir = script_dir / directory
    print(current_dir)
    if not os.path.exists(current_dir):
        os.mkdir(directory)



# function to return the three main equations:
    # cd4-cd4, cd4-cd8 and cd8-cd8
# input:
    # Pc_val = probability of commitment
    # n4_val = true_cd4 value (number of CD4 biased sequences)
    # all_seq = all the sequences in the layer
    # chain = alpha or beta
def calculate_exp(Pc_val, n4_val, all_seq, chain):
    #Calculate expected CD4-CD4 count for given Pc and n4
    # coefficients based on the dataset that is analyzed
        # these coefficients depend on the availability of CD4 or CD8 samples
    if chain == "a" or chain == "Vb5":
        # in alpha experiment of Giorgetti et al or Vb5 dataset, they had 5 CD4 and 5 CD8 samples
        # we can statistically think about it as if we were chosing two samples out of the group
        c44 = 5*4 / (10*9) # first CD4 chosen among 5 CD4 samples from 10 samples in total, second CD4 chosen among the 4 CD4 samples left
            # (we had to exclude the first one cause we already chosen it) from the 9 samples left
            # hence first 5/10 and second 4/9 
        c48 = 5*5 / (10*9)
        c88 = 5*4 / (10*9)
    # different number of samples for beta chain
    elif chain == "b":
        # in beta experiment of Giorgetti et al dataset, they had 13 CD4 samples and 10 CD8 samples
        c44 = 13*12 / (23*22)
        c48 = 13*10 / (23*22)
        c88 = 10*9 / (23*22)

    elif chain == "DeGreef":
        # in deGreef et al dataset they had 3 CD4 samples and CD8 samples
        c44 = 3*2/(6*5)
        c48 = 3*3/(6*5)
        c88 = 3*2/(6*5)


    n8_val = all_seq - n4_val
    
    # For CD4 biased sequences
    p4_44_val = Pc_val**2 * c44 # both of them are cd4 samples, hence both chose correctly
    p4_48_val = 2*Pc_val*(1-Pc_val) * c48 # one cd4, one cd8, hence one chose correctly, the other not
    p4_88_val = (1-Pc_val)**2 * c88 # two cd8s, both of them chose wrong
    p4_total_val = p4_44_val + p4_48_val + p4_88_val
    
    # For CD8 biased sequences
    p8_44_val = (1-Pc_val)**2 * c44
    p8_48_val = 2*Pc_val*(1-Pc_val) * c48
    p8_88_val = Pc_val**2 * c88
    p8_total_val = p8_44_val + p8_48_val + p8_88_val
    
    # three equations with cd4-cd4, cd4-cd8 and cd8-cd8 on the left side
    # on the right site: 
        # example: n4_val * (p4_44_val/p4_total_val) + n8_val * (p8_44_val/p8_total_val)
            # first term n4_val * (p4_44_val/p4_total_val)
                # (number of cd4 biased sequences) * (probability that both will chose cd4-cd4 if the seq was cd4 biased ~ both are correct) 
            # second term n8_val * (p8_44_val/p8_total_val)
                # (number of cd8 biased sequences) * (probability that it will chose cd8-cd8 even though the seq was cd4 biased ~ both are wrong) 
    exp_44 = n4_val * (p4_44_val/p4_total_val) + n8_val * (p8_44_val/p8_total_val)
    exp_48 = n4_val * (p4_48_val/(p4_44_val+p4_48_val+p4_88_val)) + n8_val * (p8_48_val/(p8_44_val+p8_48_val+p8_88_val))
    exp_88 = n4_val * (p4_88_val/(p4_44_val+p4_48_val+p4_88_val)) + n8_val * (p8_88_val/(p8_44_val+p8_48_val+p8_88_val))
    
    return exp_44,exp_48,exp_88


def create_heatmap(all_seq, chain, observed, parameter):

    # Create grid - use continuous values for both axes
    n4_values = np.linspace(0, all_seq, 1000)  # Continuous n4 values
    pc_values = np.linspace(0.5, 1.0, 1000)  # Continuous Pc values starting from 0.5 (we expect the probability of commitment to be higher than 50%)
        # cause commitment 50% means that it is random, cause there are two choices - either cd4 or cd8
    # initial set-up
    cd4_cd4_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    cd4_cd8_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    cd8_cd8_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    
    # we need three grids, each one for each category corresponding to cd4-cd4, cd4-cd8 or cd8-cd8
    # in these grids we'll store the the info how far the solution is from the observed one
    diff_grids = [cd4_cd4_diff_grid,cd4_cd8_diff_grid,cd8_cd8_diff_grid]

    for i, n4 in enumerate(n4_values):
        for j, pc in enumerate(pc_values):
            # calculate with given pc and n4=true_cd4 
            val44,val48,val88 = calculate_exp(pc, n4, all_seq, chain)
            # calculate how much the results differ from observed values
            cd4_cd4_diff_grid[i, j] = abs(observed[0] - val44)
            cd4_cd8_diff_grid[i, j] = abs(observed[1] - val48)
            cd8_cd8_diff_grid[i, j] = abs(observed[2] - val88)
    
    # we do the countour plots for true_cd4 and true_cd8
    if parameter == "true_cd4":
        y_axis = n4_values
        org = "lower"
        ylabel = "True CD4"
    else: 
        if parameter == "true_cd8":
            y_axis = all_seq - n4_values # we change the y axis
            org = "upper"
            ylabel = "True CD8"
        else:
            raise Exception("Parameter not valid")

    # Create heatmap
    # Use colormap where 0 (perfect match) is red and large differences are dark blue/green
    for i in range(0,len(diff_grids)):
        contour = plt.contour(pc_values, y_axis, diff_grids[i], 
                         levels=15, colors='black', linewidths=3)
        plt.clabel(contour, inline=True, fontsize=28)#, #fmt=f'{observed[i]}')
        contour = plt.contour(pc_values, y_axis, diff_grids[i], 
                         levels=[0.05], colors='black', linewidths=3)
        plt.clabel(contour, inline=True, fontsize=14)
        im = plt.imshow(diff_grids[i], 
                        aspect='auto',
                        origin=org,
                        extent=[pc_values.min(), pc_values.max(), y_axis.min(), y_axis.max()],
                    cmap='coolwarm_r', 
                    )
        print(diff_grids[i].min(), diff_grids[i].max())
        cbar = plt.colorbar(im, label='Difference |Observed - Expected|')
        cbar.ax.tick_params(labelsize=28)
        plt.xlabel('Pc', fontsize=28)
        plt.ylabel(ylabel, fontsize=28)
        plt.xticks(fontsize=28)
        plt.yticks(fontsize=28)
        plt.tight_layout()
        # name and save the figure
        name = "07_"+str(chain)+"_0"+str(i+1)+"_"+parameter+"_heatmap.png"
        plt.savefig(name)
        plt.close()