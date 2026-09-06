################################################################################
### 05 PARAMETER RECOVERY ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

## The final sequences that are not unique (filtered in first step and while creating count table)
# have been collected and seperated into layers based on in how many samples they are present.
# In experiment alpha, where they edited alpha chain, they had 15 samples from spleen, 10 CD4 and 5 CD8.
# Each sequence have been assigned to one of the 14 layers (ranging from 2 to 15) based on in how many samples it is present
# In current python script, we focus only on the 2nd layer, hence sequences present in exactly two samples

### This script validate the solution provided by depth_2_calculation using parameter recovery technique and estimates confidence intervals
# It is based on parameter recovery method, where we input our results, we calculate the probabilities, generate new dataset and then calculate Pc and true_cd4
# again. This serves as a control that we get the same result once we substitute to the equation. Moreover, we can get calculated ci

# importing libraries
import numpy as np
import pandas as pd
from scipy import stats
#from pathlib import Path
# import T cell lineage model class from different script
from T_cell_lineage_model_class_nsolve import TcellLineageModel
from T_cell_lineage_model_class_nsolve import fit_model_to_data

# script_dir = Path(__file__).parent
# print(script_dir)
# datatables_dir = script_dir.parent / "datatables"
# plots_dir = script_dir.parent / "plots"
import os

script_dir_os = os.path.dirname(os.path.abspath(__file__))
datatables_dir_os = os.path.join(script_dir_os, "datatables")
plots_dir_os = os.path.join(script_dir_os, "plots")
for directory in [datatables_dir_os,plots_dir_os]:
    if not os.path.exists(directory):
        os.makedirs(directory)
# Test with known parameters

# true_Pc = 0.915186
# true_n4 = 372
# true_n8 = 14
# total_seq = 386
# chain = "a"
# model = TcellLineageModel(true_Pc, true_n4, total_seq,chain)

# print(f"True parameters: Pc = {true_Pc}, n4 = {true_n4}")
recovery_statistics = []

# function to create 100 simulated dataset based on probabilities calculated by substitution of our obtained results
def parameter_recovery(chain,recovery_statistics):
    # results from 2nd layer divided by chain
    # for alpha chain experiment
    if chain == "a":
        true_Pc = 0.934144763175453 # calculated probability of commitment (how precice the cell fate desicion is)
        true_n4 = 304.7 # number of CD4 biased sequences
        true_n8 = 10.3 # number of CD8 biased sequences
        total_seq = 315 # all sequences in 2nd layer
        chain = "a"
        model = TcellLineageModel(true_Pc, true_n4, total_seq,chain)
    elif chain =="b":
        # for beta chain experiment
        true_Pc = 0.8038138880931018 # calculated probability of commitment (how precice the cell fate desicion is)
        true_n4 = 94.9 # number of CD4 biased sequences
        true_n8 = 3.1 # number of CD8 biased sequences
        total_seq = 98 # all sequences in 2nd layer
        model = TcellLineageModel(true_Pc, true_n4, total_seq,chain) 
    elif chain == "Vb5":
        true_Pc = 0.7225935629263628
        true_n4 = 34766.7
        true_n8 = 4754.3
        total_seq = 39521
        model = TcellLineageModel(true_Pc,true_n4,total_seq,chain)
    else:
        true_Pc = 0.9780546333214499
        true_n4 = 39630.3
        true_n8 = 22490.7
        total_seq = 62121
        model = TcellLineageModel(true_Pc, true_n4, total_seq,chain)

    print(f"True parameters: Pc = {true_Pc}, n4 = {true_n4}")
    recovery_results = []
    n_simulations = 100 # number of simulations

    for i in range(n_simulations):
        # simulating the dataset based on calculated probabilities
        # returns 
            # cd4_cd4, number of sequences present in two CD4 samples
            # cd4_cd8, number of sequences present in one CD4 and one CD8 sample
            # cd8_cd8, number of sequences present in two CD8 samples
        cd4_cd4, cd4_cd8, cd8_cd8 = model.simulate_dataset(seed=i)
        print(cd4_cd4, cd4_cd8, cd8_cd8)
        # fit model to the simulated data
        first_solution = fit_model_to_data(cd4_cd4, cd8_cd8,total_seq,chain)
        if len(first_solution) != 0: # if we received any solution
            fitted_Pc, fitted_n4 = first_solution[0]
        if fitted_Pc is not None:
            # one row in our result table
            recovery_results.append({
            'sim_id': i,
            # original commitment and true_cd4 parameters calculated before
            'true_Pc': true_Pc, 
            'true_n4': true_n4, 
            # fitted Pc and true_cd4 on the new simulated dataset
            'fitted_Pc': fitted_Pc, 
            'fitted_n4': fitted_n4,
            # numbers of each category among the 2nd layer
            'cd4_cd4': cd4_cd4, # cd4 samples only
            'cd4_cd8':cd4_cd8, # one cd4 and one cd8 sample
            'cd8_cd8': cd8_cd8, # cd8 samples only
            # the difference between the true parameter and the fitted parameter
            'error_Pc': abs(fitted_Pc - true_Pc),
            'error_n4': abs(fitted_n4 - true_n4),
            'error_n8': abs(total_seq-fitted_n4- true_n8)
            })
    # creating a table from 100 rows obtained from simulation results
    df_recovery = pd.DataFrame(recovery_results)
    print(df_recovery)
    print(f"Successful recoveries: {len(df_recovery)}/{n_simulations}")
    print(f"Mean Pc error: {df_recovery['error_Pc'].mean():.4f} ± {df_recovery['error_Pc'].std():.4f}")
    print(f"Mean true_cd4 error: {df_recovery['error_n4'].mean():.4f} ± {df_recovery['error_n4'].std():.4f}")
    print(f"Mean true_cd8 error: {df_recovery['error_n8'].mean():.4f} ± {df_recovery['error_n8'].std():.4f}")
    name = "datatables/05_parameter_recovery_"+chain+"_pc_estimation.csv"
    df_recovery.to_csv(name, index=False)

    # statistics for Pc
    fitted_Pc = df_recovery['fitted_Pc']
    confidence_level = 0.95
    mean = np.mean(fitted_Pc)
    sem = stats.sem(fitted_Pc) 
    # confidence intervals of Pc
    ci = stats.t.interval(confidence_level, len(fitted_Pc)-1, loc=mean, scale=sem)
    print(f"Mean: {mean}")
    print(f"95% Confidence Interval fitted_Pc: [{ci[0]:.8f}, {ci[1]:.8f}]")

    # statistics for true_cd4
    fitted_n4 = df_recovery['fitted_n4']
    mean4 = np.mean(fitted_n4)
    sem4 = stats.sem(fitted_n4)  
    # confidence intervals of true_cd4
    ci4 = stats.t.interval(confidence_level, len(fitted_n4)-1, loc=mean4, scale=sem4)
    print(f"Mean: {mean4}")
    print(f"95% Confidence Interval fitted_n4: [{ci4[0]:.8f}, {ci4[1]:.8f}]")

    # one row containing statistics results for this chain
    recovery_statistics.append({
    'chain': chain, 
    'success_recoveries_from_100': len(df_recovery),
    # Results for commitment
    'mean_Pc':mean, 
    'mean_Pc_error': df_recovery['error_Pc'].mean(),
    'Pc_std': df_recovery['error_Pc'].std(),
    'Pc_ci_1':ci[0], 
    'Pc_ci_2':ci[1],
    # Results for true_cd4
    'mean_true_cd4':mean4,
    'mean_error_true_cd4': df_recovery['error_n4'].mean(),
    'true_cd4_std': df_recovery['error_n4'].std(),
    'true_cd4_ci_1':ci4[0], 
    'true_cd4_ci_2':ci4[1]
    })

    print(f"Mean Pc error: {df_recovery['error_Pc'].mean():.4f} ± {df_recovery['error_Pc'].std():.4f}")
    print(f"Mean n4 error: {df_recovery['error_n4'].mean():.4f} ± {df_recovery['error_n4'].std():.4f}")
    print(f"Mean n8 error: {df_recovery['error_n8'].mean():.4f} ± {df_recovery['error_n8'].std():.4f}")

# call the function to simulate the datasets and calculate fitted Pc and fitted true_cd4
#parameter_recovery("a",recovery_statistics) # for alpha chain
#parameter_recovery("b",recovery_statistics) # for beta chain
#parameter_recovery("Vb5",recovery_statistics)
#df_recovery_statistics = pd.DataFrame(recovery_statistics)
#df_recovery_statistics.to_csv("datatables/05_parameter_recovery_pc_statistics.csv", index=False)
