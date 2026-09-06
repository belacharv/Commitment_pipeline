########################################################################################
### 06 CONFIDENCE INTERVALS IN SIMULATIONS WITH CHANGING PROBABILITY OF COMMITMENT ###
########################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

# This script works with the 2nd layer (sequences present in exactly two samples)

# It uses the same technique as in parameter recovery, but it adjusts the input values.
# In the first iteration, the parameter recovery method is used on real calculated values,
# then it enters to for loop where we adjust Pc for 0.005 in each step
# and then observe in which confidence intervals belongs our calculated value 
# The point is to specify the accuracy of the predicted calculated Pc and true_cd4

# Output is an datatable with all the CIs for each Pc 
# This table is later used in script 05.2

import numpy as np
#import matplotlib.pyplot as plt
import sympy as sp
import pandas as pd
from scipy import stats
from T_cell_lineage_model_class_nsolve import TcellLineageModel
from T_cell_lineage_model_class_nsolve import fit_model_to_data
import os
# Parameter recovery
if not os.path.exists("datatables/simulations_parameter_recovery_ci"):
    os.makedirs("datatables/simulations_parameter_recovery_ci")
os.chdir("datatables")
def simulate(pc, true_n4,model,total_seq,chain,cd4cd4):
    # Generate multiple datasets and see if we can recover parameters
    recovery_results = []
    n_simulations = 100

    for i in range(n_simulations):
        # Simulate data
        cd4_cd4, cd4_cd8, cd8_cd8 = model.simulate_dataset(seed=i)
        print(cd4_cd4, cd4_cd8, cd8_cd8)
        # Fit model to simulated data
        solutions = fit_model_to_data(cd4_cd4, cd8_cd8, total_seq,chain)
    
        if len(solutions) >= 1:
        # Take the first valid solution
            fitted_Pc, fitted_n4 = solutions[0]        
            recovery_results.append({
                'sim_id': i,
                'true_Pc': pc,
                'true_n4': true_n4,
                'fitted_Pc': fitted_Pc,
                'fitted_n4': fitted_n4,
                'cd4_cd4': cd4_cd4,
                'cd8_cd8': cd8_cd8,
                'error_Pc': abs(fitted_Pc - pc),
                'error_n4': abs(fitted_n4 - true_n4),
                'error_cd4_cd4': abs(cd4_cd4 - cd4cd4)
            })
        else:
            print("*"*50)
            print(f"No valid solutions found for simulation {i}")
            print("*"*50)
            print(cd4_cd4, cd4_cd8, cd8_cd8)
            print("*"*50)

    return recovery_results

def print_results(pc,true_n4,i,model,total_seq,chain,cd4cd4):
    recovery_results = simulate(pc,true_n4,model,total_seq,chain,cd4cd4)
    df_recovery = pd.DataFrame(recovery_results)
    n_simulations = 100
    print(df_recovery)
    print(f"Successful recoveries: {len(df_recovery)}/{n_simulations}")
    CIs = []

    if len(df_recovery) > 0:
        print(f"Mean Pc error: {df_recovery['error_Pc'].mean():.4f} ± {df_recovery['error_Pc'].std():.4f}")
        print(f"Mean n4 error: {df_recovery['error_n4'].mean():.4f} ± {df_recovery['error_n4'].std():.4f}")
        print(f"Mean cd4_cd4 error: {df_recovery['error_cd4_cd4'].mean():.4f} ± {df_recovery['error_cd4_cd4'].std():.4f}")
        name = str(i)+"_"+str(chain)+"_pc_"+str(pc)+"_simulations.csv"
        df_recovery.to_csv(name, index=False)

        fitted_Pc = df_recovery['fitted_Pc']

        confidence_level = 0.95
        mean = np.mean(fitted_Pc)
        sem = stats.sem(fitted_Pc)  # Standard error of the mean
        ci = stats.t.interval(confidence_level, len(fitted_Pc)-1, loc=mean, scale=sem)


        print(f"Mean: {mean}")
        print(f"95% Confidence Interval fitted_Pc: [{ci[0]:.8f}, {ci[1]:.8f}]")

        fitted_n4 = df_recovery['fitted_n4']
        mean4 = np.mean(fitted_n4)
        sem4 = stats.sem(fitted_n4)  # Standard error of the mean
        ci4 = stats.t.interval(confidence_level, len(fitted_n4)-1, loc=mean4, scale=sem4)

        print(f"Mean: {mean4}")
        print(f"95% Confidence Interval fitted_n4: [{ci4[0]:.8f}, {ci4[1]:.8f}]")

        cd4_cd4 = df_recovery['cd4_cd4']

        mean44 = np.mean(cd4_cd4)
        sem44 = stats.sem(cd4_cd4)  # Standard error of the mean
        ci44 = stats.t.interval(confidence_level, len(cd4_cd4)-1, loc=mean44, scale=sem44)

        print(f"Mean: {mean44}")
        print(f"95% Confidence Interval cd4_cd4: [{ci44[0]:.8f}, {ci44[1]:.8f}]")
        CIs.append({
            'true_pc':pc,
            'mean_pc':mean,
            'ci_pc_1':ci[0],
            'ci_pc_2':ci[1],
            'mean_n4':mean4,
            'ci_n4_1':ci4[0],
            'ci_n4_2':ci4[1],
            'mean_cd4_cd4':mean44,
            'ci_44_1':ci44[0],
            'ci_44_2':ci44[1],

        })

    else:
        print("No successful parameter recoveries found!")
    
    df_CIs = pd.DataFrame(CIs)
    return df_CIs

#print_results()
def main(chains):
    
# Test with known parameters
    #chains = ["a","b"]
    for chain in chains:
        os.chdir("simulations_parameter_recovery_ci")
        input_row = generate_input_row(chain)
        true_Pc = input_row[0]
        pc = input_row[3]
        true_n4 = input_row[1]
        true_n8 = input_row[2]
        total_seq = int(true_n4 + true_n8)
        cd4cd4 = input_row[4]
        ci = []
        df_ci = pd.DataFrame(ci)

        print(f"Simulations of Pc = {true_Pc}")
        print("="*60)
        model = TcellLineageModel(true_Pc, true_n4, total_seq, chain)
        df_new = print_results(true_Pc, true_n4,0,model,total_seq,chain,cd4cd4) 
        df_ci = pd.concat([df_ci,df_new],ignore_index=True)

    # we change pc every
        for i in range(1,11):
            print(f"Simulations of Pc = {pc}")
            print("="*60)
            model = TcellLineageModel(pc, true_n4, total_seq,chain)
            print(f"True parameters: Pc = {pc}, n4 = {true_n4}")
            df_ci_new_row = print_results(pc,true_n4,i,model,total_seq,chain,cd4cd4)
            df_ci = pd.concat([df_ci,df_ci_new_row],ignore_index=True)
            pc = pc - 0.001 # if we want to decrease, we need to adjust this
    
        pc = input_row[3]

        for i in range(12,22):
            print(f"Simulations of Pc = {pc}")
            print("="*60)
            model = TcellLineageModel(pc, true_n4, total_seq,chain)
            print(f"True parameters: Pc = {pc}, n4 = {true_n4}")
            df_ci_new_row = print_results(pc,true_n4,i,model,total_seq,chain,cd4cd4)
            df_ci = pd.concat([df_ci,df_ci_new_row],ignore_index=True)
            pc = pc + 0.001 # if we want to decrease, we need to adjust this
        os.chdir("..")
        name_final_table = "06_"+str(chain)+"_ci_simulations_table_new.csv"
        df_ci.to_csv(name_final_table, index=False)

def generate_input_row(chain):
    if chain == "a":
        true_Pc = 0.934144763175453 # calculated probability of commitment (how precice the cell fate desicion is)
        true_n4 = 305 # number of CD4 biased sequences
        true_n8 = 10 # number of CD8 biased sequences
        total_seq = 315 # all sequences in 2nd layer
        chain = "a"
        pc = 0.934
        cd4cd4 = 258
    elif chain =="b":
        # for beta chain experiment
        true_Pc = 0.8038138880931018 # calculated probability of commitment (how precice the cell fate desicion is)
        true_n4 = 94.9 # number of CD4 biased sequences
        true_n8 = 3.1 # number of CD8 biased sequences
        total_seq = 98 # all sequences in 2nd layer
        pc = 0.803
        cd4cd4 = 94
    elif chain == "Vb5":
        true_Pc = 0.7225935629263628
        true_n4 = 34766.7
        true_n8 = 4754.3
        total_seq = 39521
        pc = 0.723
        cd4cd4 = 16832
    elif chain == "DeGreef":
        true_Pc = 0.9780546333214499
        true_n4 = 39630.3
        true_n8 = 22490.7
        total_seq = 62121
        pc = 0.978
        cd4cd4 = 37124
    params = [true_Pc,true_n4,true_n8,pc,cd4cd4]
    return params

#main(["a","b"])