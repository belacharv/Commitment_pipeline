################################################################################
### 08 CALCULATE PC AND TRUE_CD4 WITH 2 EQ IN DIFFERENT LAYERS ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

# Current script focuses on the layer from 2 to max(n_cd4,n_cd8)
# For the two variables (Pc and true_cd4), we always use two equations
# instead of nomenclature used for 2nd layer (cd4-cd4,cd4-cd8,cd8-cd8),
    # we use nomenclature based on number of cd8 samples (since we have fewer of them)
    # so for 2nd the categories will be: cd8_0, cd8_1 and cd8_2
    # where cd8_0 means that there are no cd8 samples, cd4-cd4
    # for cd8_1, there's one cd8 sample - cd4-cd8
    # for cd8_2, there are two cd8 samples, so it's cd8-cd8


import sympy as sp
import pandas as pd
from pathlib import Path
import math
import os
###################################################################################################################
#### 0. Input
# input_row = [maxdepth, s4, s8, all_seq] + cd8_values

# input function that for each depth return information about counts specific for that depth(=layer)
    # input_row = [maxdepth, s4, s8, all_seq] + cd8_values
# input:
    # depth - the current clonal abundance, analogous to the layer in graph
    # chain - which dataset is analyzed
    # dataTableDir - the directory where the datatable with results in each depth from 03 step is stored
def input_function(depth,chain,dataTableDir):  
    
    # read the distribution table from 03 R step of each group in each layer
    if chain == "alpha":
        s4 = 5 # cd4 samples
        s8 = 5 # cd8 samples
        df = pd.read_csv(dataTableDir/'03_a_seq_counts_distribution_layers.csv')
    elif chain == "beta":
        s4 = 13
        s8 = 10
        df = pd.read_csv(dataTableDir/'03_b_seq_counts_distribution_groups_layers.csv')
    elif chain == "Vb5":
        s4 = 5
        s8 = 5
        df = pd.read_csv(dataTableDir/"03_d_seq_counts_distribution_layers.csv")
    elif chain == "DeGreef":
        s4 = 3
        s8 = 3
        df = pd.read_csv(dataTableDir/"03_c_seq_counts_distribution_layers.csv")
    else:
        raise Exception("wrong input")
    # we start from 2nd layer (seq present in 2 samples), we need to adjust for beginning of the list
    i = depth - 2
    maxdepth = df.iloc[i, 0]
    all_seq = df.iloc[i, 1]
    cd8_values = df.iloc[i, 2:].tolist()

    # vital row storing all the info specific for each layer
        # maxdepth refers to current layer
        # all_seq is Sequence count in that particular layer
        # cd8_values are group divided on the amount of cd8s in that group
            # cd8_0 (cd4_only), cd8_1, cd8_2 ... cd8_max (all the samples of cd8 present)
                # cd8_max means that in lower layers than s8 it's cd8_only
                # in higher layers it's all the samples of cd8 present
    input_row = [maxdepth, s4, s8, all_seq] + cd8_values
    print(input_row)
    return input_row


###################################################################################################################
#### 1. Generating coefficients ####
# each probability is influenced by the s4 and s8 constants

## numerator ##
# coef_generator takes care of the balance in current number of cd4 and cd8 that were used and which were not
# input:
    # currentS4 - current number of CD4 samples
    # currentS8 - current number of CD8 samples
    # sample - string cd4 or cd8
def coeff_generator(currentS4, currentS8, sample):
    if sample == "cd4":
        return currentS4 - 1, currentS8, currentS4 # we note that we used one cd4 sample
    else: # the sample is cd8
        return currentS4, currentS8 - 1, currentS8 # we note that we used one cd8 sample
    
# needs empty coefficient [] before starts
# function to generate all the coefficients
# based on recursive aproach depth-first search to search graphs (dfs)
# generates the numenator combinations of the coefficients 
# input:
    # currentS4, currentS8 - current number of cd4 or cd8 samples
    # current depth - depth of the graph we are in ~ how many samples have been used
    # maxDepth - maximal depth of the graph ~ how many samples do we need in total
    # coefficients - list of all the coef
def dfs_coefficients(currentS4, currentS8, depth, maxDepth, coefficients, coeff=1, cd8_count=0):
    # if we reached the desired number of samples used
    if depth == maxDepth:
        # Use cd8_count as the index to store coefficient for that category
        if cd8_count >= len(coefficients):
            coefficients.append(coeff)
        # If this category already exists, we already have its coefficient
        return 
    
    # if there are still some cd4 samples available
    if currentS4 > 0:
        # take one sample from cd4 (update currentS4)
        newS4, newS8, c = coeff_generator(currentS4, currentS8, "cd4")
        # go to the next depth of the graph
        dfs_coefficients(newS4, newS8, depth + 1, maxDepth, coefficients, coeff * c, cd8_count)
    
    if currentS8 > 0:
        newS4, newS8, c = coeff_generator(currentS4, currentS8, "cd8")
        dfs_coefficients(newS4, newS8, depth + 1, maxDepth, coefficients, coeff * c, cd8_count + 1)
    
    return



# Generate all coefficients up to max depth reachable in the dataset with cd4_only or cd8_only groups
## denominator ##
def coeff_rac_produce(s4,s8,depth,coefficients):
    # denominator is the same no matter the numerator combination
    samples = s4+s8
    denominator = 1
    for i in range(depth):
        denominator = denominator * (samples-i)

    ### list of racional numbers ###
    coeff_rac =[]
    for i in range(depth+1):
        coeff_rac.append(sp.Rational(coefficients[i],denominator))
    return coeff_rac

###################################################################################################################
### 2. PASCAL'S TRIANGLE ###
## We use Pascal's triangle values in this code for the binomial theorem
# then, we use these values in equations later
# i.e. 
       # 1
      # 1 1
     # 1 2 1    # for 2nd layer (1 cd4-cd4, 2 cd4-cd8, 1 cd8-cd8)
    # 1 3 3 1   # for 3rd layer

print("*"*30)
from math import factorial
def pascal(depth):
    coeff_pascal = []
    for j in range(depth+1):
        newCoeff = factorial(depth)//(factorial(j)*factorial(depth-j))
        print(newCoeff,end=" ")
        coeff_pascal.append(newCoeff)
    print("\n")
    return coeff_pascal

#################################################################################################################
#### 3. PROBABILITIES ####
#Pc, n4 = sp.symbols('Pc n4')


# For CD4  biased sequences
# generating different probabilities for each group
# assigning correct coefficients with Pc and (1-Pc)
# Pc is a probability of the correct choice
# (1-Pc) is the probability of a wrong choice
    # these two coefficient behave as the binomial theorem
    # for 2nd layer the situations are:
        # if the seq was cd4-biased
        # cd4-cd4 is correct twice, it's Pc^2 * (1-Pc)^0
        # cd4-cd8 is once correct, once wrong, Pc^1 * (1-Pc)^1
        # cd8-cd8 is wrong twice, Pc^0 * (1-Pc)^2
# c_pascal are coefficients for current depth (for example for 2nd layer it's 1 2 1)
# c_rac are the coeff for each situation
    # for 2nd layer it's cd4-cd4 10/15*9/14, cd4-cd8 10/15*5/14, cd8-cd8 5/15*4/14
def cd4_bias(c_pascal,c_rac,depth):
    Pc, n4 = sp.symbols('Pc n4')
    p4 = [] # stores all the possible categories 
    p4_all = 0
    for i in range(depth,-1,-1):
        newTerm = c_rac[depth-i] * Pc**i * (1-Pc)**(depth-i) * c_pascal[depth-i]
        p4.append(newTerm)
        p4_all = p4_all + newTerm
    return p4, p4_all

# For CD8 biased sequences
def cd8_bias(c_pascal,c_rac,depth):
    Pc, n4 = sp.symbols('Pc n4')
    p8 = []
    p8_all = 0
    for i in range(depth,-1,-1):
        newTerm = c_rac[depth-i] * Pc**(depth-i) * (1-Pc)**i * c_pascal[depth-i]
        p8.append(newTerm)
        p8_all = p8_all + newTerm
    return p8,p8_all




####################################################################################################################
#### 4.  EQUATIONS  ####


def eqSolving(eq1,eq2,depth,label):
    # variables, n4 = true_cd4
    Pc, n4 = sp.symbols('Pc n4')
    # Substitute from one equation to another
    n4_expr = sp.solve(eq1, n4)[0] # expression for true_cd4 from the first equation
    # Substitute into eq2
    eq2_sub = sp.simplify(eq2.subs(n4, n4_expr))
    #solutions1 = [] 
    df_sol = pd.DataFrame([])
    #labels = ["cd8_0, cd8_1, ..., cd8_s8"]
    i = 0
    solutions = [] # storing all the solutions here
    try:
    # solve returns a list of solutions
        Pc_solutions = sp.solve(eq2_sub, Pc)
        for pc_sol in Pc_solutions:
        # if the solution is real number and finite
            if pc_sol.is_real and pc_sol.is_finite:
                n4_sol = n4_expr.subs(Pc, pc_sol) # substitute for true_cd4
                Pc_rounded = round(float(pc_sol), 4)
                if 0.5 <= float(pc_sol) <= 1:
                    if df_sol.empty or not any(round(pc, 4) == Pc_rounded for pc in df_sol['Pc_val']):
                        solutions.append((float(pc_sol), float(n4_sol))) # store newly found solution
                        # new solution
                        sol_new = [{ # add new row with new solution
                            'n4': solutions[1], # calculated true_cd4
                            'Pc_val': solutions[0], # calculated Pc
                            'depth':depth, # current depth
                            'equation':label, # which two equations do we use to calculate this
                            'expected_cd4':-1, # will change later once tested
                            'expected_cd8':-1}]
                        df_new = pd.DataFrame(sol_new)
                        # add newly found solution to the data table
                        df_sol = pd.concat([df_sol,df_new],ignore_index=True)
    except:
        pass

    print("Solutions 1 (Pc, n4):")
    # trying to find new solutions using function nsolve instead solve
    for guess in [0.5,0.6,0.7,0.8,0.9]: # Pc expected in range 0.5 - 1 
        try:
            Pc_sol = sp.nsolve(eq2_sub, Pc,guess) 
            n4_sol = n4_expr.subs(Pc, Pc_sol)
            sol_Pc_float,sol_n4_float = float(Pc_sol),float(n4_sol)
            Pc_rounded = round(sol_Pc_float, 4)
            if 0.5 <= sol_Pc_float <= 1: 
                if df_sol.empty or not any(round(pc, 4) == Pc_rounded for pc in df_sol['Pc_val']):
                    print(f"  Pc = {sol_Pc_float:.10f}")
                    if Pc_sol in solutions: # if this solution has been already found by solve
                        pass
                    else: # new solution
                        sol3 = [{ # add new row with new solution
                            'n4': sol_n4_float, # calculated true_cd4
                            'Pc_val': sol_Pc_float, # calculated Pc
                            'depth':depth, # current depth
                            'equation':label, # which two equations do we use to calculate this
                            'expected_cd4':-1, # will change later once tested
                            'expected_cd8':-1}]
                        df_new = pd.DataFrame(sol3)
                        # add newly found solution to the data table
                        df_sol = pd.concat([df_sol,df_new],ignore_index=True)
        except:
            pass # failed attempt
    i = i +1
    print("*"*50)
    print(df_sol)
    print("*"*50)
    return solutions,df_sol
# main function calling the other functions in a proper order
# input:
    # chains - datasets
    # dir - root directory
def main(chains,dir):
    df_final = pd.DataFrame([])
    #print("Write the name of TCR chain: ")
    #chain = input()

    script_dir = Path(__file__).parent
    print(script_dir)
    print(dir)
    datatables_dir = dir / "datatables"
    plots_dir = dir / "plots"
    for directory in [datatables_dir,plots_dir]:
        if not os.path.exists(directory):
            os.makedirs(directory)
    #chains = ["alpha","beta"]
    #chains = ["Vb5"]
    for chain in chains:
        print(chain)
        if chain == "DeGreef":
            depth = 4
        else:
            depth = 6
        for i in range(2,depth):
            # needed initial set-up for dfs
            coefficients = []
            # Generate all coefficients up to depth 5
            input_r = input_function(i,chain,datatables_dir)
            #input_row = [maxdepth, s4, s8, all_seq] + cd8_values
            maxdepth = input_r[0]
            s4,s8 = input_r[1],input_r[2]
            all_seq = input_r[3]
            # extract the values of each cd8 category
            obs_val = input_r[4:]
            for j in range(len(obs_val)):
                if math.isnan(obs_val[j]):
                    obs_val[j] = None # instead of "-" in the line use None
                else:
                    obs_val[j] = int(obs_val[j]) # the actual numbers convert to integers
            # before ['337', '39', '10', '-', '-', '-'] - numbers are stored as string
            # after [337, 39, 10, None, None, None]
            print(obs_val)
            # generating numerator coefficients
            dfs_coefficients(s4, s8, 0, maxdepth,coefficients)
            print(coefficients)
            Pc, n4 = sp.symbols('Pc n4') # variables
            coeff_pascal = pascal(maxdepth) # pascal's triangle coeff 
            coeff_rac = coeff_rac_produce(s4,s8,maxdepth,coefficients) # generating whole coeff
            p4,p4_all = cd4_bias(coeff_pascal,coeff_rac,maxdepth) # probability eqs for cd4 biased seq
            p8,p8_all = cd8_bias(coeff_pascal,coeff_rac,maxdepth) # probability eqs for cd8 biased seq
            obs_val_index = 0
            equations = []
            # Preparing the equations
            # we want to search all the categories for current depth till we reach None
                # [337, 39, 10, None, None, None] for 2nd layer, so it will run 3 times
            while obs_val_index < len(obs_val) and obs_val[obs_val_index] != None:
                right_side_eq = obs_val[obs_val_index] # what it needs to equal, our category
                # current equation
                eq = sp.Eq(n4 * (p4[obs_val_index] / p4_all) + (all_seq - n4) * (p8[obs_val_index] / p8_all), right_side_eq)
                equations.append(eq) # add the equation to our list of equations for current depth
                obs_val_index += 1
            # Solving the equations
            for cd8_val_eq in range(1,i+1): # for all the possible eq
                method_name = "CD8_0, CD8_" + str(cd8_val_eq) # label determining which two eq we used
                solutions,df_sol = eqSolving(equations[0],equations[cd8_val_eq],maxdepth,method_name)
            # test the result by substituting back
                if df_sol.empty != True:
                    cd4_counts, cd8_counts = testing(df_sol,coefficients,input_r)
        
                    # Add the counts as new columns to df_sol
                    df_sol['cd4_cd4_calculated'] = cd4_counts
                    df_sol['cd8_cd8_calculated'] = cd8_counts
                    df_final = pd.concat([df_final,df_sol],ignore_index=True)

        print("\n" + "*"*60)
        print(df_final)
        filename = "08_"+chain[0]+"_depths_pc_n4_calculated.csv" # save our result table for each chain
        df_final.to_csv(datatables_dir/filename)

def testing(df_sol,coeffs,input_r):
# TESTING
    # input_row = [maxdepth, s4, s8, all_seq] + cd8_values
    maxdepth = input_r[0]
    s4,s8 = input_r[1],input_r[2]
    #cd4_cd4,cd8_cd8,cd8_1_mix = input_r[3],input_r[4],input_r[5]
    all_seq = input_r[3]
    obs_val = input_r[4:]
    for i in range(len(obs_val)):
        if obs_val[i] == "-" or obs_val[i] == "" or pd.isna(obs_val[i]):
            obs_val[i] = None
        else:
            obs_val[i] = float(obs_val[i])   
    cd4_cd4 = obs_val[0] if obs_val[0] is not None else 0  # CD8_0 (CD4_only)
    cd8_cd8 = obs_val[maxdepth] if len(obs_val) > maxdepth and obs_val[maxdepth] is not None else 0  # CD8_maxdepth

    Pc, n4 = sp.symbols('Pc n4')
    coeff_pascal = pascal(maxdepth)
    coeff_rac = coeff_rac_produce(s4,s8,maxdepth,coeffs)
    p4,p4_all = cd4_bias(coeff_pascal,coeff_rac,maxdepth)
    p8,p8_all = cd8_bias(coeff_pascal,coeff_rac,maxdepth)
    print("\n" + "="*50)
    print("TESTING SOLUTIONS")
    print("="*50)
    cd4_cd4_counts = []
    cd8_cd8_counts = []
    n4_val_all,Pc_val_all = df_sol['n4'], df_sol['Pc_val']
    df_sol['expected_cd4'], df_sol['expected_cd8'] = cd4_cd4, cd8_cd8 # add expected values
    #i = 0
    for i in range(len(df_sol['n4'])):
        n4_val,Pc_val = n4_val_all[i],Pc_val_all[i]
        print(f"\nTesting Solution {i+1}: Pc = {Pc_val:.6f}, n4 = {n4_val:.6f}")
    
    # Substitute the solution values into the probability expressions
        p4_val = []
        p4_all_val = 0
        p8_val = []
        p8_all_val = 0
        for j in range(maxdepth+1):
            p4_val.append(p4[j].subs(Pc, Pc_val)) # substitute the calculated values to get exact probabilites
            p4_all_val = p4_all_val + p4_val[j]

        for k in range(maxdepth+1):
            p8_val.append(p8[k].subs(Pc, Pc_val))
            p8_all_val = p8_all_val + p8_val[k]
        # calculate the fitted count based on the probabilities and fitted true_cd4 and fitted Pc
        cd4_cd4_count = n4_val * (p4_val[0] / p4_all_val) + (all_seq - n4_val) * (p8_val[0] / p8_all_val)
        cd8_cd8_count = n4_val * (p4_val[maxdepth] / p4_all_val) + (all_seq - n4_val) * (p8_val[maxdepth] / p8_all_val)
        cd4_cd4_counts.append(float(cd4_cd4_count))
        cd8_cd8_counts.append(float(cd8_cd8_count))

        print(f"  CD4-CD4 count: {float(cd4_cd4_count)} (expected: {cd4_cd4})")
        print(f"  CD8-CD8 count: {float(cd8_cd8_count)} (expected: {cd8_cd8})")
        print(f"  Error in CD4-CD4: {abs(float(cd4_cd4_count) - cd4_cd4)}")
        print(f"  Error in CD8-CD8: {abs(float(cd8_cd8_count) - cd8_cd8)}")

    return cd4_cd4_counts, cd8_cd8_counts