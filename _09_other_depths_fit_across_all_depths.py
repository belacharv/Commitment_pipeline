################################################################################
### 08 CALCULATE PC AND TRUE_CD4 FITTING ALL THE EQ FROM ALL THE LAYERS ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate


# Current script focuses on the layer from 2 to 5
# It generates all the equations for each layer and then try to fit Pc (probability of commitment) 
    # and true_cd4 (number of the real CD4-biased sequences in the dataset, labeled as n4 for short)
    # to all of them at the same time
# it uses different methods for that and then compare results from each method


import sympy as sp
import pandas as pd
from scipy.optimize import least_squares
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import differential_evolution, basinhopping
from scipy.optimize import curve_fit
from _08_other_depths_fit_using_two_eq import cd4_bias
from _08_other_depths_fit_using_two_eq import cd8_bias
from _08_other_depths_fit_using_two_eq import input_function
from _08_other_depths_fit_using_two_eq import dfs_coefficients
from _08_other_depths_fit_using_two_eq import coeff_rac_produce
from _08_other_depths_fit_using_two_eq import pascal
import math
from pathlib import Path
## previous three steps are in the same as in other_depths_all_the_eq_considered
####################################################################################################################
#### 4.  EQUATIONS  ####

# function calculating residues based on given guesses - MODIFIED FOR ALL DEPTHS
def residuals_all_depths(params, all_depths_data):
    """
    params: [Pc, n4_depth2, n4_depth3, n4_depth4, n4_depth5]
    all_depths_data: list of dicts containing data for each depth
    """
    Pc_val = params[0]
    res = []
    
    for current_depth, depth_data in enumerate(all_depths_data):
        n4_val = params[current_depth + 1]
        
        # creating mathematical functions for each scenario (ranging from CD4 max CD8 0 to CD4 0 CD8 max)
        # when the sequence is cd4 biased (p4) or cd8 biased (p8)
        # p4_all and p8_all are functions suming up all the equations from p4 or p8
        p4_funcs = depth_data['p4_funcs']
        p8_funcs = depth_data['p8_funcs']
        p4_all_func = depth_data['p4_all_func']
        p8_all_func = depth_data['p8_all_func']
        all_seq = depth_data['all_seq']
        obs_val = depth_data['obs_val']
        
        # parameters observed from the initial guesses
        p4_all_val = p4_all_func(Pc_val)
        p8_all_val = p8_all_func(Pc_val)
        # for each of the observed values (=number of the sequences adressing that scenario)   
        # for example for depth 2 in Giorgetti el al. alpha it's
            # CD4-CD4 258
            # CD4-CD8 47
            # CD8-CD8 10  
        # calculate residuals
        for i, obs in enumerate(obs_val):
            if obs is not None:
                p4_val = p4_funcs[i](Pc_val)
                p8_val = p8_funcs[i](Pc_val)
                predicted = n4_val * (p4_val / p4_all_val) + (all_seq - n4_val) * (p8_val / p8_all_val)
                res.append(abs(predicted - obs))
    return res

# if residuals are not returned as squared
def objective(params, all_depths_data):
    res = residuals_all_depths(params, all_depths_data)
    return np.sum(np.array(res)**2)  # Sum of squared residuals

#########################################################################################################
## different algorithms that could be used for the fit ##

# applying least squares formula method
# fit the variables to the equations and then calculates residuals
# for more information: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html 
def least_squares_fit(initial_guess,all_depths_data,method,lower_bounds,upper_bounds):
    result = least_squares(
        residuals_all_depths,
            initial_guess,
            args=(all_depths_data,),
            # used lists created according to all_seq from each parameter 
            bounds=(lower_bounds, upper_bounds),
            method= method
        )            
    return result       

# Minimize a scalar function of the variables using the L-BFGS-B algorithm 
    # (Limited-memory Broyden–Fletcher–Goldfarb–Shanno algorithm)
    # uses Hessian matrix, but it doesn't store the full matrix Hk
        # it constructs Hk implicitly using a small number of vectors for the last few iterations
# for more information: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html 
def minimize_fit(initial_guess,all_depths_data):
    result = minimize(
                    objective,
                    initial_guess,
                    args=(all_depths_data,),
                    method='L-BFGS-B',  # Supports bounds
                    bounds=[(0.5, 1)] + [(0, d['all_seq']) for d in  all_depths_data]
                )
    return result

# Differential evolution algorithm made by Storn and Price
    # Finds the global minimum of a multivariate function
# for more information: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html
def diff_evolution(all_depths_data):
    result = differential_evolution(
        objective,
        bounds=[(0.5, 1)] + [(0, d['all_seq']) for d in all_depths_data],
        args=(all_depths_data,),
        strategy='best1bin', # default
        maxiter=1000
    )
    return result

# Basin-hopping algorithm
    # two-phase method that combines a global stepping algorithm with local minimization at each step.
# for more information: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html 
def basinhopping_fit(initial_guess,all_depths_data):
    result = basinhopping(
        objective,
        initial_guess,
        minimizer_kwargs={
           'method': 'L-BFGS-B', # algorithm explained earlier
          'bounds': [(0.5, 1)] + [(0, d['all_seq']) for d in all_depths_data],
          'args': (all_depths_data,)
        }
    )
    return result


### Curve_fit method ###
## A curve fitting method studies the relationship between predictors (independent variables) and response variable (dependent variable)
# The principle relies on construction of a mathematical function that fits the best to a series of data points
# parameters needed for a curve fit:
    # f - model function that takes independent variables and parameters to fit
    # xdata - independent variables
    # ydata -  dependent variables

# constructing the ydata vector with all the observed values from all the depths
# stores index references to compute predictions in the curve-fit model function
def build_curvefit_data(all_depths_data):   
    ydata = []
    ref_index = []  # stores (current_depth, cd8_number) tuple where 
        # current depth is the layer we're in
        # cd8_number belogns to each of the scenarios (cd8_0,cd8_1...cd8_max)
    # in current depth
    for current_depth, depth_data in enumerate(all_depths_data):
        # extract the number of events for each scenario (cd8_0,cd8_1...cd8_max)
        for cd8_number, obs in enumerate(depth_data["obs_val"]):
            if obs is not None: # if we haven't reached our limit
                ydata.append(obs) # add value to the list of observed values
                ref_index.append((current_depth, cd8_number)) # store the indeces

    return np.array(ydata, dtype=float), ref_index

# construct the vector with predicted results
# result of our fit
def model_curve_fit(ref_index, all_depths_data, params):
    Pc_val = params[0] # independent parameter
    y_pred = [] # list of predicted results, will be compared with ydata and determinate the residuals

    for current_depth, cd8_number in ref_index:
        depth_data = all_depths_data[current_depth] # extract data
        n4_val = params[current_depth + 1] # extract true_cd4 parameter for current depth
        all_seq = depth_data["all_seq"]

        # create model functions
        p4_val = depth_data["p4_funcs"][cd8_number](Pc_val)
        p8_val = depth_data["p8_funcs"][cd8_number](Pc_val)
        p4_all_val = depth_data["p4_all_func"](Pc_val)
        p8_all_val = depth_data["p8_all_func"](Pc_val)
        # result value of our prediction
        predicted = n4_val * (p4_val / p4_all_val) + (all_seq - n4_val) * (p8_val / p8_all_val)
        y_pred.append(predicted)

    return np.array(y_pred, dtype=float)


def eqSolving_all_depths(all_depths_data, chain, method):
    #Fit Pc and n4 values for all depths simultaneously
    solutions = []
    df_sol = pd.DataFrame()
    # parameters used
        # Pc, n4_depth2, n4_depth3, n4_depth4, n4_depth5
    
    # Bounds for each depth and for P
    lower_bounds = [0]  # Pc lower bound
    upper_bounds = [1]  # Pc upper bound
    for depth_data in all_depths_data:
        all_seq = depth_data['all_seq']
        lower_bounds.append(0)
        upper_bounds.append(all_seq)
    
    ## Initial guess for all the parameters
    # Initial guesses for least_squares function, we expect the Pc somewhere in range 0.7-0.95
    Pc_guesses = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]   
    # for each of these Pc values we guess all the other parameters (true_cd4 for each depth) 
    for Pc_guess in Pc_guesses:
        initial_guess = [Pc_guess]
        for depth_data in all_depths_data:
            initial_guess.append(depth_data['all_seq'] * Pc_guess) # 
        
        try:
            # apply each different method to fit the two parameters for multiple equations
            if method == "minimize":
                result = minimize_fit(initial_guess,all_depths_data)
            elif method == "diff_ev":
                result = diff_evolution(all_depths_data)
            elif method == "basinhopping":
                result = basinhopping_fit(initial_guess,all_depths_data)
            elif method == "curve_fit":
                # initialize the arrays with observed values (ydata) and indeces for each event
                    # ref_index stores current depth and the index of a scenario
                ydata, ref_index = build_curvefit_data(all_depths_data)
                bounds_lower = lower_bounds
                bounds_upper = upper_bounds

                popt, pcov = curve_fit(
                    lambda _, *params: model_curve_fit(ref_index, all_depths_data, params),
                    xdata=np.zeros(len(ydata)), 
                    ydata=ydata,
                    p0=initial_guess,
                    bounds=(bounds_lower, bounds_upper),
                    maxfev=20000
                )
                result = type("obj", (), {"success": True, "x": popt, "fun": ydata - model_curve_fit(ref_index, all_depths_data, popt)})

            else:
                result = least_squares_fit(initial_guess,all_depths_data,method,lower_bounds,upper_bounds)

            if result.success:
                params = result.x
                Pc_val = params[0]
                
                # Check if this solution is already found
                is_duplicate = False
                for existing_sol in solutions:
                    if all(abs(existing_sol[i] - params[i]) < 1e-3 for i in range(len(params))):
                        is_duplicate = True
                        print(f"duplicate {existing_sol}")
                        break
                # if the result is new and Pc value makes sense
                if not is_duplicate and 0 <= Pc_val <= 1:
                    solutions.append(tuple(float(p) for p in params))
                    
                    # Build solution dictionary
                    sol_dict = {'Pc_val': float(Pc_val)}
                    
                    for current_depth, depth_data in enumerate(all_depths_data):
                        depth = depth_data['depth']
                        n4_val = params[current_depth + 1]
                        all_seq = depth_data['all_seq']
                        sol_dict[f'n4_depth{depth}'] = float(n4_val)
                        sol_dict[f'all_seq_depth{depth}'] = all_seq
                        sol_dict[f'depth{depth}'] = depth
                    
                    
                    sol_dict['chain'] = chain
                    sol_dict['method'] = method
                    if method == "minimize" or method == "diff_ev" or method == "basinhopping":
                        residual_sum = float(result.fun)   # minimize already returns sum-of-squared-residuals
                    else:
                        residual_sum = float(np.sum(np.array(result.fun)**2))
                    sol_dict['residual_sum'] = residual_sum
                    
                    df_new = pd.DataFrame([sol_dict])
                    df_sol = pd.concat([df_sol, df_new], ignore_index=True)
        except Exception as e:
            print(f"Failed for Pc guess {Pc_guess}: {e}")
            pass
    
    print(f"Solutions for {chain} chain (all depths):")
    for s in solutions:
        print(f"  Pc = {s[0]:.6f}, n4 values = {s[1:]}")
    
    print("*"*50)
    print(df_sol)
    print("*"*50)
    return solutions, df_sol

# input function processing one row
def process_input(input_r,all_depths_data):
    coefficients = []
    visited = set()
    # input row looks like this
            #input_row = [maxdepth, s4, s8, all_seq] + cd8_values
    # assign each variable from input row
    maxdepth = input_r[0]
    s4, s8 = input_r[1], input_r[2]
    all_seq = input_r[3]
    
    # number of cd8 in each category
    # we now got ['258','47','10','-','-','-'] and we want to convert the row to integers
    obs_val = input_r[4:]
    for j in range(len(obs_val)):
        if math.isnan(obs_val[j]):
            obs_val[j] = None # instead of "-" in the line use None
        else:
            obs_val[j] = int(obs_val[j]) # the actual numbers convert to integers

    dfs_coefficients(s4, s8, 0, maxdepth, coefficients, visited)
    print(list(coefficients))
    Pc, n4 = sp.symbols('Pc n4')
    coeff_pascal = pascal(maxdepth)
    coeff_rac = coeff_rac_produce(s4, s8, maxdepth, coefficients)
    p4, p4_all = cd4_bias(coeff_pascal, coeff_rac, maxdepth)
    p8, p8_all = cd8_bias(coeff_pascal, coeff_rac, maxdepth)
        
    # Convert to numerical functions
    p4_funcs = [sp.lambdify(Pc, expr, 'numpy') for expr in p4]
    p8_funcs = [sp.lambdify(Pc, expr, 'numpy') for expr in p8]
    p4_all_func = sp.lambdify(Pc, p4_all, 'numpy')
    p8_all_func = sp.lambdify(Pc, p8_all, 'numpy')
    # add data from this row to our global input table  
    all_depths_data.append({
        'depth': maxdepth,
        'p4_funcs': p4_funcs,
        'p8_funcs': p8_funcs,
        'p4_all_func': p4_all_func,
        'p8_all_func': p8_all_func,
        'all_seq': all_seq,
        'obs_val': obs_val,
        'input_r': input_r,
        'coefficients': coefficients
    })

def main(chains,dir):
    script_dir = Path(__file__).parent
    print(script_dir)
    print(dir)
    datatables_dir = dir / "datatables"
    df_final = pd.DataFrame([])
    methods = ["curve_fit","trf","dogbox","minimize","diff_ev","basinhopping"]
    for chain in chains:
        all_depths_data = []
        if chain == "DeGreef":
            maxdepth = 4
        else:
            maxdepth = 6
        for i in range(2, maxdepth):
            # Generate all coefficients up to depth 5
            input_r = input_function(i, chain,datatables_dir)
            process_input(input_r,all_depths_data)
            # input row looks like this
                #input_row = [maxdepth, s4, s8, all_seq] + cd8_values
        
        for method in methods:
        # Now fit all depths together
            solutions, df_sol = eqSolving_all_depths(all_depths_data, chain,method)
            print(len(solutions))
            df_final = pd.concat([df_final, df_sol], ignore_index=True)

    
    print(df_final)
    filename = "09_"+chain+"_all_depths_pc_n4_calculated.csv"
    df_final.to_csv(datatables_dir/filename)