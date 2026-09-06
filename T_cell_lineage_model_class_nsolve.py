################################################################################
### LINEAGE MODEL SIMULATIONS ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

# Current script is dedicated to generate a class in which we can simulate different dataset from 2nd layer based on our calculated Pc and true_cd4 value
# This script is then used for parameter recovery validation and calculation of confidence intervals based on simulations

# loading libraries
import numpy as np
import sympy as sp


# Lineage model class
    # serves to simulate datasets and then calculate the Pc and true_cd4 value for those datasets
class TcellLineageModel:
    # main atributes of the class
    def __init__(self, Pc, n4_total, total_sequences, chain):
        self.Pc = Pc  # Commitment precision
        self.n4_total = n4_total  # Number of CD4-biased sequences, true_cd4
        self.n8_total = total_sequences - n4_total  # Number of CD8-biased sequences, true_cd8
        self.total_sequences = total_sequences # total number of sequences in 2nd layer
        
        # Coefficients from our data set-up
        # for Giogetti et al alpha and for current study experiment on Vβ5
        # where the set-up was 5 cd4 samples and 5 cd8 samples
        if chain == "a" or chain == "Vb5":
            self.c44 = 5*4/(10*9)  
            self.c48 = 5*5/(10*9) 
            self.c88 = 5*4/(10*9)   

        # for Giogetti et al beta
        # where the set-up was 13 cd4 samples and 10 cd8 samples
        if chain == "b":
            self.c44 = 13*12/(23*22)
            self.c48 = 13*10/(23*22)
            self.c88 = 10*9/(23*22)

        # for DeGreef where the set-up was 3 cd4 samples and 3 cd8 samples
        if chain == "DeGreef":
            self.c44 = 3*2/(6*5)
            self.c48 = 3*3/(6*5)
            self.c88 = 3*2/(6*5)
    
    # simulating the dataset based on the Pc and true_cd4 given
    def simulate_dataset(self, seed=None):
        """
        Each sequence is either CD4-biased or CD8-biased, determining the probabilities
        for the pair outcomes (CD4-CD4, CD4-CD8, CD8-CD8).
        """
        if seed is not None:
            np.random.seed(seed)
        
        # final number of seq in each category
        cd4_cd4_count = 0
        cd4_cd8_count = 0
        cd8_cd8_count = 0
        
        # Simulating all the sequences in 2nd layer, where each appears exactly in two samples
        for seq_idx in range(self.total_sequences):
            
            # Determine if this sequence is CD4-biased or CD8-biased
            if seq_idx < self.n4_total: # if we haven't assigned all the cd4 biased seq yet
                # CD4-biased sequence
                p_44 = self.Pc**2 * self.c44 # correct choice
                p_48 = 2*self.Pc*(1-self.Pc) * self.c48 # one mistake
                p_88 = (1-self.Pc)**2 * self.c88 # two mistakes
            else: # we have assigned all the cd4 biased seq, the rest will be cd8 biased
                # CD8-biased sequence
                p_44 = (1-self.Pc)**2 * self.c44 # two mistakes
                p_48 = 2*self.Pc*(1-self.Pc) * self.c48 # one mistake
                p_88 = self.Pc**2 * self.c88 # correct choice
            
            # Probabilities for each scenario
            total_p = p_44 + p_48 + p_88
            probs = [p_44/total_p, p_48/total_p, p_88/total_p]
            # Sample outcome for this seq (assigning to each category based on calculated probabilities)
            outcome = np.random.choice([0, 1, 2], p=probs)  # 0=CD4-CD4, 1=CD4-CD8, 2=CD8-CD8
            # raise number of seq in that particular category based on outcome
            if outcome == 0:
                cd4_cd4_count += 1
            elif outcome == 1:
                cd4_cd8_count += 1
            else:
                cd8_cd8_count += 1
        
        # Verify total count that it equals our total_sequences number
        total_count = cd4_cd4_count + cd4_cd8_count + cd8_cd8_count
        assert total_count == self.total_sequences, f"Total count {total_count} != {self.total_sequences}"        
        return cd4_cd4_count, cd4_cd8_count, cd8_cd8_count # return the combination of numbers within each category
    
# # Function to calculate the Pc and true_cd4 from the dataset simulated in the TcellLineageModel class

def fit_model_to_data(cd4_cd4_obs, cd8_cd8_obs, total_sequences, chain):
    """Fit the model to observed data using scipy optimization"""
    from scipy.optimize import fsolve
    
    # Set coefficients based on chain
    if chain == "a" or chain == "Vb5":
        c44 = 5*4 / (10*9)
        c48 = 5*5 / (10*9)
        c88 = 5*4 / (10*9)
    elif chain == "b":
        c44 = 13*12 / (23*22)
        c48 = 13*10 / (23*22)
        c88 = 10*9 / (23*22)
    elif chain == "DeGreef":
        c44 = 3*2/(6*5)
        c48 = 3*3/(6*5)
        c88 = 3*2/(6*5)
    
    def equations(vars):
        Pc_val, n4_val = vars
        
        # CD4-biased probabilities
        p4_44 = Pc_val**2 * c44
        p4_48 = 2*Pc_val*(1-Pc_val) * c48
        p4_88 = (1-Pc_val)**2 * c88
        p4_total = p4_44 + p4_48 + p4_88
        
        # CD8-biased probabilities
        p8_44 = (1-Pc_val)**2 * c44
        p8_48 = 2*Pc_val*(1-Pc_val) * c48
        p8_88 = Pc_val**2 * c88
        p8_total = p8_44 + p8_48 + p8_88
        
        # Equations
        eq1 = n4_val * (p4_44/p4_total) + (total_sequences - n4_val) * (p8_44/p8_total) - cd4_cd4_obs
        eq2 = n4_val * (p4_88/p4_total) + (total_sequences - n4_val) * (p8_88/p8_total) - cd8_cd8_obs
        
        return [eq1, eq2]
    
    valid_solutions = []
    
    # Try different initial guesses
    for Pc_guess in [0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
        n4_guess = total_sequences * 0.9  # Assume most are CD4-biased
        
        try:
            solution = fsolve(equations, [Pc_guess, n4_guess])
            Pc_sol, n4_sol = solution
            
            # Validate solution
            if 0.5 <= Pc_sol <= 1 and 0 <= n4_sol <= total_sequences:
                # Check not duplicate
                is_duplicate = any(
                    abs(existing_Pc - Pc_sol) < 1e-4 and abs(existing_n4 - n4_sol) < 1e-4
                    for existing_Pc, existing_n4 in valid_solutions
                )
                
                if not is_duplicate:
                    valid_solutions.append((round(Pc_sol, 6), round(n4_sol, 6)))
        except:
            continue
    
    return valid_solutions