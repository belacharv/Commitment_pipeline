################################################################################
### 04 CALCULATING PC AND TRUE CD4 IN 2ND LAYER ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate


## The final sequences that are not unique (filtered in first step and while creating count table)
# have been collected and seperated into layers based on in how many samples they are present.
# In experiment alpha, where they edited alpha chain, they had 15 samples from spleen, 10 CD4 and 5 CD8.
# Each sequence have been assigned to one of the 14 layers (ranging from 2 to 15) based on in how many samples it is present
# In current python script, we focus only on the 2nd layer, hence sequences present in exactly two samples

### This script sets up the goal of calculating Pc (probability of commitment) and true_cd4 (number of CD4-biased sequences)
# these are the two variables: Pc and true_cd4
# and we calculate it using equation of probability to see CD4-CD4 (one seq to be present in two CD4 samples) and CD8-CD8


import sympy as sp
import os

script_dir_os = os.path.dirname(os.path.abspath(__file__))
datatables_dir_os = os.path.join(script_dir_os, "datatables")
plots_dir_os = os.path.join(script_dir_os, "plots")
for directory in [datatables_dir_os,plots_dir_os]:
    if not os.path.exists(directory):
        os.makedirs(directory)

# with hashtag displayed number for alpha
print("Number of CD4 samples:",end=" ") # 5
n = int(input()) 
print("Number of CD8 samples:",end=" ") # 5
m = int(input())
print("Number of CD4-CD4 sequences:",end=" ") # 258
cd4_cd4 = int(input())
print("Number of CD8-CD8 sequences:",end=" ") # 10
cd8_cd8 = int(input())
all = n+m
print("Number of all seq:",end=" ") # 315
all_seq = int(input())
# for beta the values are
# 21, 10, 89, 4, 114 
# in this exact order

# defining our two variables
# Pc = probability of commitment
# n4 = true_cd4
Pc, n4 = sp.symbols('Pc n4') 
# constants for each combination
    # c44 = constant CD4-CD4
def generating_constants(depth):
    const = [None for i in range(depth)]
    print(const)
    return const
c44 = sp.Rational(n*(n-1), all*(all-1))  
c48 = sp.Rational(n*m, all*(all-1))   
c88 = sp.Rational(m*(m-1), all*(all-1))    
# Probability to go by one of the 3 possible cases 
    # cd4-cd4 marked as 44
    # cd4-cd8 marked as 48
    # cd8-cd8 marked as 88
# For CD4 biased sequences
p4_44 = Pc**2 * c44 # both of the T cells chose the right fate (according to the bias of the sequence)
p4_48 = 2*Pc*(1-Pc) * c48 # one was correct, one was wrong
p4_88 = (1-Pc)**2 * c88 # both were wrong
# For CD8 biased sequences
p8_44 = (1-Pc)**2 * c44
p8_48 = 2*Pc*(1-Pc) * c48
p8_88 = Pc**2 * c88

# Equations
eq1 = sp.Eq(n4 * (p4_44 / (p4_44 + p4_48 + p4_88)) + (all_seq - n4) * (p8_44 / (p8_44 + p8_48 + p8_88)), cd4_cd4) # only cd4 samples
eq2 = sp.Eq(n4 * (p4_88 / (p4_44 + p4_48 + p4_88)) + (all_seq - n4) * (p8_88 / (p8_44 + p8_48 + p8_88)), cd8_cd8) # only cd8 samples


# Substitute from one equation to another
# expression for n4 in terms of Pc
n4_expr_from_eq1 = sp.solve(eq1, n4)[0] # taking the first solution of eq1 for true_cd4
# substitute the expression for n4 into eq2
eq2_sub = sp.simplify(eq2.subs(n4, n4_expr_from_eq1))


solutions = []
# nsolve equires an initial guess around which it will look for exact solutions
for guess in [0.01, 0.9]:  # we expect Pc to be close to 1 and other solution to be close to 0
    try:
        Pc_solution = sp.nsolve(eq2_sub, Pc, guess) # solve the second equation using substitution from the first and find a solution within the range
        n4_solution = n4_expr_from_eq1.subs(Pc, Pc_solution) # substitute Pc to our equation expression for n4
        solutions.append((float(Pc_solution), float(n4_solution)))
    except:
        pass

print("Solutions (Pc, n4):")
for s in solutions:
    print(s)

# TESTING
print("\n" + "*"*50)
print("TESTING SOLUTIONS")
print("*"*50)
# for all the solutions we found
for i, (Pc_val, n4_val) in enumerate(solutions):
    print(f"\nTesting Solution {i+1}: Pc = {Pc_val:.6f}, n4 = {n4_val:.6f}")
    
    # Substitute the solution values into the probability expressions
    # CD4 biased seq
    p4_44_val = p4_44.subs(Pc, Pc_val)
    p4_48_val = p4_48.subs(Pc, Pc_val)
    p4_88_val = p4_88.subs(Pc, Pc_val)
    # CD8 biased seq
    p8_44_val = p8_44.subs(Pc, Pc_val)
    p8_48_val = p8_48.subs(Pc, Pc_val)
    p8_88_val = p8_88.subs(Pc, Pc_val)
    
    # Calculate CD4-CD4 count (should equal to cd4_cd4)
    cd4_cd4_count = n4_val * (p4_44_val / (p4_44_val + p4_48_val + p4_88_val)) + (all_seq - n4_val) * (p8_44_val / (p8_44_val + p8_48_val + p8_88_val))
    
    # Calculate CD8-CD8 count (should equal to cd8_cd8)
    cd8_cd8_count = n4_val * (p4_88_val / (p4_44_val + p4_48_val + p4_88_val)) + (all_seq - n4_val) * (p8_88_val / (p8_44_val + p8_48_val + p8_88_val))
    # printing results of the test
    print(f"  CD4-CD4 count: {float(cd4_cd4_count)} (expected: {cd4_cd4})")
    print(f"  CD8-CD8 count: {float(cd8_cd8_count)} (expected: {cd8_cd8})")
    print(f"  Error in CD4-CD4: {abs(float(cd4_cd4_count) - cd4_cd4)}")
    print(f"  Error in CD8-CD8: {abs(float(cd8_cd8_count) - cd8_cd8)}")