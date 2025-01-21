"""

3b_rouwen.py
----------

This code discretizes a Markov process into N states using the Rouwenhorst's Method.

"""

import math
import os
import numpy as np
from scipy.stats import norm
from numpy import arange,count_nonzero,expand_dims,identity,linalg,linspace,nonzero,ones,prod,tile,zeros


def rouwenhorst_ar1(N, mu, rho, sigma):
    """
    This function discretizes an AR(1) process.
    
            y(t) = mu + gamma*y(t-1) + eps(t), eps(t) ~ NID(0,sigma^2)
        =>  y(t) = 0.5 + 0.85*y(t-1) + eps(t), eps(t) ~ NID(0,1)
    
    Input:
        mu    : Intercept of AR(1).
        rho   : Persistence of AR(1).
        sigma : Standard deviation of error term.
        N     : Number of states.
        m     : Parameter such that m time the unconditional std. dev. of the AR(1) is equal to the largest grid point.
        
    Output:
        states    : Grid for the AR(1) process.
        P2 : Transition probability matrix.
    """
    
    # Step 1: Construct state vector
    # Calculate unconditional mean and standard deviation
    unconditional_mean = mu / (1 - rho)
    unconditional_std = sigma / np.sqrt(1 - rho**2)
    
    # Create equally spaced grid
    state_min = unconditional_mean - np.sqrt(N-1) * unconditional_std
    state_max = unconditional_mean + np.sqrt(N-1) * unconditional_std
    states = linspace(state_min, state_max, N)
    
    # Step 2: Set up initial transition matrix (2x2)
    p = q = (1 + rho) / 2
    P2 = np.array([[p, 1-p],
                [1-q, q]])
    
    if N == 2:
        return states, P2
    
    # Step 3: Recursive construction of transition matrices
    for n in range(3, N+1):
        P_prev = P2
        P = zeros((n, n))
        
        # Construct using the recursive formula
        P[:-1,:-1] += p * P_prev
        P[:-1,1:] += (1-p) * P_prev
        P[1:,:-1] += (1-q) * P_prev
        P[1:,1:] += q * P_prev
        
        # Scale rows (except first and last)
        P[1:-1,:] /= 2
        
        P2 = P
    
    return states, P2

# Parameters
N = 7
mu = 0.5
rho = 0.85
sigma = 1

# N, mu, rho, sigma
states, P2 = rouwenhorst_ar1(N, mu, rho, sigma)
# Print the results
print("Grid points (states):")
print(states)
print("\nTransition probability matrix (P2):")
print(P2) 

## Define the output folder and file
output_folder = "Output"
output_file = os.path.join(output_folder, "rouwen_output.txt")

# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Save the results to a text file
with open(output_file, "w") as file:
    file.write("Grid Points (states):\n")
    file.write(str(states) + "\n\n")
    file.write("Transition Probability Matrix (P2):\n")
    for row in P2:
        file.write(" ".join(f"{val:.6f}" for val in row) + "\n")

print("Grid points (states):")
print(states)
print("\nTransition probability matrix (P2):")
print(P2)
print(f"Results saved to {output_file}")