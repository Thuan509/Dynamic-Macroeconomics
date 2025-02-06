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
    
    # Step 1: Calculate parameters
    p = (1 + rho) / 2
    q = p
    
    # Step 2: Create state vector
    y_mean = mu / (1 - rho)  # unconditional mean
    y_std = sigma / np.sqrt(1 - rho**2)  # unconditional standard deviation
    
    state_min = y_mean - np.sqrt(N-1) * y_std
    state_max = y_mean + np.sqrt(N-1) * y_std
    states = np.linspace(state_min, state_max, N)
    
    # Step 3: Build transition matrix recursively
    # Start with 2x2 matrix
    P = np.array([[p, 1-p], [1-q, q]])
    
    # Expand to NxN
    for n in range(3, N+1):
        P_old = P
        P = np.zeros((n, n))
        
        # Fill new transition matrix
        P[:-1,:-1] += p * P_old
        P[:-1,1:] += (1-p) * P_old
        P[1:,:-1] += (1-q) * P_old
        P[1:,1:] += q * P_old
        
        # Scale rows (except first and last)
        P[1:-1,:] /= 2
        return states, P

# Parameters
N = 7
mu = 0.5
rho = 0.85
sigma = 1

# N, mu, rho, sigma
states, P = rouwenhorst_ar1(N, mu, rho, sigma)
# Print the results
print("Grid points (states):")
print(states)
print("\nTransition probability matrix (P2):")
print(P) 

## Define the output folder and file
output_folder = "Output"
output_file = os.path.join(output_folder, "rouwen_output.txt")

# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Save the results to a text file
with open(output_file, "w") as file:
    file.write("Grid Points (states):\n")
    file.write(str(states) + "\n\n")
    file.write("Transition Probability Matrix (P):\n")
    for row in P:
        file.write(" ".join(f"{val:.6f}" for val in row) + "\n")

print("Grid points (states):")
print(states)
print("\nTransition probability matrix (P):")
print(P)
print(f"Results saved to {output_file}")