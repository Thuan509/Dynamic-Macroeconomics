"""
Question 3b
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
        states : Grid for the AR(1) process.
        P      : Transition probability matrix.
    """
    
    # Step 1: Calculate parameters
    p = (1 + rho) / 2
    q = p
    
    # Step 2: Create state vector
    y_mean = mu / (1 - rho)  # unconditional mean
    y_std = sigma / np.sqrt(1 - rho**2)  # unconditional standard deviation
    state_min = y_mean - np.sqrt(N-1) * y_std
    state_max = y_mean + np.sqrt(N-1) * y_std
    states = np.linspace(state_min, state_max, N) # Define the grid
    
    if N == 1:
        return np.array([0.0]), np.array([[1.0]])

    mat = np.array([[p, 1 - p], [1 - q, q]])  # Base case (N=2)

    for n in range(3, N + 1):
        new_mat = np.zeros((n, n))
        new_mat[:-1, :-1] += p * mat
        new_mat[:-1, 1:] += (1 - p) * mat
        new_mat[1:, :-1] += (1 - q) * mat
        new_mat[1:, 1:] += q * mat
        new_mat[1:-1, :] /= 2
        mat = new_mat

    return states, mat

# Parameters
N = 7
mu = 0.5
rho = 0.85
sigma = 1

# N, mu, rho, sigma
states, P = rouwenhorst_ar1(N, mu, rho, sigma)
# Print the results
print("The shape of the state vector:", states.shape)
print("Grid points (states):")
print(states)
print("The shape of the transition probability matrix:", P.shape)
print("\nTransition probability matrix (P):")
print(P) 


## Define the output folder and file
output_folder = "Output"
output_file = os.path.join(output_folder, "rouwen_output.txt")

# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Save the results to a text file
with open(output_file, "w") as file:
    file.write(f"Shape of the state vector: {states.shape}\n")
    file.write("Grid Points (states):\n")
    file.write(str(states) + "\n\n")

    file.write(f"Shape of the transition probability matrix: {P.shape}\n")
    file.write("Transition Probability Matrix (P):\n")
    for row in P:
        file.write(" ".join(f"{val:.6f}" for val in row) + "\n")

# Print results
print(f"Shape of the state vector: {states.shape}")
print("Grid Points (states):")
print(states)

print(f"\nShape of the transition probability matrix: {P.shape}")
print("Transition Probability Matrix (P):")
print(P)

print(f"\nResults saved to {output_file}")