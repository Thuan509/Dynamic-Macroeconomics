"""
Question 3c
markov.py
-------

This code simulate and plot the Markov Chain for the AR(1) process.


"""

import os
import numpy as np
import matplotlib.pyplot as plt
from numpy import cumsum,linspace,nonzero,ones,zeros
from numpy.random import uniform


# Rouwenhorst method for discretizing AR(1)
def rouwenhorst_ar1(mu, rho, sigma, N):
    """
    Discretizes an AR(1) process using Rouwenhorst's method.
    
    Parameters:
        mu    : Intercept
        rho   : Persistence of AR(1)
        sigma : Standard deviation of shocks
        N     : Number of states
    
    Returns:
        states : Grid for AR(1) process
        P      : Transition probability matrix
    """
    p = (1 + rho) / 2
    q = p

    # Define state space
    y_mean = mu / (1 - rho)  # Unconditional mean
    y_std = sigma / np.sqrt(1 - rho**2)  # Unconditional standard deviation
    state_min = y_mean - np.sqrt(N-1) * y_std
    state_max = y_mean + np.sqrt(N-1) * y_std
    states = np.linspace(state_min, state_max, N)

    if N == 1:
        return np.array([0.0]), np.array([[1.0]])

    # Transition matrix using Rouwenhorst's method
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
mu = 0.5       # Intercept
rho = 0.85     # Persistence
sigma = 1      # Standard deviation of shocks
N = 7          # Number of states
T = 50         # Simulation length

# Generate state space and transition matrix
states, P = rouwenhorst_ar1(mu, rho, sigma, N)

#%% Simulation function.

def simulate(grid, pmat, T, seed=2025):
    """
    Simulates a Markov chain given a grid of points from the discretized process and the transition matrix.

    Parameters:
        grid : (K, N) array
            Grid of discretized state points.
        pmat : (N, N) array
            Transition probability matrix.
        T : int
            Number of periods to simulate.
        seed : int, optional (default=2025)
            Seed for reproducibility.

    Returns:
        y : (K, T) array
            Simulated series of the Markov chain.
    """

    #%% Set seed for reproducibility.
    np.random.seed(seed)

    #%% Initialize.
    N = grid.shape[1]  # Number of states.
    
    # Draw initial state from a uniform distribution over states
    state0 = np.random.choice(N, p=np.ones(N)/N)  

    #%% Simulate.
    cmat = np.cumsum(pmat, axis=1)  # CDF matrix.
    y = np.zeros((grid.shape[0], T))  # Container for simulated values.

    for t in range(T):
        y[:, t] = grid[:, state0]  # Store current state value
        rand_val = np.random.uniform()  # Draw random value
        state0 = np.searchsorted(cmat[state0, :], rand_val)  # Find next state

    #%% Plot simulation.
    plt.figure(figsize=(10, 5))
    plt.plot(range(T), y[0, :], linestyle='-', color='b', alpha=0.7)
    plt.xlabel("Time Period")
    plt.ylabel("State Value")
    plt.title("Simulated Markov Chain (50 periods)")
    plt.grid(True)
    # Save plot
    output_folder = "Output"
    os.makedirs(output_folder, exist_ok=True)
    plot_path = os.path.join(output_folder, "markov_chain_simulation.png")
    plt.savefig(plot_path)
    plt.show()

    print(f"Simulation completed. Plot saved at {plot_path}")

    return y

simulate(simulate(states.reshape(1, -1), P, T))
