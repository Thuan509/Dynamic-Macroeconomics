"""
Question 3d
multi_markov.py
-------

This code simulates and plots the Markov Chain for the AR(1) process for different persistence values.

"""

import os
import numpy as np
import matplotlib.pyplot as plt

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
    state_min = y_mean - np.sqrt(N - 1) * y_std
    state_max = y_mean + np.sqrt(N - 1) * y_std
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

# Simulation function
def simulate(grid, pmat, T, seed=2025):
    """
    Simulates a Markov chain given a grid of points from the discretized process and the transition matrix.

    Parameters:
        grid : (1, N) array
            Grid of discretized state points.
        pmat : (N, N) array
            Transition probability matrix.
        T : int
            Number of periods to simulate.
        seed : int, optional (default=2025)
            Seed for reproducibility.

    Returns:
        y : (T,) array
            Simulated series of the Markov chain.
    """
    np.random.seed(seed)

    N = grid.shape[1]  # Number of states

    # Draw initial state from a uniform distribution
    state0 = np.random.choice(N, p=np.ones(N) / N)

    cmat = np.cumsum(pmat, axis=1)  # CDF matrix
    y = np.zeros(T)  # Store simulated values

    for t in range(T):
        y[t] = grid[0, state0]  # Store current state value
        rand_val = np.random.uniform()  # Draw random value
        state0 = np.searchsorted(cmat[state0, :], rand_val)  # Find next state

    return y

# Parameters
sigma = 1       # Standard deviation of shocks
N = 7           # Number of states
T = 50          # Simulation length
rho_values = [0.75, 0.85, 0.95, 0.99]  # Different persistence values
mu = 0.5

# Set seed for reproducibility
np.random.seed(2025)

# Dictionary to store simulations
simulations = {}

# Store initial state spaces
state_spaces = {}

# Run simulations for different values of rho
for rho in rho_values:
    states, P = rouwenhorst_ar1(mu, rho, sigma, N)
    y_sim = simulate(states.reshape(1, -1), P, T)  # Reshape for consistency
    simulations[rho] = y_sim  # Store simulation result
    state_spaces[rho] = states  # Store initial state space

    # Print the initial state space
    print(f"Initial state space for gamma = {rho}:")
    print(states, "\n")

# Plot all simulations on one graph
plt.figure(figsize=(10, 5))
for rho in rho_values:
    plt.plot(range(T), simulations[rho], linestyle="-", label=f"$\gamma$ = {rho}")

plt.xlabel("Time Period")
plt.ylabel("State")
plt.title("Simulated Markov Chains for Different $\gamma$ Values")
plt.legend()
plt.grid(True)

# Save plot
output_folder = "Output"
os.makedirs(output_folder, exist_ok=True)
plot_path = os.path.join(output_folder, "multi_markov_simulation.png")
plt.savefig(plot_path)
plt.show()

# Save state spaces to a text file
state_space_file = os.path.join(output_folder, "state_spaces.txt")

with open(state_space_file, "w", encoding="utf-8") as f:
    for rho, states in state_spaces.items():
        f.write(f"gamma = {rho}\n")
        f.write(np.array2string(states, precision=4) + "\n\n")

print(f"State spaces saved in {state_space_file}")        
print(f"Simulation completed. Plot saved at {plot_path}")
