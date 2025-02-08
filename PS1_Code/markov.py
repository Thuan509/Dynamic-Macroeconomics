"""

markov.py
-------

This code simulate and plot the Markov Chain for the AR(1) process.


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

# Set seed for reproducibility
np.random.seed(2025)

# Draw initial state uniformly from the state vector
initial_state_idx = np.random.choice(range(N))  # Uniformly select index
mc = [states[initial_state_idx]]

# Simulate Markov Chain
current_state_idx = initial_state_idx

for t in range(1, T):
    u = np.random.uniform(0, 1)  # Draw from uniform [0,1]
    cumulative_probs = np.cumsum(P[current_state_idx])  # CDF of transition row

    # Determine the next state based on cumulative probabilities
    next_state_idx = np.where(u <= cumulative_probs)[0][0]

    mc.append(states[next_state_idx])
    current_state_idx = next_state_idx

# Convert to NumPy array for plotting
mc = np.array(mc)

# Plot the Markov Chain simulation
plt.figure(figsize=(10, 5))
plt.step(range(T), mc, linestyle="-", color="b")
plt.xlabel("Time Period")
plt.ylabel("State")
plt.title("Simulated Markov Chain (50 periods)")
plt.grid(True)

# Save plot
output_folder = "Output"
os.makedirs(output_folder, exist_ok=True)
plot_path = os.path.join(output_folder, "markov_chain_simulation.png")
plt.savefig(plot_path)
plt.show()

print(f"Simulation completed. Plot saved at {plot_path}")