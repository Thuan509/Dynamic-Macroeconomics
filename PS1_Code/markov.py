import numpy as np
from numpy.random import uniform
import random

def simulate_markov_chain(P, states, T, initial_state_idx):
    """
    Simulate Markov chain for T periods
    P: transition matrix
    states: state vector
    T: number of periods
    initial_state_idx: index of initial state
    """
    path = np.zeros(T)
    current_state = initial_state_idx
    
    for t in range(T):
        path[t] = states[current_state]
        # Draw next state using transition probabilities
        current_state = np.random.choice(len(states), p=P[current_state])
    
    return path

# Set seed
np.random.seed(2025)
random.seed(2025)

# Parameters
T = 50  # number of periods
gamma_values = [0.75, 0.85, 0.95, 0.99]
mu = 0.5
sigma = 1

# Initialize storage for all paths
all_paths = {}

# Draw initial state (uniform distribution)
initial_state_idx = np.random.randint(0, 7)  # 7 states

# Generate paths for each gamma
for gamma in gamma_values:
    # Calculate states and transition matrix using Rouwenhorst method
    unconditional_mean = mu / (1 - gamma)
    unconditional_std = sigma / np.sqrt(1 - gamma**2)
    
    # Create state vector
    N = 7
    state_min = unconditional_mean - np.sqrt(N-1) * unconditional_std
    state_max = unconditional_mean + np.sqrt(N-1) * unconditional_std
    states = np.linspace(state_min, state_max, N)
    
    # Create transition matrix (simplified for example)
    p = q = (1 + gamma) / 2
    P = np.zeros((N, N))
    # ... (transition matrix construction code as before)
    
    # Simulate path
    path = simulate_markov_chain(P, states, T, initial_state_idx)
    all_paths[gamma] = path.tolist()

# Print first few values of each path
for gamma, path in all_paths.items():
    print(f"Gamma = {gamma}: {path[:5]}")