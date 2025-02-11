import numpy as np
import quantecon as qe
import matplotlib.pyplot as plt

# Parameters
rho = 0.85      # Persistence of AR(1) process
sigma = 1       # Standard deviation of shocks
N = 7           # Number of states
T = 50          # Simulation length
mu = 0.5        # Intercept

# Generate Markov Chain using Quantecon
mc = qe.markov.approximation.rouwenhorst(N, rho, sigma, mu=mu)

# Extract state values (grid points) and transition matrix
states = mc.state_values  # Shift by mu if necessary
P = mc.P


# Set seed for reproducibility
np.random.seed(2025)

# Step 1: Draw initial state uniformly from the state vector
initial_state_idx = np.random.choice(range(N))  # Uniformly select index
chain = [states[initial_state_idx]]

# Step 2: Simulate Markov Chain
current_state_idx = initial_state_idx

for t in range(1, T):
    u = np.random.uniform(0, 1)  # Draw from uniform [0,1]
    cumulative_probs = np.cumsum(P[current_state_idx])  # CDF of transition row

    # Determine the next state based on cumulative probabilities
    next_state_idx = np.where(u <= cumulative_probs)[0][0]

    chain.append(states[next_state_idx])
    current_state_idx = next_state_idx

# Convert to NumPy array for plotting
chain = np.array(chain)

# Plot the Markov Chain simulation
plt.figure(figsize=(10, 5))
plt.step(range(T), chain, linestyle="-", color="b", where='mid')
plt.xlabel("Time Period")
plt.ylabel("State")
plt.title("Simulated Markov Chain (50 periods)")
plt.grid(True)
plt.show()
