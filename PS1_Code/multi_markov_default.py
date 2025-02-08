import numpy as np
import quantecon as qe
import matplotlib.pyplot as plt

# Parameters
sigma = 1       # Standard deviation of shocks
N = 7           # Number of states
T = 50          # Simulation length
rho_values = [0.75, 0.85, 0.95, 0.99]  # Different persistence values

# Set seed for reproducibility
np.random.seed(2025)

# Initialize dictionary to store simulations
simulations = {}

# Run simulations for different values of rho
for rho in rho_values:
    mc = qe.markov.approximation.rouwenhorst(N, rho, sigma)
    states = mc.state_values
    P = mc.P

    # Step 1: Draw initial state uniformly from the state vector
    initial_state_idx = np.random.choice(range(N))  
    chain = [states[initial_state_idx]]

    # Step 2: Simulate Markov Chain
    current_state_idx = initial_state_idx

    for t in range(1, T):
        u = np.random.uniform(0, 1)  # Draw from uniform [0,1]
        cumulative_probs = np.cumsum(P[current_state_idx])  # CDF of transition row
        next_state_idx = np.where(u <= cumulative_probs)[0][0]
        chain.append(states[next_state_idx])
        current_state_idx = next_state_idx

    simulations[rho] = np.array(chain)

# Plot all simulations on one graph
plt.figure(figsize=(12, 6))
for rho in rho_values:
    plt.step(range(T), simulations[rho], linestyle="-", label=f"gamma = {rho}")

plt.xlabel("Time Period")
plt.ylabel("State")
plt.title("Simulated Markov Chains for Different $\gamma$ Values")
plt.legend()
plt.grid(True)
plt.show()
