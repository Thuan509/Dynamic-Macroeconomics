import numpy as np
import quantecon as qe

# Define parameters
gamma = 0.85  # Persistence parameter (γ1)
sigma = 1  # Standard deviation of shocks
N = 7  # Number of states
mu = 0.5

# Discretize AR(1) using Rouwenhorst’s method
mc = qe.markov.approximation.rouwenhorst(N, gamma, sigma, mu)

print(mc.P.shape)

print(mc.state_values.shape)

print(mc.state_values)
print(mc.P)
