import numpy as np

# Define parameters
gamma1 = 0.85  # Persistence
N = 7  # Number of states
mu = 0.5  # Mean term
sigma = 1  # Standard deviation

# Compute state space limits
y_min = mu - (sigma * np.sqrt(N-1)) / (1 - gamma1)
y_max = mu + (sigma * np.sqrt(N-1)) / (1 - gamma1)

# Create state vector (equally spaced grid)
state_vector = np.linspace(y_min, y_max, N)

# Define Rouwenhorst transition matrix function
def rouwenhorst(N, p, q):
    if N == 2:
        return np.array([[p, 1 - p],
                         [1 - q, q]])
    else:
        Pi_Nm1 = rouwenhorst(N-1, p, q)  # Recursive call
        
        # Construct new transition matrix
        Pi_N = p * np.block([[Pi_Nm1, np.zeros((N-1, 1))],
                             [np.zeros((1, N-1)), 0]]) + \
               (1 - p) * np.block([[np.zeros((N-1, 1)), Pi_Nm1],
                                   [0, np.zeros((1, N-1))]]) + \
               q * np.block([[0, np.zeros((1, N-1))],
                             [np.zeros((N-1, 1)), Pi_Nm1]]) + \
               (1 - q) * np.block([[np.zeros((1, N-1)), 0],
                                   [Pi_Nm1, np.zeros((N-1, 1))]])

        # Normalize middle rows
        Pi_N[1:-1, :] /= 2
        return Pi_N

# Compute transition matrix
p = (1 + gamma1) / 2
q = p
transition_matrix = rouwenhorst(N, p, q)

print(state_vector)
print(transition_matrix)