"""

3c_markov.py
-------

This code call the functions needed to simulate and plot the Markov Chain for the AR(1) process.


"""

#%%
import os
import gc
gc.set_threshold(1,1,1)

import matplotlib.pyplot as plt
import numpy as np

# 
import functions
from functions import tauchen_ar1,simulate

#%%

"""

                            Discretize AR(1).

"""

mu    = 0.50000                                                                 # Intercept.
rho   = 0.85000                                                                 # Persistence.
sigma = 1.00000                                                                 # Std. dev. of error term.
N     = 7                                                                       # Number of grid points.
m     = 3                                                                       # For Tauchen (1986).

ar    = {'tau': tauchen_ar1} 
ar_y  = {}
ar_pi = {}

for i in ar:
    if i != "tau":
        ar_y[i], ar_pi[i] = ar[i](mu,rho,sigma,N)
    else:
        ar_y[i], ar_pi[i] = ar[i](mu,rho,sigma,N,m)
        
del mu,rho,sigma,N,m
gc.collect()

#%%

"""

      Simulate the AR(1)'s based on discretized processes.

"""

seed = 2025
T = 50

sim = {}

for i in ar_y.keys():
    np.random.seed(seed)
    sim["ar_"+i] = simulate(ar_y[i],ar_pi[i],T)
    
del seed
gc.collect()

#%%

"""

    Plot series.

"""
# Ensure the output folder exists
output_folder = "Output"
os.makedirs(output_folder, exist_ok=True)

time = range(0,T)

for i in sim:
    
    if sim[i].shape[0] == 1:
        
        plt.figure()
        plt.plot(time,np.squeeze(sim[i]))
        plt.xlabel('Time')
        plt.ylabel('Y')    
        plt.title("AR(1): "+i)
        plot_path = os.path.join(output_folder, f"{i}_plot.png")
        plt.savefig(plot_path)
        plt.show()

del i
gc.collect()

#%%

"""

                    Estimate the AR(1) series using simulated data.

"""

# Path for the AR(1) estimation output file
ar1_estimation_file = os.path.join(output_folder, "ar1_estimation.txt")

# Estimate the AR(1) series using simulated data
big1 = np.ones((T-1, 1))
bhat = {}

for i in sim:
    if sim[i].shape[0] == 1:
        y = np.expand_dims(sim[i][0, 1:T], axis=1)
        x = np.expand_dims(sim[i][0, 0:T-1], axis=1)
        bigX = np.concatenate((big1, x), axis=1)

        bhat[i] = np.linalg.inv(bigX.T @ bigX) @ bigX.T @ y

    else:
        y = sim[i][:, 1:T].T
        x = sim[i][:, 0:T-1].T
        bigX = np.concatenate((big1, x), axis=1)

        bhat[i] = (np.linalg.inv(bigX.T @ bigX) @ bigX.T @ y).T

# Save the results to a text file
with open(ar1_estimation_file, "w") as file:
    file.write("AR(1) Estimation Results:\n")
    for key, value in bhat.items():
        file.write(f"\nKey: {key}\n")
        file.write(f"Estimates:\n{value}\n")

print(f"AR(1) estimation results saved to {ar1_estimation_file}")

# Clean up variables
del i, T, ar, x, y, big1, bigX
gc.collect()
