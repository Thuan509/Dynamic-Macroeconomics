"""

3b_tauchen.py
----------

This code discretizes a Markov process into N states using the method in Tauchen (EL,1986).

"""

#%% Imports.
import os
from scipy.stats import norm
from numpy import arange,count_nonzero,expand_dims,identity,linalg,linspace,nonzero,ones,prod,tile,zeros

#%% Tauchen's Method for AR(1).

def tauchen_ar1(mu,rho,sigma,N,m):
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
        y    : Grid for the AR(1) process.
        pmat : Transition probability matrix.
        
    """
    
    #%% Construct equally spaced grid.
    
    ar_mean = mu/(1-rho)
    ar_sd = sigma/((1-rho**2)**(1/2))
    
    y1 = ar_mean-(m*ar_sd)                                                                          # Smallest grid point.
    yn = ar_mean+(m*ar_sd)                                                                          # Largest grid point.
     
    y,d = linspace(y1,yn,N,endpoint=True,retstep=True)                                              # Equally spaced grid. Include endpoint and record stepsize.
        
    #%% Compute transition probability matrix from state j (row) to k (column).
    
    ymatk = tile(expand_dims(y,axis=0),(N,1))
    ymatj = mu+rho*ymatk.T
    
    pmat = norm.cdf(ymatk,loc=ymatj-(d/2),scale=sigma)-norm.cdf(ymatk,loc=ymatj+(d/2),scale=sigma)  # Transition probabilities to state 2, ..., N-1.
    pmat[:,0] = norm.cdf(y[0],loc=mu+rho*y-(d/2),scale=sigma)                                       # Transition probabilities to state 1.
    pmat[:,N-1] = 1-norm.cdf(y[N-1],loc=mu+rho*y+(d/2),scale=sigma)                                 # Transition probabilities to state N.
    
    #%% Output.
    
    y = expand_dims(y,axis=0)
    
    if count_nonzero(pmat.sum(axis=1)<0.99) > 0:
        raise Exception("Some columns of transition matrix don't sum to 1.") 
        
    return y,pmat

y, pmat = tauchen_ar1(0.5, 0.85, 1, 7, 3)

# Print the results
print("Grid points (y):")
print(y)
print("\nTransition probability matrix (pmat):")
print(pmat) 

## Define the output folder and file
output_folder = "Output"
output_file = os.path.join(output_folder, "tauchen_output.txt")

# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Save the results to a text file
with open(output_file, "w") as file:
    file.write("Grid Points (y):\n")
    file.write(str(y) + "\n\n")
    file.write("Transition Probability Matrix (pmat):\n")
    for row in pmat:
        file.write(" ".join(f"{val:.6f}" for val in row) + "\n")

print("Grid points (y):")
print(y)
print("\nTransition probability matrix (pmat):")
print(pmat)
print(f"Results saved to {output_file}")