"""

functions.py
----------

This code include the functions to discretizes an AR(1) process and simulate and plot the Markov Chain.

"""

#%% Imports.

from scipy.stats import norm
from numpy import arange,count_nonzero,expand_dims,identity,linalg,linspace,nonzero,ones,prod,tile,zeros,cumsum
from numpy.random import uniform

#%% Tauchen's Method for AR(1).

def tauchen_ar1(mu,rho,sigma,N,m):
    """
    
    This function discretizes an AR(1) process.
    
            y(t) = mu + rho*y(t-1) + eps(t), eps(t) ~ NID(0,sigma^2)
    
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


"""

                        Simulate Markov chain.

"""

def simulate(grid,pmat,T):
    """
    
    This code simulates a Markov chain given a grid of points from the 
    discretized process and the associated transition matrix.
    
    Input:
        grid : K x N Grid of discretized points.
        pmat : Transition probability matrix.
        seed : Set seed for rng.
        T    : Number of periods to simulate.
        
    Output:
        y : Simulated series.

    """
    
    #%% Initialize.
     
    cmat = cumsum(pmat,axis=1)                                                  # CDF matrix.
    
    if grid.shape[1]%2 == 0:                                                    # Initial state is the mean.
        state0 = int((grid.shape[1]/2)-1)
    else:
        state0 = int((grid.shape[1]-1)/2)
    
    y = zeros((grid.shape[0],T*2))
    
    for i in range(0,T*2):
        y[:,i] = grid[:,state0]
        state1 = cmat[state0,uniform()<=cmat[state0,:]]
        state0 = nonzero(cmat[state0,:]==state1[0])[0][0]
    
    y = y[:,T:T*2]
    
    #%% Output.
    
    return y