"""

solve.py
--------
This code solves the model.

"""

#%% Imports from Python
import numpy as np
from numpy import argmax, expand_dims, inf, squeeze, tile, zeros, seterr
from numpy.linalg import norm
from types import SimpleNamespace
from model import util
import time

seterr(all='ignore')

#%% Solve the model using VFI.
def plan_allocations(myClass):
    '''
    Solves the stochastic growth model with inelastic labor (n=1),
    tax distortions, and government spending.
    
    Input:
        myClass : Model class with parameters, grids, and utility function.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Solving the Model by Value Function Iteration with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')
    
    # Namespace for optimal policy functions.
    setattr(myClass, 'sol', SimpleNamespace())
    sol = myClass.sol

    # Model parameters, grids and functions.
    par = myClass.par  # Parameters

    beta = par.beta  # Discount factor
    sigma = par.sigma  # CRRA
    alpha = par.alpha  # Capital share
    delta = par.delta  # Depreciation rate
    tau_k = par.tau_k_ss  # Capital tax rate
    tau_n = par.tau_n_ss  # Labor tax rate
    g_ss = par.g_ss  # Government spending

    klen = par.klen  # Grid size for capital
    kgrid = par.kgrid  # Capital grid

    Alen = par.Alen  # Grid size for A
    Agrid = par.Agrid[0]  # Productivity grid
    pmat = par.pmat  # Transition matrix for A

    kmat = tile(expand_dims(kgrid, axis=1), (1, Alen))  # Capital matrix
    Amat = tile(expand_dims(Agrid, axis=0), (klen, 1))  # Productivity matrix

    #util = par.util  # Utility function

    t0 = time.time()

    print('--------------------------------------Iterating on Bellman Eq.------------------------------------\n')

    # Value Function Iteration
    y0 = Amat * (kmat**alpha)  # Output with n=1
    i0 = delta * kmat  # In steady state, k = k' = k*
    c0 = (1 - tau_k) * (y0 - i0) - g_ss  # Consumption adjusted for taxes and government spending
    c0[c0 < 0.0] = 0.0
    v0 = np.where(c0 > 0, util(c0, g_ss, par.sigma), -np.inf) / (1.0 - par.beta)
  # Initial value function
    v0[c0 <= 0.0] = -inf  # Prevent negative consumption

    crit = 1e-6
    maxiter = 10000
    diff = 1
    iter = 0

    while (diff > crit) and (iter < maxiter):  # Iterate on the Bellman Equation until convergence
        v1 = zeros((klen, Alen))  # New value function
        k1 = zeros((klen, Alen))  # Capital policy function

        for p in range(0, klen):  # Loop over capital states
            for j in range(0, Alen):  # Loop over productivity states

                # Compute output with tax distortions
                y = Agrid[j] * (kgrid[p]**alpha)  # Output (n=1)
                i = kgrid - ((1 - delta) * kgrid[p])  # Investment: i = k' - (1 - delta) k
                c = (1 - tau_k) * (y - i) - g_ss  # Consumption adjusted for government spending
                c[c < 0.0] = 0.0

                # Solve the maximization problem
                ev = squeeze(v0 @ pmat[j, :].T)  # Expected value function
                vall = util(c, sigma) + beta * ev  # Bellman equation
                vall[c <= 0.0] = -inf  # Prevent negative consumption
                v1[p, j] = max(vall)  # Maximized value function
                k1[p, j] = kgrid[argmax(vall)]  # Optimal k'

        diff = norm(v1 - v0)  # Check convergence
        v0 = v1  # Update value function
        iter += 1  # Update counter

        # Print iteration progress
        if iter % 25 == 0:
            print('Iteration: ', iter, '.\n')

    t1 = time.time()
    print('Elapsed time is ', t1 - t0, ' seconds.')
    print('Converged in ', iter, ' iterations.')

    # Compute macro variables and policy functions
    sol.y = Amat * (kmat**alpha)  # Output (n=1)
    sol.k = k1  # Capital policy function
    sol.i = k1 - ((1.0 - delta) * kmat)  # Investment policy function
    sol.c = (1 - tau_k) * (sol.y - sol.i) - g_ss  # Adjusted consumption
    sol.c[sol.c < 0.0] = 0.0
    sol.v = v1  # Value function
    sol.v[sol.c <= 0.0] = -inf

#%% Utility function (Modified for inelastic labor)
def util(c, g, sigma):
    """ CRRA utility function (no labor choice since n=1) """
    # Government spending 
    ug = (g**(1 - sigma)) / (1 - sigma)
    # Consumption
    if sigma == 1.0:
        return np.log(c) + np.log(g)
    else:
        uc =  (c**(1 - sigma)) / (1 - sigma)
    # Total
    u = uc + ug
    return u
