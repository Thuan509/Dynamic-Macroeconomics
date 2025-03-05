"""

solve.py
--------
This code solves the model.

"""

#%% Imports from Python
from numpy import argmax,expand_dims,inf,squeeze,tile,zeros,seterr
from numpy.linalg import norm
from scipy.optimize import fminbound
from types import SimpleNamespace
import time
seterr(all='ignore')

#%% Solve the model using VFI.
def plan_allocations(myClass):
    '''
    
    This function solves the stochastic growth model.
    
    Input:
        myClass : Model class with parameters, grids, and utility function.
        
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Solving the Model by Value Function Iteration')
    print('--------------------------------------------------------------------------------------------------\n')
    
    # Namespace for optimal policy funtions.
    setattr(myClass,'sol',SimpleNamespace())
    sol = myClass.sol

    # Model parameters, grids and functions.
    
    par = myClass.par # Parameters.

    beta = par.beta # Discount factor.
    sigma = par.sigma # CRRA.
    gamma = par.gamma # Weight on leisure.
    nu = par.nu # Frisch elasticity.

    alpha = par.alpha # Capital's share of income.
    delta = par.delta # Depreciation rate

    klen = par.klen # Grid size for k.
    kgrid = par.kgrid # Grid for k (state and choice).

    Alen = par.Alen # Grid size for A.
    Agrid = par.Agrid[0] # Grid for A.
    pmat = par.pmat # Grid for A.

    kmat = tile(expand_dims(kgrid,axis=1),(1,Alen)) # k for each value of A.
    Amat = tile(expand_dims(Agrid,axis=0),(klen,1)) # A for each value of k.

    util = par.util # Utility function.
    n0 = zeros((klen,klen,Alen))# Container for n.

    t0 = time.time()

    # Solve for labor choice.
    print('--------------------------------------Solving for Labor Supply------------------------------------\n')
    for h1 in range(0,klen): # Loop over k state.
        for h2 in range(0,klen): # Loop over k choice.
            for h3 in range(0,Alen): # Loop over A state.
                # Intratemporal condition.
                foc = lambda n: intra_foc(n,kgrid[h2],Agrid[h3],kgrid[h1],alpha,delta,sigma,nu,gamma)
                n0[h1,h2,h3] = fminbound(foc,0.0,1.0)
        
        # Print counter.
        if h1%25 == 0:
            print('Capital State: ',h1,'.\n')

    print('--------------------------------------Iterating on Bellman Eq.------------------------------------\n')

    # Value Function Iteration.
    y0 = Amat*(kmat**alpha)*(squeeze(n0[:,1,:])**(1.0-alpha)) # Given combinations of k and A and the value of n associated with the lowest possible k'.
    i0 = delta*kmat # In steady state, k=k'=k*.
    c0 = y0-i0 # Steady-state consumption.
    c0[c0<0.0] = 0.0
    v0 = util(c0,squeeze(n0[:,1,:]),sigma,nu,gamma)/(1.0-beta) # Guess of value function for each value of k.
    v0[c0<=0.0] = -inf # Set the value function to negative infinity number when c <= 0.
                    
    crit = 1e-6;
    maxiter = 10000;
    diff = 1;
    iter = 0;

    while (diff > crit) and (iter < maxiter): # Iterate on the Bellman Equation until convergence.
    
        v1 = zeros((klen,Alen)) # Container for V.
        k1 = zeros((klen,Alen)) # Container for k'.
        n1 = zeros((klen,Alen)) # Container for n.

        for p in range(0,klen): # Loop over the k-states.
            for j in range(0,Alen): # Loop over the A-states.

                # Macro variables.
                y = Agrid[j]*(kgrid[p]**alpha)*(squeeze(n0[p,:,j])**(1.0-alpha)) # Output in the "last period", given combinations of k and A and the value of n associated with the lowest possible k'.
                i = kgrid-((1-delta)*kgrid[p]) # Possible values for investment, i=k'-(1-delta)k, when choosing k' from kgrid and given k.
                c = y-i # Possible values for consumption, c = y-i, given y and i.
                c[c<0.0] = 0.0

                # Solve the maximization problem.
                ev = squeeze(v0@pmat[j,:].T); #  The next-period value function is the expected value function over each possible next-period A, conditional on the current state j.
                vall = util(c,n0[p,:,j],sigma,nu,gamma) + beta*ev # Compute the value function for each choice of k', given k.
                vall[c<=0.0] = -inf # Set the value function to negative infinity number when c <= 0.
                v1[p,j] = max(vall) # Maximize: vmax is the maximized value function; ind is where it is in the grid.
                k1[p,j] = kgrid[argmax(vall)] # Optimal k'.
                n1[p,j] = n0[p,argmax(vall),j] # Choice of n given k,k', and A.
        
        diff = norm(v1-v0) # Check convergence.
        v0 = v1; # Update guess.

        iter = iter + 1; # Update counter.
        
        # Print counter.
        if iter%25 == 0:
            print('Iteration: ',iter,'.\n')

    t1 = time.time()
    print('Elapsed time is ',t1-t0,' seconds.')
    print('Converged in ',iter,' iterations.')

    # Macro variables, value, and policy functions.
    sol.y = Amat*(kmat**alpha)*(n1**(1.0-alpha)) # Output.
    sol.k = k1 # Capital policy function.
    sol.n = n1 # Labor supply policy function.
    sol.i = k1-((1.0-delta)*kmat) # Investment policy function.
    sol.c = sol.y-sol.i # Consumption policy function.
    sol.c[sol.c<0.0] = 0.0
    sol.v = v1 # Value function.
    sol.v[sol.c<=0.0] = -inf

#%% Intra-temporal conditions for labor.
def intra_foc(n,kp,A,k,alpha,delta,sigma,nu,gamma):

    c = (A*(k**alpha)*(n**(1.0-alpha)))+((1.0-delta)*k-kp)
    mpl = A*(1.0-alpha)*((k/n)**alpha)

    # Leisure.
    un = -gamma*(1.0-n)**(1.0/nu)
    
    # Consumption.
    if sigma == 1.0:
        uc = 1.0/c # Log utility.
    else:
        uc = c**(-sigma) # CRRA utility.
    
    # Total.
    ucn = uc*mpl + un

    return ucn
