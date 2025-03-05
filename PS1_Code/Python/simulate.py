"""

simulate.py
-----------
This code simulates the model.

"""

#%% Imports from Python
from numpy import cumsum,linspace,squeeze,where,zeros
from numpy.random import choice,rand,seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    
    This function simulates the stochastic growth model.
    
    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
        
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model')
    print('--------------------------------------------------------------------------------------------------\n')
    
    # Namespace for simulation.
    setattr(myClass,'sim',SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids and functions.
    
    par = myClass.par # Parameters.
    sol = myClass.sol # Policy functions.

    sigma = par.sigma # CRRA.
    gamma = par.gamma # Weight on leisure.
    nu = par.nu # Frisch Elasticity.
    util = par.util # Utility function.
    seed_sim = par.seed_sim # Seed for simulation.

    klen = par.klen # Capital grid size.
    Alen = par.Alen # Productivity grid size.
    kgrid = par.kgrid # Capital today (state).
    Agrid = par.Agrid[0] # Productivity today (state).
    pmat = par.pmat # Productivity today (state).
    
    yout = sol.y # Production function.
    kpol = sol.k # Policy function for capital.
    cpol = sol.c # Policy function for consumption.
    ipol = sol.i # Policy function for investment.
    npol = sol.n # Policy function for labor supply.

    T = par.T # Time periods.
    Asim = zeros(par.T*2) # Container for simulated productivity.
    ysim = zeros(par.T*2) # Container for simulated output.
    ksim = zeros(par.T*2) # Container for simulated capital stock.
    csim = zeros(par.T*2) # Container for simulated consumption.
    nsim = zeros(par.T*2) # Container for simulated labor supply.
    isim = zeros(par.T*2) # Container for simulated investment.
    usim = zeros(par.T*2) # Container for simulated utility.
            
    # Begin simulation.
    
    seed(seed_sim)

    pmat0 = matrix_power(pmat,1000)
    pmat0 = pmat0[0,:] # % Stationary distribution.
    cmat = cumsum(par.pmat,axis=1) # CDF matrix.

    A0_ind = choice(linspace(0,Alen,Alen,endpoint=False,dtype=int),1,p=pmat0) # Index for initial productivity.
    k0_ind = choice(linspace(0,klen,klen,endpoint=False,dtype=int),1) # Index for initial capital stock.

    Asim[0] = Agrid[A0_ind] # Productivity in period 1.
    ysim[0] = yout[k0_ind,A0_ind] # Output in period 1 given k0 and A0.
    csim[0] = cpol[k0_ind,A0_ind] # Consumption in period 1 given k0 and A0.
    ksim[0] = kpol[k0_ind,A0_ind] # Capital choice for period 2 given k0.
    nsim[0] = npol[k0_ind,A0_ind] # Labor supply in period 1 given k0 and A0.
    isim[0] = ipol[k0_ind,A0_ind] # Investment in period 1 given k0 and A0.
    usim[0] = util(csim[0],nsim[0],sigma,nu,gamma) # Utility in period 1 given k0 and A0.

    A1_ind = where(rand(1)<=squeeze(cmat[A0_ind,:])) # Draw productivity for next period.
    At_ind = A1_ind[0][0]

    # Simulate endogenous variables.

    for j in range(1,T*2): # Time loop.
        kt_ind = where(ksim[j-1]==kgrid); # Capital choice in the previous period is the state today. Find where the latter is on the grid.
        Asim[j] = Agrid[At_ind] # Productivity in period t.
        ysim[j] = yout[kt_ind,At_ind] # Output in period t.
        csim[j] = cpol[kt_ind,At_ind] # Consumption in period t.
        nsim[j] = npol[kt_ind,At_ind] # Labor supply for period t.
        ksim[j] = kpol[kt_ind,At_ind] # Capital stock for period t+1.
        isim[j] = ipol[kt_ind,At_ind] # Investment in period t.
        usim[j] = util(csim[j],nsim[j],sigma,nu,gamma) # Utility in period t.
        A1_ind = where(rand(1)<=squeeze(cmat[At_ind,:])) # Draw next state.
        At_ind = A1_ind[0][0] # State next period.

    # Burn the first half.
    sim.Asim = Asim[T:2*T+1] # Simulated productivity.
    sim.ysim = ysim[T:2*T+1] # Simulated output.
    sim.ksim = ksim[T:2*T+1] # Simulated capital choice.
    sim.csim = csim[T:2*T+1] # Simulated consumption.
    sim.nsim = nsim[T:2*T+1] # Simulated labor supply.
    sim.isim = isim[T:2*T+1] # Simulated investment.
    sim.usim = usim[T:2*T+1] # Simulated utility.