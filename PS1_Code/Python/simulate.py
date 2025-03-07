"""

simulate.py
-----------
This code simulates the model.

"""

#%% Imports from Python
import numpy as np
from numpy import cumsum,linspace,squeeze,where,zeros
from numpy.random import choice,rand,seed, normal
from numpy.linalg import matrix_power
from types import SimpleNamespace
from model import util

#%% Simulate the model.
def grow_economy(myClass):
    """
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    """

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulating the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions
    par = myClass.par
    sol = myClass.sol

    sigma = par.sigma  # CRRA risk aversion
    seed_sim = par.seed_sim  # Seed for simulation

    alpha = par.alpha  # Capital share
    delta = par.delta  # Depreciation rate

    # Tax and government spending processes
    rho_tau_k = par.rho_tau_k  # Persistence of capital tax rate
    sigma_tau_k = par.sigma_tau_k  # Std dev of capital tax shock
    rho_g = par.rho_g  # Persistence of government spending
    sigma_g = par.sigma_g  # Std dev of government spending shock

    # Initialize simulation variables
    T = par.T  # Number of time periods
    tau_k_sim = zeros(T * 2)  # Capital tax rate
    g_sim = zeros(T * 2)  # Government spending

    klen = par.klen  # Capital grid size
    Alen = par.Alen  # Productivity grid size
    kgrid = par.kgrid  # Capital today (state)
    Agrid = par.Agrid[0]  # Productivity today (state)
    pmat = par.pmat  # Transition matrix for productivity

    yout = sol.y  # Production function
    kpol = sol.k  # Capital policy function
    cpol = sol.c  # Consumption policy function
    ipol = sol.i  # Investment policy function

    # Initialize capital tax and government spending at steady state
    tau_k_sim[0] = par.tau_k_ss  # Initial capital tax rate
    g_sim[0] = par.g_ss  # Initial government spending

    # Initialize other economic variables
    Asim = zeros(T * 2)  # Productivity
    ysim = zeros(T * 2)  # Output
    ksim = zeros(T * 2)  # Capital stock
    csim = zeros(T * 2)  # Consumption
    isim = zeros(T * 2)  # Investment
    usim = zeros(T * 2)  # Utility

    # Set up random seed for reproducibility
    seed(seed_sim)

    # Compute stationary distribution for productivity
    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix

    # Initial states
    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index

    # Assign initial values
    Asim[0] = Agrid[A0_ind]  # Productivity
    ksim[0] = kgrid[k0_ind]  # Capital
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    usim[0] = util(csim[0], sigma)  # Utility

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = np.exp(par.rho * np.log(Asim[j - 1]) + normal(0, par.sigma_eps))

        # Capital tax follows AR(1) process
        tau_k_sim[j] = rho_tau_k * tau_k_sim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        g_sim[j] = rho_g * g_sim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        csim[j] = (1 - tau_k_sim[j]) * (ysim[j] - isim[j]) - g_sim[j]  # Consumption with taxes and government spending
        usim[j] = util(csim[j], sigma)  # Utility function

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]  # Update state for next period

    # Store results in simulation object
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_k_sim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = g_sim[T:2 * T + 1]  # Government spending