"""

simulate.py
-----------
This code simulates the model.

"""

#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
#%% Imports from Python
from numpy import cumsum, linspace, squeeze, where, zeros, exp, log
from numpy.random import choice, rand, normal, seed
from numpy.linalg import matrix_power
from types import SimpleNamespace

#%% Simulate the model.
def grow_economy(myClass):
    '''
    This function simulates the stochastic growth model with inelastic labor supply (n=1).

    Input:
        myClass : Model class with parameters, grids, utility function, and policy functions.
    '''

    print('\n--------------------------------------------------------------------------------------------------')
    print('Simulate the Model with Inelastic Labor Supply (n=1)')
    print('--------------------------------------------------------------------------------------------------\n')

    # Namespace for simulation.
    setattr(myClass, 'sim', SimpleNamespace())
    sim = myClass.sim

    # Model parameters, grids, and functions.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.

    sigma = par.sigma  # CRRA.
    util = par.util  # Utility function.
    seed_sim = par.seed_sim  # Seed for simulation.

    alpha = par.alpha  # Capital share in production.
    delta = par.delta  # Depreciation rate.
    
    rho_A = par.rho_A  # Persistence of productivity shock
    sigma_A = par.sigma_A  # Std deviation of productivity shock

    rho_tau_k = par.rho_tau_k  # Persistence of capital tax shock
    sigma_tau_k = par.sigma_tau_k  # Std deviation of capital tax shock

    rho_g = par.rho_g  # Persistence of government spending shock
    sigma_g = par.sigma_g  # Std deviation of government spending shock

    klen = par.klen  # Capital grid size.
    Alen = par.Alen  # Productivity grid size.
    kgrid = par.kgrid  # Capital today (state).
    Agrid = par.Agrid[0]  # Productivity today (state).
    pmat = par.pmat  # Productivity transition matrix.

    yout = sol.y  # Production function.
    kpol = sol.k  # Policy function for capital.
    cpol = sol.c  # Policy function for consumption.
    ipol = sol.i  # Policy function for investment.

    T = par.T  # Time periods.

    # Containers for simulation results
    Asim = zeros(par.T * 2)  # Productivity
    ysim = zeros(par.T * 2)  # Output
    ksim = zeros(par.T * 2)  # Capital stock
    csim = zeros(par.T * 2)  # Consumption
    isim = zeros(par.T * 2)  # Investment
    usim = zeros(par.T * 2)  # Utility
    tau_ksim = zeros(par.T * 2)  # Capital tax
    gsim = zeros(par.T * 2)  # Government spending

    # Initialize simulation
    seed(seed_sim)

    pmat0 = matrix_power(pmat, 1000)
    pmat0 = pmat0[0, :]  # Stationary distribution.
    cmat = cumsum(par.pmat, axis=1)  # CDF matrix.

    A0_ind = choice(linspace(0, Alen, Alen, endpoint=False, dtype=int), 1, p=pmat0)  # Initial productivity index.
    k0_ind = choice(linspace(0, klen, klen, endpoint=False, dtype=int), 1)  # Initial capital index.

    # Initial values
    Asim[0] = Agrid[A0_ind]  # Initial productivity
    ksim[0] = kgrid[k0_ind]  # Initial capital
    tau_ksim[0] = par.tau_k_ss  # Initial capital tax (steady state)
    gsim[0] = par.g_ss  # Initial government spending (steady state)

    # Compute initial output, investment, and consumption
    ysim[0] = yout[k0_ind, A0_ind]  # Output
    isim[0] = ipol[k0_ind, A0_ind]  # Investment
    csim[0] = cpol[k0_ind, A0_ind]  # Consumption
    usim[0] = util(csim[0], 1, sigma, 0, 0)  # Utility (n=1 is fixed)

    A1_ind = where(rand(1) <= squeeze(cmat[A0_ind, :]))  # Draw next productivity state
    At_ind = A1_ind[0][0]

    # Simulation loop
    for j in range(1, T * 2):
        # Productivity follows log-AR(1) process
        Asim[j] = exp(rho_A * log(Asim[j - 1]) + normal(0, sigma_A))

        # Capital tax follows AR(1) process
        tau_ksim[j] = rho_tau_k * tau_ksim[j - 1] + normal(0, sigma_tau_k)

        # Government spending follows AR(1) process
        gsim[j] = rho_g * gsim[j - 1] + normal(0, sigma_g)

        # Find closest capital index in the grid
        kt_ind = where(ksim[j - 1] == kgrid)[0][0]

        # Update economic variables
        ysim[j] = yout[kt_ind, At_ind]  # Output
        csim[j] = (1 - tau_ksim[j]) * (ysim[j] - isim[j]) - gsim[j]  # Consumption with tax and government spending
        isim[j] = ipol[kt_ind, At_ind]  # Investment
        ksim[j] = kpol[kt_ind, At_ind]  # Capital stock for next period
        usim[j] = util(csim[j], 1, sigma, 0, 0)  # Utility

        # Ensure non-negative consumption
        if csim[j] < 0:
            csim[j] = 0

        # Draw next productivity state
        A1_ind = where(rand(1) <= squeeze(cmat[At_ind, :]))
        At_ind = A1_ind[0][0]

    # Burn first half of the simulation
    sim.Asim = Asim[T:2 * T + 1]  # Productivity
    sim.ysim = ysim[T:2 * T + 1]  # Output
    sim.ksim = ksim[T:2 * T + 1]  # Capital
    sim.csim = csim[T:2 * T + 1]  # Consumption
    sim.isim = isim[T:2 * T + 1]  # Investment
    sim.usim = usim[T:2 * T + 1]  # Utility
    sim.tau_ksim = tau_ksim[T:2 * T + 1]  # Capital tax rate
    sim.gsim = gsim[T:2 * T + 1]  # Government spending
