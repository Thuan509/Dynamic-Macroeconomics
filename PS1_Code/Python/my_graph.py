"""

my_graph.py
-----------
This code plots the value and policy functions.

"""

#%% Imports from Python
from matplotlib.pyplot import close, figure, plot, xlabel, ylabel, title, savefig, show
from numpy import linspace

#%% Plot the model functions and simulations.
def track_growth(myClass):
    '''
    This function plots the model functions and simulations.
    
    Input:
        myClass : Model class with parameters, grids, utility function, policy functions, and simulations.
    '''

    # Model parameters, policy and value functions, and simulations.
    par = myClass.par  # Parameters.
    sol = myClass.sol  # Policy functions.
    sim = myClass.sim  # Simulations.

    # Production function.
    figure(1)
    plot(par.kgrid, sol.y)
    xlabel('$k_{t}$')
    ylabel('$y_{t}$') 
    title('Production Function')
    savefig(par.figout + "\\ypol.png")

    # Capital policy function.
    figure(2)
    plot(par.kgrid, sol.k)
    xlabel('$k_{t}$')
    ylabel('$k_{t+1}$') 
    title('Capital Policy Function')
    savefig(par.figout + "\\kpol.png")

    # Consumption policy function.
    figure(3)
    plot(par.kgrid, sol.c)
    xlabel('$k_{t}$')
    ylabel('$c_{t}$') 
    title('Consumption Policy Function')
    savefig(par.figout + "\\cpol.png")

    # Investment policy function.
    figure(4)
    plot(par.kgrid, sol.i)
    xlabel('$k_{t}$')
    ylabel('$i_{t}$') 
    title('Investment Policy Function')
    savefig(par.figout + "\\ipol.png")

    # Value function.
    figure(5)
    plot(par.kgrid, sol.v)
    xlabel('$k_{t}$')
    ylabel('$V_t(k_t)$') 
    title('Value Function')
    savefig(par.figout + "\\vfun.png")

    # Time vector for simulations
    tgrid = linspace(1, par.T, par.T, dtype=int)

    # Simulated output.
    figure(6)
    plot(tgrid, sim.ysim)
    xlabel('Time')
    ylabel('$y^{sim}_t$') 
    title('Simulated Output')
    savefig(par.figout + "\\ysim.png")

    # Simulated capital choice.
    figure(7)
    plot(tgrid, sim.ksim)
    xlabel('Time')
    ylabel('$k^{sim}_{t+1}$') 
    title('Simulated Capital Choice')
    savefig(par.figout + "\\ksim.png")

    # Simulated consumption.
    figure(8)
    plot(tgrid, sim.csim)
    xlabel('Time')
    ylabel('$c^{sim}_{t}$') 
    title('Simulated Consumption')
    savefig(par.figout + "\\csim.png")

    # Simulated investment.
    figure(9)
    plot(tgrid, sim.isim)
    xlabel('Time')
    ylabel('$i^{sim}_{t}$') 
    title('Simulated Investment')
    savefig(par.figout + "\\isim.png")

    # Simulated utility.
    figure(10)
    plot(tgrid, sim.usim)
    xlabel('Time')
    ylabel('$u^{sim}_t$') 
    title('Simulated Utility')
    savefig(par.figout + "\\usim.png")

    # Simulated productivity.
    figure(11)
    plot(tgrid, sim.Asim)
    xlabel('Time')
    ylabel('$A^{sim}_t$') 
    title('Simulated Productivity')
    savefig(par.figout + "\\Asim.png")

    # Simulated capital tax rate.
    figure(12)
    plot(tgrid, sim.tau_ksim)
    xlabel('Time')
    ylabel('Capital Tax Rate ($\tau_k$)') 
    title('Simulated Capital Tax Rate')
    savefig(par.figout + "\\tau_ksim.png")

    # Simulated government spending.
    figure(13)
    plot(tgrid, sim.gsim)
    xlabel('Time')
    ylabel('Government Spending ($G_t$)') 
    title('Simulated Government Spending')
    savefig(par.figout + "\\gsim.png")

    #show()
    #close('all')