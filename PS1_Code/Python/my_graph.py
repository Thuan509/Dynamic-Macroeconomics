"""

my_graph.py
-----------
This code plots the value and policy functions.

"""

#%% Imports from Python
from matplotlib.pyplot import close,figure,plot,xlabel,ylabel,title,savefig,show
from numpy import linspace

#%% Plot the model functions and simulations.
def track_growth(myClass):
    '''
    
    This function plots the model functions and simulations.
    
    Input:
        myClass : Model class with parameters, grids, utility function, policy functions, and simulations.
        
    '''

    # Model parameters, policy and value functions, and simulations.
    par = myClass.par # Parameters.
    sol = myClass.sol # Policy functions.
    sim = myClass.sim # Simulations.

    # Production function.

    figure(1)
    plot(par.kgrid,sol.y)
    xlabel('$k_{t}$')
    ylabel('$y_{t}$') 
    title('Production Function')
    
    figname = myClass.par.figout+"\\ypol.png"
    print(figname)
    savefig(figname)

    # Plot capital policy function.
    
    figure(2)
    plot(par.kgrid,sol.k)
    xlabel('$k_{t}$')
    ylabel('$k_{t+1}$') 
    title('Capital Policy Function')
    
    figname = myClass.par.figout+"\\kpol.png"
    savefig(figname)

    # Plot consumption policy function.
    
    figure(3)
    plot(par.kgrid,sol.c)
    xlabel('$k_{t}$')
    ylabel('$c_{t}$') 
    title('Consumption Policy Function')
    
    figname = myClass.par.figout+"\\cpol.png"
    savefig(figname)

    # Plot investment policy function.
    
    figure(4)
    plot(par.kgrid,sol.i)
    xlabel('$k_{t}$')
    ylabel('$i_{t}$') 
    title('Investment Policy Function')
    
    figname = myClass.par.figout+"\\ipol.png"
    savefig(figname)

    # Plot labor supply policy function.
    
    figure(5)
    plot(par.kgrid,sol.n)
    xlabel('$k_{t}$')
    ylabel('$n_t$') 
    title('Labor Supply Policy Function')

    figname = myClass.par.figout+"\\npol.png"
    savefig(figname)
    
    # Plot value function.
    
    figure(6)
    plot(par.kgrid,sol.v)
    xlabel('$k_{t}$')
    ylabel('$V_t(k_t)$') 
    title('Value Function')

    figname = myClass.par.figout+"\\vfun.png"
    savefig(figname)
    
    # Plot simulated output.
    
    tgrid = linspace(1,par.T,par.T,dtype=int)

    figure(7)
    plot(tgrid,sim.ysim)
    xlabel('Time')
    ylabel('$y^{sim}_t$') 
    title('Simulated Output')

    figname = myClass.par.figout+"\\ysim.png"
    savefig(figname)
    
    # Plot simulated capital choice.
    
    figure(8)
    plot(tgrid,sim.ksim)
    xlabel('Time')
    ylabel('$k^{sim}_{t+1}$') 
    title('Simulated Capital Choice')

    figname = myClass.par.figout+"\\ksim.png"
    savefig(figname)
    
    # Plot simulated consumption.
    
    figure(9)
    plot(tgrid,sim.csim)
    xlabel('Time')
    ylabel('$c^{sim}_{t}$') 
    title('Simulated Consumption')

    figname = myClass.par.figout+"\\csim.png"
    savefig(figname)
    
    # Plot simulated investment.
    
    figure(10)
    plot(tgrid,sim.isim)
    xlabel('Time')
    ylabel('$i^{sim}_{t}$') 
    title('Simulated Investment')

    figname = myClass.par.figout+"\\isim.png"
    savefig(figname)
    
    # Plot simulated utility.
    
    figure(11)
    plot(tgrid,sim.usim)
    xlabel('Time')
    ylabel('$u^{sim}_t$') 
    title('Simulated Utility')

    figname = myClass.par.figout+"\\usim.png"
    savefig(figname)

    # Plot simulated productivity.
    
    figure(12)
    plot(tgrid,sim.Asim)
    xlabel('Time')
    ylabel('$A^{sim}_t$') 
    title('Simulated Productivity')

    figname = myClass.par.figout+"\\Asim.png"
    savefig(figname)

    # Plot simulated labor supply.
    
    figure(13)
    plot(tgrid,sim.nsim)
    xlabel('Time')
    ylabel('$n^{sim}_t$') 
    title('Simulated Labor Supply')

    figname = myClass.par.figout+"\\nsim.png"
    savefig(figname)

    #show()
    #close('all')